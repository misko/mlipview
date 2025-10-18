from __future__ import annotations

import asyncio
import json
from typing import Optional

import ray
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from ray import serve
from ray.serve.schema import LoggingConfig

from fairchem_local_server2.pb_adapter import (
    decode_action_from_bytes,
    encode_result_to_bytes,
)
from fairchem_local_server2.session_models import (
    ClientAction,
    ServerResult,
    SessionState,
)
from fairchem_local_server2.worker_pool import WorkerPool
from fairchem_local_server.model_runtime import (
    MODEL_NAME,
    TASK_NAME,
    UMA_DEPLOYMENT_NAME,
    _PredictDeploy,
    health_snapshot,
    install_predict_handle,
)

app = FastAPI(title="UMA Serve WS API", debug=False)


@serve.deployment(
    ray_actor_options={"num_gpus": 0},
    logging_config={"enable_access_log": False},
)
@serve.ingress(app)
class WSIngress:
    def __init__(self, predict_handle=None, pool_size: int = 8):
        # UMA handle is optional; if present, install for UMA calculator usage
        if predict_handle is not None:
            install_predict_handle(predict_handle)
        self._uma_handle = predict_handle
        self._pool = WorkerPool(size=pool_size, uma_handle=self._uma_handle)

    @app.get("/serve/health")
    def health(self):
        snap = health_snapshot()
        snap["status"] = "ok"
        return snap

    @app.get("/uma/stats")
    async def uma_stats(self):
        try:
            if self._uma_handle is None:
                return {"status": "error", "error": "UMA not installed"}
            stats = await self._uma_handle.get_stats.remote()  # type: ignore
            return {"status": "ok", "stats": stats}
        except Exception as e:
            return {"status": "error", "error": str(e)}

    @app.post("/uma/reset")
    async def uma_reset(self):
        try:
            if self._uma_handle is None:
                return {"status": "error", "error": "UMA not installed"}
            await self._uma_handle.reset_stats.remote()  # type: ignore
            return {"status": "ok"}
        except Exception as e:
            return {"status": "error", "error": str(e)}

    @app.websocket("/ws")
    async def ws(self, ws: WebSocket):
        await ws.accept()
        state = SessionState()
        max_unacked = 10
        last_ack = 0
        pending = 0

        # Concurrency: producer (simulation loop) and consumer (message loop)
        stop_evt = asyncio.Event()

        async def send_result(payload: ServerResult):
            # Try protobuf first, fallback to JSON
            b = encode_result_to_bytes(payload)
            if b is not None:
                await ws.send_bytes(b)
            else:
                await ws.send_text(payload.json())

        async def sim_loop():
            nonlocal pending
            try:
                while not stop_evt.is_set():
                    if not state.running or state.atomic_numbers == []:
                        await asyncio.sleep(0.02)
                        continue

                    # Backpressure: don't exceed last_ack + 10
                    if state.server_seq - state.client_ack >= max_unacked:
                        await asyncio.sleep(0.005)
                        continue

                    worker = self._pool.any()
                    if state.sim_type == "md":
                        v_in = state.velocities if state.velocities else None
                        fut = worker.run_md.remote(
                            atomic_numbers=state.atomic_numbers,
                            positions=state.positions,
                            velocities=v_in,
                            cell=state.cell,
                            steps=1,
                            temperature=state.params.temperature,
                            timestep_fs=state.params.timestep_fs,
                            friction=state.params.friction,
                            calculator=state.params.calculator,
                        )
                    elif state.sim_type == "relax":
                        fut = worker.run_relax.remote(
                            atomic_numbers=state.atomic_numbers,
                            positions=state.positions,
                            cell=state.cell,
                            steps=1,
                            fmax=state.params.fmax,
                            max_step=state.params.max_step,
                            calculator=state.params.calculator,
                        )
                    else:
                        await asyncio.sleep(0.01)
                        continue

                    # Await result (Ray ObjectRef) without blocking event loop
                    res = await asyncio.to_thread(ray.get, fut)
                    # Ray returns pydantic dict
                    state.positions = res["positions"]
                    state.velocities = (
                        res.get("velocities") or state.velocities
                    )
                    state.forces = res.get("forces") or state.forces

                    state.server_seq += 1
                    out = ServerResult(
                        seq=state.server_seq,
                        client_seq=state.client_seq,
                        positions=state.positions,
                        velocities=state.velocities,
                        forces=state.forces,
                        cell=state.cell,
                    )
                    await send_result(out)
            except Exception:
                # On any error, end the session loop
                stop_evt.set()

        async def recv_loop():
            nonlocal last_ack
            try:
                while True:
                    # Accept binary (protobuf) or text (JSON)
                    msg = None
                    try:
                        data = await ws.receive()
                    except WebSocketDisconnect:
                        break
                    # data is a dict with keys like 'type' and one of
                    # 'text' or 'bytes'
                    evtype = data.get("type")
                    if evtype == "websocket.disconnect":
                        break
                    if evtype == "websocket.receive":
                        if "bytes" in data and data["bytes"] is not None:
                            m = decode_action_from_bytes(data["bytes"] or b"")
                            if m is not None:
                                msg = m
                        elif "text" in data and data["text"] is not None:
                            raw = data["text"]
                            try:
                                msg = ClientAction.parse_raw(raw)
                            except Exception:
                                # accept plain JSON dicts or ping strings
                                if raw == "ping":
                                    await ws.send_text("pong")
                                    continue
                                try:
                                    msg = ClientAction.parse_obj(
                                        json.loads(raw)
                                    )
                                except Exception as e:
                                    await ws.send_text(
                                        json.dumps(
                                            {
                                                "error": f"bad message: {e}"
                                            }
                                        )
                                    )
                                    continue
                    if msg is None:
                        continue

                    # ack handling
                    state.client_seq = max(state.client_seq, int(msg.seq))
                    if msg.ack is not None:
                        state.client_ack = max(state.client_ack, int(msg.ack))

                    if msg.type == "init_system":
                        if msg.atomic_numbers is None or msg.positions is None:
                            await ws.send_text(
                                json.dumps(
                                    {"error": "init_system missing fields"}
                                )
                            )
                            continue
                        state.atomic_numbers = list(msg.atomic_numbers)
                        state.positions = list(msg.positions)
                        state.velocities = (
                            list(msg.velocities) if msg.velocities else []
                        )
                        state.cell = msg.cell
                        state.forces = []
                        await send_result(
                            ServerResult(
                                seq=state.server_seq,
                                client_seq=state.client_seq,
                                positions=state.positions,
                                velocities=state.velocities or [],
                                forces=state.forces or [],
                                cell=state.cell,
                                message="initialized",
                            )
                        )
                    elif msg.type == "update_positions":
                        if msg.positions is None:
                            continue
                        state.positions = list(msg.positions)
                        # Reset v/f on any changed atoms
                        state.velocities = []
                        state.forces = []
                    elif msg.type == "start_simulation":
                        if msg.simulation_type is None:
                            await ws.send_text(
                                json.dumps(
                                    {"error": "missing simulation_type"}
                                )
                            )
                            continue
                        state.sim_type = msg.simulation_type
                        if msg.simulation_params is not None:
                            state.params = msg.simulation_params
                        state.running = True
                    elif msg.type == "stop_simulation":
                        state.running = False
                    elif msg.type == "ping":
                        # Treat 'ping' as a heartbeat or ACK carrier. If it's
                        # carrying an ACK (msg.ack present), do not respond to
                        # avoid spurious non-protobuf frames. Otherwise, reply.
                        if msg.ack is None:
                            await ws.send_text("pong")
            except WebSocketDisconnect:
                pass
            except RuntimeError:
                # e.g., "Cannot call receive once a disconnect message..."
                pass
            finally:
                stop_evt.set()

        # Run concurrently until either loop exits
        await asyncio.gather(sim_loop(), recv_loop())


def _detect_default_ngpus() -> int:
    try:
        import torch  # type: ignore
        if torch.cuda.is_available():
            # one replica per visible GPU
            return max(1, torch.cuda.device_count())
        # No CUDA -> default to 0 UMA replicas
        return 0
    except Exception:
        return 0


def _detect_default_ncpus() -> int:
    try:
        import multiprocessing as mp

        return max(1, mp.cpu_count())
    except Exception:
        return 8


def deploy(
    ngpus: Optional[int] = None,
    ncpus: Optional[int] = None,
    nhttp: Optional[int] = None,
):
    """Deploy UMA GPU replicas and WS ingress.

    - UMA predictor replicas: one per GPU
    - ASE worker pool size: ncpus (CPU-bound)
    """
    replica_count = (
        int(ngpus) if ngpus is not None else _detect_default_ngpus()
    )
    pool_size = int(ncpus) if ncpus is not None else _detect_default_ncpus()
    ingress_replicas = max(1, int(nhttp) if nhttp is not None else 1)
    dag = None
    if replica_count > 0:
        uma = _PredictDeploy.options(
            name=UMA_DEPLOYMENT_NAME,
            num_replicas=int(replica_count),
            ray_actor_options={"num_gpus": 1},
        ).bind(MODEL_NAME, TASK_NAME)
        dag = WSIngress.options(num_replicas=ingress_replicas).bind(
            uma, pool_size
        )
    else:
        # LJ-only or client-provided calculator flows
        dag = WSIngress.options(num_replicas=ingress_replicas).bind(
            None, pool_size
        )
    serve.run(
        dag,
        name="ws_app",
        route_prefix="/",
        logging_config=LoggingConfig(
            enable_access_log=False,
            log_level="ERROR",
        ),
    )
    return dag


if __name__ == "__main__":
    import argparse
    import time

    parser = argparse.ArgumentParser(description="Start UMA WS Serve app")
    parser.add_argument("--ngpus", type=int, default=None)
    parser.add_argument("--ncpus", type=int, default=None)
    args = parser.parse_args()

    deploy(args.ngpus, args.ncpus, None)
    while True:
        time.sleep(60)
