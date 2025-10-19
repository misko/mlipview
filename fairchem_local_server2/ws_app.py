from __future__ import annotations

import asyncio
from typing import Optional

import numpy as np
import ray
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from ray import serve
from ray.serve.schema import LoggingConfig

from fairchem_local_server2 import session_pb2 as pb
from fairchem_local_server2.session_models import SessionState, SimulationParams
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
        print(
            (
                f"[ws:init] pool_size={pool_size} "
                f"uma_handle={'yes' if predict_handle is not None else 'no'}"
            ),
            flush=True,
        )

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
        try:
            print("[ws] accepted connection", flush=True)
        except Exception:
            pass
        state = SessionState()
        max_unacked = 10
        last_ack = 0
        pending = 0
        # switch to JSON I/O if client sends JSON text frames
        json_mode = False

        # Concurrency: producer (simulation loop) and consumer (message loop)
        stop_evt = asyncio.Event()

        def _np_from_vec3_list(vecs):
            if not vecs:
                return None
            out = []
            for v in vecs:
                try:
                    out.append([float(v.v[0]), float(v.v[1]), float(v.v[2])])
                except Exception:
                    # assume plain [x,y,z]
                    out.append([float(v[0]), float(v[1]), float(v[2])])
            return np.array(out, dtype=float)

        def _mat3_to_np(m: Optional[pb.Mat3]):
            if m is None or len(m.m) != 9:
                return None
            return np.array(
                [
                    [m.m[0], m.m[1], m.m[2]],
                    [m.m[3], m.m[4], m.m[5]],
                    [m.m[6], m.m[7], m.m[8]],
                ],
                dtype=float,
            )

        def _fill_vec3_repeated(dst, arr: Optional[np.ndarray]):
            if arr is None:
                return
            for row in np.asarray(arr, dtype=float):
                v = pb.Vec3()
                v.v.extend([float(row[0]), float(row[1]), float(row[2])])
                dst.append(v)

        async def _send_result_bytes(
            *,
            seq: int,
            client_seq: Optional[int],
            user_interaction_count: Optional[int],
            sim_step: Optional[int],
            positions: Optional[np.ndarray],
            velocities: Optional[np.ndarray],
            forces: Optional[np.ndarray],
            cell: Optional[np.ndarray],
            message: Optional[str] = None,
            energy: Optional[float] = None,
        ):
            if json_mode:
                # JSON payload for browser clients without protobuf
                out = {
                    "seq": int(seq),
                    "client_seq": int(client_seq or 0),
                }
                if user_interaction_count is not None:
                    out["userInteractionCount"] = int(user_interaction_count)
                if sim_step is not None:
                    out["simStep"] = int(sim_step)
                if positions is not None:
                    out["positions"] = np.asarray(positions).tolist()
                if velocities is not None:
                    out["velocities"] = np.asarray(velocities).tolist()
                if forces is not None:
                    out["forces"] = np.asarray(forces).tolist()
                if cell is not None:
                    out["cell"] = np.asarray(cell).tolist()
                if energy is not None:
                    out["energy"] = float(energy)
                if message:
                    out["message"] = message
                await ws.send_json(out)
            else:
                # Protobuf binary
                msg = pb.ServerResult()
                msg.seq = int(seq)
                if client_seq is not None:
                    msg.client_seq = int(client_seq)
                # Optional correlation fields
                try:
                    if user_interaction_count is not None:
                        msg.user_interaction_count = int(user_interaction_count)
                    if sim_step is not None:
                        msg.sim_step = int(sim_step)
                except Exception:
                    pass
                _fill_vec3_repeated(msg.positions, positions)
                _fill_vec3_repeated(msg.forces, forces)
                _fill_vec3_repeated(msg.velocities, velocities)
                if cell is not None:
                    m = pb.Mat3()
                    flat = [
                        float(cell[0, 0]),
                        float(cell[0, 1]),
                        float(cell[0, 2]),
                        float(cell[1, 0]),
                        float(cell[1, 1]),
                        float(cell[1, 2]),
                        float(cell[2, 0]),
                        float(cell[2, 1]),
                        float(cell[2, 2]),
                    ]
                    m.m.extend(flat)
                    msg.cell.CopyFrom(m)
                if message:
                    msg.message = message
                # optional energy field
                try:
                    if energy is not None:
                        msg.energy = float(energy)
                except Exception:
                    pass
                await ws.send_bytes(msg.SerializeToString())

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
                        v_in = (
                            state.velocities if state.velocities is not None else None
                        )
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
                    state.positions = np.array(res["positions"], dtype=float)
                    v_out = res.get("velocities")
                    state.velocities = (
                        np.array(v_out, dtype=float)
                        if v_out is not None
                        else state.velocities
                    )
                    f_out = res.get("forces")
                    state.forces = (
                        np.array(f_out, dtype=float)
                        if f_out is not None
                        else state.forces
                    )
                    # Prefer per-step final energy if provided by the worker
                    try:
                        _fe = res.get("final_energy")
                        energy_out = float(_fe) if _fe is not None else None
                    except Exception:
                        energy_out = None

                    state.server_seq += 1
                    # Advance server-side sim_step counter per produced frame
                    try:
                        state.sim_step = int(state.sim_step or 0) + 1
                    except Exception:
                        state.sim_step = 1
                    await _send_result_bytes(
                        seq=state.server_seq,
                        client_seq=state.client_seq,
                        user_interaction_count=(state.user_interaction_count or None),
                        sim_step=(state.sim_step or None),
                        positions=state.positions,
                        velocities=state.velocities,
                        forces=state.forces,
                        cell=state.cell,
                        energy=energy_out,
                    )
                    try:
                        if state.server_seq % 10 == 0:
                            print(
                                (
                                    f"[ws] sent frame seq={state.server_seq} "
                                    f"sim_step={state.sim_step}"
                                ),
                                flush=True,
                            )
                    except Exception:
                        pass
            except Exception:
                # On any error, end the session loop
                stop_evt.set()

        async def recv_loop():
            nonlocal last_ack
            nonlocal json_mode
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
                            try:
                                m = pb.ClientAction()
                                m.ParseFromString(
                                    data["bytes"]
                                )  # type: ignore[arg-type]
                                msg = m
                            except Exception:
                                msg = None
                        elif "text" in data and data["text"] is not None:
                            # Try JSON mode: parse a JSON object and map
                            # to fields
                            import json

                            try:
                                obj = json.loads(data["text"])
                                # Toggle JSON mode so sender uses JSON frames
                                json_mode = True
                                # Minimal mapping for required fields
                                _type = str(obj.get("type") or "").lower()
                                if _type == "ping":
                                    # ack passthrough
                                    a = int(obj.get("ack") or 0)
                                    if a:
                                        state.client_ack = max(state.client_ack, a)
                                    # don't generate a reply in JSON mode
                                    continue

                                # Construct a pseudo message namespace
                                # using dicts
                                class _Msg:
                                    pass

                                m = _Msg()
                                setattr(m, "seq", int(obj.get("seq") or 0))
                                setattr(m, "type", _type)
                                setattr(
                                    m,
                                    "atomic_numbers",
                                    list(
                                        map(
                                            int,
                                            obj.get("atomic_numbers") or [],
                                        )
                                    ),
                                )
                                setattr(m, "positions", obj.get("positions") or [])
                                _vel = obj.get("velocities") or []
                                setattr(m, "velocities", _vel)
                                setattr(m, "cell", obj.get("cell"))
                                setattr(m, "ack", int(obj.get("ack") or 0))
                                sp = obj.get("simulation_params") or {}
                                setattr(m, "simulation_params", sp)
                                setattr(
                                    m,
                                    "simulation_type",
                                    obj.get("simulation_type"),
                                )
                                msg = m
                            except Exception:
                                # Not JSON: treat as ping
                                if data["text"] == "ping":
                                    await ws.send_text("pong")
                                continue
                    if msg is None:
                        continue

                    # ack handling
                    state.client_seq = max(
                        state.client_seq, int(getattr(msg, "seq", 0) or 0)
                    )
                    # ack handling
                    try:
                        if hasattr(msg, "HasField") and msg.HasField("ack"):
                            state.client_ack = max(state.client_ack, int(msg.ack))
                        else:
                            a = int(getattr(msg, "ack", 0) or 0)
                            if a:
                                state.client_ack = max(state.client_ack, a)
                    except Exception:
                        pass
                    # Correlation fields (optional)
                    try:
                        # Only access attributes if present in schema
                        if hasattr(msg, "user_interaction_count"):
                            u = int(getattr(msg, "user_interaction_count", 0) or 0)
                            if u > (state.user_interaction_count or 0):
                                state.user_interaction_count = u
                        if hasattr(msg, "sim_step"):
                            s = int(getattr(msg, "sim_step", 0) or 0)
                            # Keep last seen client step before run starts
                            if not state.running and s > 0:
                                state.sim_step = s
                    except Exception:
                        pass

                    if (
                        hasattr(pb.ClientAction.Type, "INIT_SYSTEM")
                        and getattr(msg, "type", None)
                        == pb.ClientAction.Type.INIT_SYSTEM
                    ) or getattr(msg, "type", None) in ("init_system", "init"):
                        if len(msg.atomic_numbers) == 0 or len(msg.positions) == 0:
                            # Cannot initialize without atoms and positions;
                            # ignore
                            continue
                        state.atomic_numbers = list(msg.atomic_numbers)
                        state.positions = _np_from_vec3_list(msg.positions)
                        v_in = _np_from_vec3_list(msg.velocities)
                        state.velocities = v_in
                        state.cell = _mat3_to_np(
                            msg.cell if msg.HasField("cell") else None
                        )
                        state.forces = None
                        try:
                            cs = "on" if state.cell is not None else "off"
                            print(
                                (
                                    "[ws] INIT_SYSTEM "
                                    f"natoms={len(state.atomic_numbers)} "
                                    f"cell={cs}"
                                ),
                                flush=True,
                            )
                        except Exception:
                            pass
                        await _send_result_bytes(
                            seq=state.server_seq,
                            client_seq=state.client_seq,
                            user_interaction_count=(
                                state.user_interaction_count or None
                            ),
                            sim_step=(state.sim_step or None),
                            positions=state.positions,
                            velocities=(
                                state.velocities
                                if state.velocities is not None
                                else np.zeros(
                                    (len(state.atomic_numbers), 3),
                                    dtype=float,
                                )
                            ),
                            forces=np.zeros(
                                (len(state.atomic_numbers), 3),
                                dtype=float,
                            ),
                            cell=state.cell,
                            message="initialized",
                        )
                    elif (
                        hasattr(pb.ClientAction.Type, "SIMPLE_CALCULATE")
                        and getattr(msg, "type", None)
                        == pb.ClientAction.Type.SIMPLE_CALCULATE
                    ) or getattr(msg, "type", None) == "simple_calculate":
                        # one-shot energy/forces; use current state
                        if state.atomic_numbers == [] or state.positions is None:
                            continue
                        try:
                            worker = self._pool.any()
                            res = await asyncio.to_thread(
                                ray.get,
                                worker.run_simple.remote(
                                    atomic_numbers=state.atomic_numbers,
                                    positions=state.positions,
                                    cell=state.cell,
                                    properties=["energy", "forces"],
                                    calculator=state.params.calculator,
                                ),
                            )
                            results = res.get("results", {})
                            E = results.get("energy")
                            F = results.get("forces")
                            if F is not None:
                                F_np = np.array(F, dtype=float)
                            else:
                                F_np = None
                            # server_seq unchanged for one-shot; echo fields
                            await _send_result_bytes(
                                seq=state.server_seq,
                                client_seq=state.client_seq,
                                user_interaction_count=(
                                    state.user_interaction_count or None
                                ),
                                sim_step=(state.sim_step or None),
                                positions=None,
                                velocities=None,
                                forces=F_np,
                                cell=None,
                                energy=(float(E) if E is not None else None),
                                message="simple_calculate",
                            )
                        except Exception:
                            # swallow to keep session alive
                            pass
                    elif (
                        hasattr(pb.ClientAction.Type, "UPDATE_POSITIONS")
                        and getattr(msg, "type", None)
                        == pb.ClientAction.Type.UPDATE_POSITIONS
                    ) or getattr(msg, "type", None) == "update_positions":
                        if len(msg.positions) == 0:
                            continue
                        state.positions = _np_from_vec3_list(msg.positions)
                        # Reset v/f on any changed atoms
                        state.velocities = None
                        state.forces = None
                        try:
                            print("[ws] UPDATE_POSITIONS", flush=True)
                        except Exception:
                            pass
                    elif (
                        hasattr(pb.ClientAction.Type, "START_SIMULATION")
                        and getattr(msg, "type", None)
                        == pb.ClientAction.Type.START_SIMULATION
                    ) or getattr(msg, "type", None) == "start_simulation":
                        if not msg.HasField("simulation_type"):
                            continue
                        state.sim_type = (
                            "md"
                            if (msg.simulation_type == pb.ClientAction.SimType.MD)
                            else "relax"
                        )
                        if msg.HasField("simulation_params"):
                            sp = msg.simulation_params
                            state.params = SimulationParams(
                                calculator=sp.calculator or "uma",
                                temperature=float(sp.temperature),
                                timestep_fs=float(sp.timestep_fs),
                                friction=float(sp.friction),
                                fmax=float(sp.fmax),
                                max_step=float(sp.max_step),
                                optimizer=(sp.optimizer or "bfgs"),
                            )
                        state.running = True
                        try:
                            print(
                                (
                                    f"[ws] START_SIM type={state.sim_type} "
                                    f"T={state.params.temperature} "
                                    f"dt={state.params.timestep_fs} "
                                    f"friction={state.params.friction}"
                                ),
                                flush=True,
                            )
                        except Exception:
                            pass
                    elif (
                        hasattr(pb.ClientAction.Type, "STOP_SIMULATION")
                        and getattr(msg, "type", None)
                        == pb.ClientAction.Type.STOP_SIMULATION
                    ) or getattr(msg, "type", None) == "stop_simulation":
                        state.running = False
                        try:
                            print("[ws] STOP_SIM", flush=True)
                        except Exception:
                            pass
                    elif (
                        hasattr(pb.ClientAction.Type, "PING")
                        and (getattr(msg, "type", None) == pb.ClientAction.Type.PING)
                    ) or getattr(msg, "type", None) == "ping":
                        # In protobuf-only mode, we don't echo unless needed.
                        try:
                            if hasattr(msg, "ack"):
                                a = int(getattr(msg, "ack", 0) or 0)
                                if a:
                                    state.client_ack = max(state.client_ack, a)
                                    print(f"[ws] ACK {a}", flush=True)
                        except Exception:
                            pass
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
    http_port: Optional[int] = None,
):
    """Deploy UMA GPU replicas and WS ingress.

    - UMA predictor replicas: one per GPU
    - ASE worker pool size: ncpus (CPU-bound)
    """
    replica_count = int(ngpus) if ngpus is not None else _detect_default_ngpus()
    pool_size = int(ncpus) if ncpus is not None else _detect_default_ncpus()
    ingress_replicas = max(1, int(nhttp) if nhttp is not None else 1)
    dag = None
    # If a specific HTTP port is requested, start Serve with that port
    if http_port is not None:
        try:
            serve.start(
                detached=False,
                http_options={
                    "host": "0.0.0.0",
                    "port": int(http_port),
                },
            )
        except Exception:
            # If Serve is already started, ignore error
            pass
    if replica_count > 0:
        uma = _PredictDeploy.options(
            name=UMA_DEPLOYMENT_NAME,
            num_replicas=int(replica_count),
            ray_actor_options={"num_gpus": 1},
        ).bind(MODEL_NAME, TASK_NAME)
        dag = WSIngress.options(num_replicas=ingress_replicas).bind(uma, pool_size)
    else:
        # LJ-only or client-provided calculator flows
        dag = WSIngress.options(num_replicas=ingress_replicas).bind(None, pool_size)
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
