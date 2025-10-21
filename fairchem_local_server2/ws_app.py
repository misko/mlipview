from __future__ import annotations

import asyncio
import os
from collections import deque
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
        # Debug flags
        self._ws_debug = os.environ.get("WS_DEBUG", "0") in (
            "1",
            "true",
            "TRUE",
            "True",
        )
        self._log_calls = os.environ.get("WS_LOG_CALLS", "0") in (
            "1",
            "true",
            "TRUE",
            "True",
        )
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
        # JSON fallback removed: protobuf-only WS

        # Concurrency: producer (simulation loop) and consumer (message loop)
        stop_evt = asyncio.Event()

        # ---- minimal raw-frame queue you can peek into ----
        msg_buf: deque[pb.ClientAction] = deque()
        buf_has_data = asyncio.Event()
        reader_done = False

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
            stress: Optional[np.ndarray] = None,
            simulation_stopped: Optional[bool] = None,
        ):
            if self._ws_debug:
                try:
                    print(
                        f"[ws:tx] seq={seq} client_seq={client_seq} "
                        f"uic={user_interaction_count} step={sim_step} "
                        f"has={{pos:{positions is not None}, "
                        f"vel:{velocities is not None}, "
                        f"forces:{forces is not None}, "
                        f"cell:{cell is not None}}} "
                        f"energy={energy if energy is not None else 'na'} "
                        f"msg={message or ''}",
                        flush=True,
                    )
                except Exception:
                    pass
            # Protobuf binary only (ensure this runs regardless of debug flag)
            msg = pb.ServerResult()
            msg.seq = int(seq)
            if client_seq is not None:
                msg.client_seq = int(client_seq)
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
            if stress is not None:
                try:
                    if isinstance(stress, np.ndarray) and stress.shape == (3, 3):
                        sm = pb.Mat3()
                        sflat = [
                            float(stress[0, 0]),
                            float(stress[0, 1]),
                            float(stress[0, 2]),
                            float(stress[1, 0]),
                            float(stress[1, 1]),
                            float(stress[1, 2]),
                            float(stress[2, 0]),
                            float(stress[2, 1]),
                            float(stress[2, 2]),
                        ]
                        sm.m.extend(sflat)
                        msg.stress.CopyFrom(sm)
                except Exception:
                    pass
            if message:
                msg.message = message
            if simulation_stopped is True:
                try:
                    # Optional field: mark simulation stopped
                    msg.simulation_stopped = True
                except Exception:
                    pass
            # optional energy field
            try:
                if energy is not None:
                    msg.energy = float(energy)
            except Exception:
                pass
            await ws.send_bytes(msg.SerializeToString())

        async def sim_loop():
            nonlocal pending
            import time

            stall_notice_last = 0.0
            stall_notice_interval = 0.25  # seconds
            try:
                while not stop_evt.is_set():
                    if state.user_input_atomic_numbers is not None:
                        state.atomic_numbers = state.user_input_atomic_numbers
                        state.user_input_atomic_numbers = None
                    if state.user_input_positions is not None:
                        state.positions = state.user_input_positions
                        state.user_input_positions = None
                    if state.user_input_velocities is not None:
                        state.velocities = state.user_input_velocities
                        state.user_input_velocities = None
                    if state.user_input_cell is not None:
                        state.cell = state.user_input_cell
                        state.user_input_cell = None

                    if not state.running or state.atomic_numbers == []:
                        await asyncio.sleep(0.02)
                        continue

                    # Backpressure: don't exceed last_ack + 10
                    if state.server_seq - state.client_ack >= max_unacked:
                        # Send a periodic informational frame so clients
                        # know we're waiting on ACKs instead of stalling
                        now = time.time()
                        if now - stall_notice_last >= stall_notice_interval:
                            stall_notice_last = now
                            # Explicit debug for WAITING_FOR_ACK send
                            delta = state.server_seq - state.client_ack
                            if self._ws_debug:
                                print(
                                    (
                                        "[ws] WAITING_FOR_ACK send "
                                        f"server_seq={state.server_seq} "
                                        f"client_ack={state.client_ack} "
                                        f"delta={delta}"
                                    ),
                                    flush=True,
                                )
                            await _send_result_bytes(
                                seq=state.server_seq,
                                client_seq=state.client_seq,
                                user_interaction_count=(
                                    state.user_interaction_count or None
                                ),
                                sim_step=(state.sim_step or None),
                                positions=state.positions,
                                velocities=state.velocities,
                                forces=state.forces,
                                cell=None,
                                energy=None,
                                message="WAITING_FOR_ACK",
                            )
                        await asyncio.sleep(0.02)
                        continue

                    # Snapshot the user_interaction_count that is active at the
                    # beginning of this simulation step. This value will be
                    # echoed on the produced frame even if newer interactions
                    # arrive before the step completes.
                    uic_snapshot = state.user_interaction_count

                    worker = self._pool.any()
                    if state.sim_type == "md":
                        if self._log_calls:
                            print(
                                (
                                    f"[ws:sim][MD] steps=1 "
                                    f"T={state.params.temperature} "
                                    f"dt={state.params.timestep_fs} "
                                    f"friction={state.params.friction} "
                                    f"natoms={len(state.atomic_numbers)}"
                                ),
                                flush=True,
                            )
                        fut = worker.run_md.remote(
                            atomic_numbers=state.atomic_numbers,
                            positions=state.positions,
                            velocities=state.velocities,
                            cell=state.cell,
                            steps=1,
                            temperature=state.params.temperature,
                            timestep_fs=state.params.timestep_fs,
                            friction=state.params.friction,
                            calculator=state.params.calculator,
                        )
                    elif state.sim_type == "relax":
                        if self._log_calls:
                            print(
                                (
                                    f"[ws:sim][RELAX] steps=1 "
                                    f"fmax={state.params.fmax} "
                                    f"max_step={state.params.max_step} "
                                    f"calc={state.params.calculator} "
                                    f"natoms={len(state.atomic_numbers)}"
                                ),
                                flush=True,
                            )
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
                    energy_out = float(res.get("final_energy"))

                    state.server_seq += 1
                    # Advance server-side sim_step counter per produced frame
                    state.sim_step = int(state.sim_step or 0) + 1
                    await _send_result_bytes(
                        seq=state.server_seq,
                        client_seq=state.client_seq,
                        user_interaction_count=uic_snapshot,
                        sim_step=(state.sim_step or None),
                        positions=state.positions,
                        velocities=state.velocities,
                        forces=state.forces,
                        cell=state.cell,
                        energy=energy_out,
                    )
                    if state.server_seq % 10 == 0 and self._ws_debug:
                        print(
                            (
                                f"[ws] sent frame seq={state.server_seq} "
                                f"sim_step={state.sim_step}"
                            ),
                            flush=True,
                        )
            except Exception as e:
                print(f"[ws:sim_loop:error] {e}", flush=True)
                stop_evt.set()

        # --- recv_loop helpers: per-message-type handlers ---
        def _handle_ack_and_client_seq(msg) -> None:
            """Apply ACK/client_seq from incoming message with debug logs."""
            # client_seq
            state.client_seq = max(state.client_seq, int(getattr(msg, "seq", 0) or 0))
            # ack
            prev_ack = int(state.client_ack)
            if hasattr(msg, "HasField") and msg.HasField("ack"):
                state.client_ack = max(
                    state.client_ack,
                    int(getattr(msg, "ack", 0) or 0),
                )
            else:
                a = int(getattr(msg, "ack", 0) or 0)
                if a:
                    state.client_ack = max(state.client_ack, a)
            if self._ws_debug:
                if state.client_ack and state.client_ack > prev_ack:
                    delta = state.server_seq - state.client_ack

                    print(
                        (
                            "[ws] ACK recv "
                            f"client_ack={state.client_ack} "
                            f"server_seq={state.server_seq} "
                            f"delta={delta}"
                        ),
                        flush=True,
                    )
                elif state.client_ack:
                    print(
                        f"[ws:rx][ack] client_ack={state.client_ack}",
                        flush=True,
                    )

        def _extract_correlation_fields(msg) -> Optional[int]:
            """Return uic_in_msg; also update UIC and sim_step."""
            uic_in_msg: Optional[int] = None
            if hasattr(msg, "user_interaction_count"):
                u = int(getattr(msg, "user_interaction_count", 0) or 0)
                uic_in_msg = u
                if u > (state.user_interaction_count or 0):
                    state.user_interaction_count = u
            if hasattr(msg, "sim_step"):
                s = int(getattr(msg, "sim_step", 0) or 0)
                if not state.running and s > 0:
                    state.sim_step = s
            return uic_in_msg

        async def _handle_user_interaction(msg, uic_in_msg: Optional[int]) -> None:
            if hasattr(msg, "atomic_numbers") and len(msg.atomic_numbers) > 0:
                state.user_input_atomic_numbers = list(msg.atomic_numbers)
            if hasattr(msg, "positions") and len(msg.positions) > 0:
                state.user_input_positions = _np_from_vec3_list(msg.positions)
            if hasattr(msg, "velocities") and len(msg.velocities) > 0:
                state.user_input_velocities = _np_from_vec3_list(msg.velocities)
            if hasattr(msg, "cell"):
                state.user_input_cell = _mat3_to_np(msg.cell)
            n = (
                0
                if state.user_input_positions is None
                else int(state.user_input_positions.shape[0])
            )
            if self._ws_debug:
                print(
                    (
                        f"[ws] USER_INTERACTION recv natoms={n} "
                        "running="
                        f"{'true' if state.running else 'false'}"
                    ),
                    flush=True,
                )
            if state.running:
                # During run: drop forces so next produced frame recomputes
                state.forces = None
                if self._ws_debug:
                    print(
                        ("[ws:state] applied USER_INTERACTION during " "running sim"),
                        flush=True,
                    )
                return
            # Idle compute and send (positions omitted)
            print("[ws] USER_INTERACTION idle -> compute", flush=True)
            # copy over user input into the state
            if state.user_input_atomic_numbers is not None:
                state.atomic_numbers = list(state.user_input_atomic_numbers)
                state.user_input_atomic_numbers = None
                assert len(state.atomic_numbers) != 0
            if state.user_input_positions is not None:
                state.positions = np.array(state.user_input_positions, dtype=float)
                state.user_input_positions = None
            if state.user_input_velocities is not None:
                state.velocities = np.array(state.user_input_velocities, dtype=float)
                state.user_input_velocities = None
            if state.user_input_cell is not None:
                state.cell = np.array(state.user_input_cell, dtype=float)
                state.user_input_cell = None

            worker = self._pool.any()
            try:
                res = await asyncio.to_thread(
                    ray.get,
                    worker.run_simple.remote(
                        atomic_numbers=(state.atomic_numbers),
                        positions=state.positions,
                        cell=state.cell,
                        properties=["energy", "forces"],
                        calculator=state.params.calculator,
                    ),
                )
            except Exception as e:
                print(
                    (
                        "[ws] USER_INTERACTION compute failed with '"
                        f"{state.params.calculator}"
                        "' (no fallback): "
                        f"{e}"
                    ),
                    flush=True,
                )
                state.server_seq += 1
                await _send_result_bytes(
                    seq=state.server_seq,
                    client_seq=state.client_seq,
                    user_interaction_count=(int(uic_in_msg) if uic_in_msg else None),
                    sim_step=(state.sim_step or None),
                    positions=None,
                    velocities=None,
                    forces=None,
                    cell=None,
                    energy=None,
                    stress=None,
                    message=("COMPUTE_ERROR: " + str(e)),
                )
                return
            results = res.get("results", {})
            E = results.get("energy")
            F = results.get("forces")
            F_np = np.array(F, dtype=float) if F is not None else None
            S = results.get("stress") if isinstance(results, dict) else None
            S_np = None
            try:
                if S is not None:
                    S_np_arr = np.asarray(S, dtype=float)
                    if S_np_arr.shape == (3, 3):
                        S_np = S_np_arr
            except Exception:
                S_np = None
            state.server_seq += 1
            await _send_result_bytes(
                seq=state.server_seq,
                client_seq=state.client_seq,
                user_interaction_count=(int(uic_in_msg) if uic_in_msg else None),
                sim_step=(state.sim_step or None),
                positions=None,
                velocities=None,
                forces=F_np,
                cell=None,
                energy=(float(E) if E is not None else None),
                stress=S_np,
            )

        def _handle_start_simulation(msg) -> None:
            if not msg.HasField("simulation_type"):
                return
            state.sim_type = (
                "md" if (msg.simulation_type == pb.ClientAction.SimType.MD) else "relax"
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
            if self._ws_debug:
                print(
                    (
                        f"[ws] START_SIM type={state.sim_type} "
                        f"T={state.params.temperature} "
                        f"dt={state.params.timestep_fs} "
                        f"friction={state.params.friction}"
                    ),
                    flush=True,
                )
                print(
                    (
                        f"[ws:state] counters uic="
                        f"{state.user_interaction_count} "
                        f"sim_step={state.sim_step}"
                    ),
                    flush=True,
                )

        async def _handle_stop_simulation() -> None:
            state.running = False
            print("[ws] STOP_SIM", flush=True)
            if self._ws_debug:
                print("[ws:state] running=false", flush=True)
            try:
                state.server_seq += 1
                await _send_result_bytes(
                    seq=state.server_seq,
                    client_seq=state.client_seq,
                    user_interaction_count=(state.user_interaction_count or None),
                    sim_step=(state.sim_step or None),
                    positions=None,
                    velocities=None,
                    forces=None,
                    cell=None,
                    energy=None,
                    message="SIMULATION_STOPPED",
                    stress=None,
                    simulation_stopped=True,
                )
            except Exception:
                pass

        def _handle_ping(msg) -> None:
            if hasattr(msg, "ack"):
                a = int(getattr(msg, "ack", 0) or 0)
                if a:
                    state.client_ack = max(state.client_ack, a)
                    if self._ws_debug:
                        print(f"[ws] ACK {a}", flush=True)

        async def _reader():
            nonlocal reader_done
            try:
                while True:
                    try:
                        data = await ws.receive()
                    except WebSocketDisconnect:
                        break

                    evtype = data.get("type")
                    if evtype == "websocket.disconnect":
                        break
                    if evtype != "websocket.receive":
                        continue

                    b = data.get("bytes")
                    if b is None:
                        # Ignore text frames for this protobuf-only server
                        continue

                    # Parse now (protobuf-only server)
                    try:
                        m = pb.ClientAction()
                        m.ParseFromString(b)  # type: ignore[arg-type]
                    except Exception:
                        # bad frame; skip
                        continue

                    msg_buf.append(m)
                    buf_has_data.set()
            except WebSocketDisconnect:
                pass
            except RuntimeError:
                pass
            finally:
                reader_done = True
                buf_has_data.set()  # wake any waiters
                stop_evt.set()

        async def recv_loop():
            nonlocal last_ack
            try:
                while True:
                    if not msg_buf:
                        if reader_done:
                            break
                        try:
                            await asyncio.wait_for(buf_has_data.wait(), timeout=0.05)
                        except asyncio.TimeoutError:
                            continue
                        finally:
                            if not msg_buf:
                                buf_has_data.clear()
                        continue

                    # Optional PEEK without consuming:
                    # peek_msg: pb.ClientAction = msg_buf[0]

                    msg: pb.ClientAction = msg_buf.popleft()
                    if not msg_buf:
                        buf_has_data.clear()

                    # shared header handling
                    _handle_ack_and_client_seq(msg)
                    uic_in_msg = _extract_correlation_fields(msg)

                    # dispatch
                    mtype = getattr(msg, "type", None)
                    if (
                        hasattr(pb.ClientAction.Type, "USER_INTERACTION")
                        and mtype == pb.ClientAction.Type.USER_INTERACTION
                    ) or mtype == "user_interaction":
                        await _handle_user_interaction(msg, uic_in_msg)
                    elif (
                        hasattr(pb.ClientAction.Type, "START_SIMULATION")
                        and mtype == pb.ClientAction.Type.START_SIMULATION
                    ) or mtype == "start_simulation":
                        _handle_start_simulation(msg)
                    elif (
                        hasattr(pb.ClientAction.Type, "STOP_SIMULATION")
                        and mtype == pb.ClientAction.Type.STOP_SIMULATION
                    ) or mtype == "stop_simulation":
                        await _handle_stop_simulation()
                    elif (
                        hasattr(pb.ClientAction.Type, "PING")
                        and mtype == pb.ClientAction.Type.PING
                    ) or mtype == "ping":
                        _handle_ping(msg)
                    else:
                        continue
            finally:
                stop_evt.set()

        # Run concurrently until either loop exits
        await asyncio.gather(sim_loop(), _reader(), recv_loop())


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
