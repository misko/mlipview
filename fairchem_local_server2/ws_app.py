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

# ---- helpers: protobuf <-> numpy -------------------------------------------


def _type_name(m) -> str:
    """Return upper-case enum name if available, else stringified type."""
    t = getattr(m, "type", None)
    enum_type = getattr(getattr(m, "__class__", None), "Type", None)
    if enum_type is not None and isinstance(t, int):
        try:
            return enum_type.Name(t)
        except Exception:
            pass
    return str(t or "").upper()


def _is_user_interaction(m) -> bool:
    return _type_name(m) == "USER_INTERACTION"


def _pb_vec3s_to_np(vecs) -> np.ndarray:
    """Vec3 { repeated double v } -> (N,3) float64 array with validation."""
    if not vecs:
        return np.empty((0, 3), dtype=np.float64)
    # Validate and copy
    out = np.empty((len(vecs), 3), dtype=np.float64)
    for i, v in enumerate(vecs):
        # v.v is a RepeatedScalarContainer; ensure length 3
        if len(v.v) != 3:
            raise ValueError(f"Vec3.v must have length 3, got {len(v.v)} at index {i}")
        out[i, 0] = float(v.v[0])
        out[i, 1] = float(v.v[1])
        out[i, 2] = float(v.v[2])
    return out


def _np_to_pb_vec3s(arr: Optional[np.ndarray], dest_repeated) -> None:
    """Fill a repeated Vec3 (with 'v' of length 3) from a (N,3) array. Clears dest first."""
    del dest_repeated[:]
    if arr is None:
        return
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim == 1:
        if a.size % 3 != 0:
            raise ValueError("Vec3 array length must be divisible by 3")
        a = a.reshape(-1, 3)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError("Vec3 array must have shape (N, 3)")
    for x, y, z in a:
        v = pb.Vec3()
        v.v.extend((float(x), float(y), float(z)))
        dest_repeated.append(v)


def _np_to_mat3(arr: Optional[np.ndarray]) -> Optional["pb.Mat3"]:
    if arr is None:
        return None
    arr = np.asarray(arr, dtype=np.float64)
    if arr.shape != (3, 3):
        raise ValueError("cell must be shape (3,3)")
    m = pb.Mat3()
    m.m.extend([float(x) for x in arr.reshape(-1)])
    return m


def _mat3_to_np(m: "pb.Mat3") -> np.ndarray:
    """Mat3 { repeated double m (len 9 row-major) } -> (3,3) float64 array."""
    a = np.fromiter(getattr(m, "m", []), dtype=np.float64, count=9)
    if a.size != 9:
        raise ValueError(f"Mat3.m must contain 9 elements, got {a.size}")
    return a.reshape(3, 3)


# ---- message coalescer ------------------------------------------------------


def _coalesce_user_interactions(msg_buf: "deque[pb.ClientAction]") -> "pb.ClientAction":
    """Consume consecutive USER_INTERACTION messages and merge with 'last write wins'."""
    uis: list[pb.ClientAction] = []
    while msg_buf and _is_user_interaction(msg_buf[0]):
        uis.append(msg_buf.popleft())

    merged = pb.ClientAction()
    if hasattr(pb.ClientAction, "Type") and hasattr(
        pb.ClientAction.Type, "USER_INTERACTION"
    ):
        merged.type = pb.ClientAction.Type.USER_INTERACTION
    else:
        merged.type = "user_interaction"  # type: ignore[assignment]

    # Accumulators
    latest_atomic_numbers: Optional[list[int]] = None
    latest_positions: Optional[np.ndarray] = None
    latest_velocities: Optional[np.ndarray] = None
    latest_cell: Optional[np.ndarray] = None
    latest_seq: Optional[int] = None
    latest_ack_val: Optional[int] = None
    latest_uic: Optional[int] = None
    latest_sim_step: Optional[int] = None

    for m in uis:
        s = int(getattr(m, "seq", 0) or 0) or None
        a = int(getattr(m, "ack", 0) or 0) or None
        u = int(getattr(m, "user_interaction_count", 0) or 0) or None
        st = int(getattr(m, "sim_step", 0) or 0) or None

        if s is not None:
            latest_seq = s
        if a is not None:
            latest_ack_val = a
        if u is not None:
            latest_uic = u
        if st is not None:
            latest_sim_step = st

        an = list(getattr(m, "atomic_numbers", []) or [])
        if an:
            latest_atomic_numbers = [int(z) for z in an]

        if getattr(m, "positions", None) and len(m.positions) > 0:
            latest_positions = _pb_vec3s_to_np(m.positions)
        if getattr(m, "velocities", None) and len(m.velocities) > 0:
            latest_velocities = _pb_vec3s_to_np(m.velocities)
        if getattr(m, "HasField", None) and m.HasField("cell"):
            latest_cell = _mat3_to_np(m.cell)

    if latest_seq is not None:
        merged.seq = int(latest_seq)
    if latest_ack_val is not None:
        merged.ack = int(latest_ack_val)
    if latest_uic is not None:
        merged.user_interaction_count = int(latest_uic)
    if latest_sim_step is not None:
        merged.sim_step = int(latest_sim_step)

    del merged.atomic_numbers[:]
    if latest_atomic_numbers is not None:
        merged.atomic_numbers.extend(int(z) for z in latest_atomic_numbers)

    _np_to_pb_vec3s(latest_positions, merged.positions)
    _np_to_pb_vec3s(latest_velocities, merged.velocities)

    cell_pb = _np_to_mat3(latest_cell)
    if cell_pb is not None:
        merged.cell.CopyFrom(cell_pb)

    return merged


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
            f"[ws:init] pool_size={pool_size} uma_handle={'yes' if predict_handle is not None else 'no'}",
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

        # Concurrency: producer (simulation loop) and consumer (message loop)
        stop_evt = asyncio.Event()

        # ---- minimal raw-frame queue you can peek into ----
        msg_buf: deque[pb.ClientAction] = deque()
        buf_has_data = asyncio.Event()
        reader_done = False

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
            # Debug log
            if self._ws_debug:
                print(
                    "[ws:tx] "
                    f"seq={seq} client_seq={client_seq} "
                    f"uic={user_interaction_count} step={sim_step} "
                    "has={"
                    f"pos:{positions is not None}, "
                    f"vel:{velocities is not None}, "
                    f"forces:{forces is not None}, "
                    f"cell:{cell is not None}, "
                    f"stress:{stress is not None}"
                    "} "
                    f"energy={energy if energy is not None else 'na'} "
                    f"msg={message or ''}",
                    flush=True,
                )

            # ---- build protobuf ----
            msg = pb.ServerResult()
            msg.seq = int(seq)
            if client_seq is not None:
                msg.client_seq = int(client_seq)
            if user_interaction_count is not None:
                msg.user_interaction_count = int(user_interaction_count)
            if sim_step is not None:
                msg.sim_step = int(sim_step)

            # Vec3 arrays
            _np_to_pb_vec3s(positions, msg.positions)
            _np_to_pb_vec3s(velocities, msg.velocities)
            _np_to_pb_vec3s(forces, msg.forces)

            # Cell / Stress (3x3)
            cell_pb = _np_to_mat3(cell)
            if cell_pb is not None:
                msg.cell.CopyFrom(cell_pb)

            stress_pb = _np_to_mat3(stress)
            if stress_pb is not None:
                msg.stress.CopyFrom(stress_pb)

            if message:
                msg.message = str(message)

            if simulation_stopped is True:
                msg.simulation_stopped = True

            if energy is not None:
                msg.energy = float(energy)

            # ---- send ----
            await ws.send_bytes(msg.SerializeToString())

        async def sim_loop():
            import time

            stall_notice_last = 0.0
            stall_notice_interval = 0.25  # seconds
            try:
                while not stop_evt.is_set():
                    # Commit any staged user inputs
                    if state.user_input_atomic_numbers is not None:
                        state.atomic_numbers = state.user_input_atomic_numbers
                        state.user_input_atomic_numbers = None
                    if state.user_input_positions is not None:
                        state.positions = np.asarray(
                            state.user_input_positions, dtype=np.float64
                        )
                        state.user_input_positions = None
                    if state.user_input_velocities is not None:
                        state.velocities = np.asarray(
                            state.user_input_velocities, dtype=np.float64
                        )
                        state.user_input_velocities = None
                    if state.user_input_cell is not None:
                        state.cell = np.asarray(state.user_input_cell, dtype=np.float64)
                        state.user_input_cell = None

                    if not state.running or state.atomic_numbers == []:
                        await asyncio.sleep(0.02)
                        continue

                    # Backpressure: don't exceed last_ack + 10
                    if state.server_seq - state.client_ack >= max_unacked:
                        now = time.time()
                        if now - stall_notice_last >= stall_notice_interval:
                            stall_notice_last = now
                            if self._ws_debug:
                                delta = state.server_seq - state.client_ack
                                print(
                                    f"[ws] WAITING_FOR_ACK send server_seq={state.server_seq} "
                                    f"client_ack={state.client_ack} delta={delta}",
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

                    # Snapshot UIC at the beginning of this step
                    uic_snapshot = state.user_interaction_count

                    worker = self._pool.any()
                    if state.sim_type == "md":
                        if self._log_calls:
                            print(
                                f"[ws:sim][MD] steps=1 T={state.params.temperature} "
                                f"dt={state.params.timestep_fs} friction={state.params.friction} "
                                f"natoms={len(state.atomic_numbers)}",
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
                                f"[ws:sim][RELAX] steps=1 fmax={state.params.fmax} "
                                f"max_step={state.params.max_step} calc={state.params.calculator} "
                                f"natoms={len(state.atomic_numbers)}",
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
                    state.positions = np.asarray(res["positions"], dtype=np.float64)

                    v_out = res.get("velocities")
                    if v_out is not None:
                        state.velocities = np.asarray(v_out, dtype=np.float64)

                    f_out = res.get("forces")
                    if f_out is not None:
                        F_np = np.asarray(f_out, dtype=np.float64)
                        if F_np.ndim == 1:
                            if F_np.size % 3 != 0:
                                # Malformed; drop instead of corrupting state
                                F_np = None
                            else:
                                F_np = F_np.reshape(-1, 3)
                        elif F_np.ndim == 2 and F_np.shape[1] != 3:
                            F_np = None
                        state.forces = F_np

                    energy_out = float(res.get("final_energy"))

                    state.server_seq += 1
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
                            f"[ws] sent frame seq={state.server_seq} sim_step={state.sim_step}",
                            flush=True,
                        )
            except Exception as e:
                print(f"[ws:sim_loop:error] {e}", flush=True)
                stop_evt.set()

        # --- recv_loop helpers: per-message-type handlers ---
        def _handle_ack_and_client_seq(msg) -> None:
            """Apply ACK/client_seq from incoming message with debug logs."""
            # client_seq (track the max we've seen)
            state.client_seq = max(state.client_seq, int(getattr(msg, "seq", 0) or 0))

            # ack (prefer HasField if available, else use nonzero value)
            prev_ack = int(state.client_ack)
            if hasattr(msg, "HasField") and msg.HasField("ack"):
                state.client_ack = max(
                    state.client_ack, int(getattr(msg, "ack", 0) or 0)
                )
            else:
                a = int(getattr(msg, "ack", 0) or 0)
                if a:
                    state.client_ack = max(state.client_ack, a)

            if self._ws_debug:
                if state.client_ack and state.client_ack > prev_ack:
                    delta = state.server_seq - state.client_ack
                    print(
                        f"[ws] ACK recv client_ack={state.client_ack} server_seq={state.server_seq} delta={delta}",
                        flush=True,
                    )
                elif state.client_ack:
                    print(f"[ws:rx][ack] client_ack={state.client_ack}", flush=True)

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
            # Stage incoming deltas into the "user_input_*" slots
            if getattr(msg, "atomic_numbers", None):
                if len(msg.atomic_numbers) > 0:
                    state.user_input_atomic_numbers = [
                        int(z) for z in msg.atomic_numbers
                    ]

            if getattr(msg, "positions", None) and len(msg.positions) > 0:
                state.user_input_positions = _pb_vec3s_to_np(msg.positions)

            if getattr(msg, "velocities", None) and len(msg.velocities) > 0:
                state.user_input_velocities = _pb_vec3s_to_np(msg.velocities)

            if getattr(msg, "HasField", None) and msg.HasField("cell"):
                state.user_input_cell = _mat3_to_np(msg.cell)

            n = (
                0
                if state.user_input_positions is None
                else int(state.user_input_positions.shape[0])
            )

            if self._ws_debug:
                print(
                    f"[ws] USER_INTERACTION recv natoms={n} running={'true' if state.running else 'false'}",
                    flush=True,
                )

            if state.running:
                # Invalidate forces so next produced frame recomputes
                state.forces = None
                if self._ws_debug:
                    print(
                        "[ws:state] applied USER_INTERACTION during running sim",
                        flush=True,
                    )
                return

            # Idle: apply staged user inputs to the live state (then clear the staging)
            if self._ws_debug:
                print("[ws] USER_INTERACTION idle -> compute", flush=True)

            if state.user_input_atomic_numbers is not None:
                state.atomic_numbers = list(state.user_input_atomic_numbers)
                state.user_input_atomic_numbers = None
                assert (
                    len(state.atomic_numbers) != 0
                ), "atomic_numbers must be non-empty"

            if state.user_input_positions is not None:
                state.positions = np.asarray(
                    state.user_input_positions, dtype=np.float64
                )
                state.user_input_positions = None

            if state.user_input_velocities is not None:
                state.velocities = np.asarray(
                    state.user_input_velocities, dtype=np.float64
                )
                state.user_input_velocities = None

            if state.user_input_cell is not None:
                state.cell = np.asarray(state.user_input_cell, dtype=np.float64)
                state.user_input_cell = None

            worker = self._pool.any()
            try:
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
            except Exception as e:
                print(
                    f"[ws] USER_INTERACTION compute failed with '{state.params.calculator}' (no fallback): {e}",
                    flush=True,
                )
                state.server_seq += 1
                await _send_result_bytes(
                    seq=state.server_seq,
                    client_seq=state.client_seq,
                    user_interaction_count=(
                        int(uic_in_msg) if uic_in_msg is not None else None
                    ),
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

            results = res.get("results", {}) if isinstance(res, dict) else {}
            E = results.get("energy")
            F = results.get("forces")
            S = results.get("stress")

            # Numpy coercions with shape checks
            F_np = None
            if F is not None:
                try:
                    F_np = np.asarray(F, dtype=np.float64)
                    if F_np.ndim == 1:
                        if F_np.size % 3 == 0:
                            F_np = F_np.reshape(-1, 3)
                        else:
                            raise ValueError("forces length is not divisible by 3")
                    if F_np.shape[1] != 3:
                        raise ValueError("forces must have shape (N, 3)")
                except Exception:
                    F_np = None  # fall back to None if malformed

            S_np = None
            if S is not None:
                try:
                    S_np = np.asarray(S, dtype=np.float64)
                    if S_np.shape != (3, 3):
                        S_np = None
                except Exception:
                    S_np = None

            state.server_seq += 1
            await _send_result_bytes(
                seq=state.server_seq,
                client_seq=state.client_seq,
                user_interaction_count=(
                    int(uic_in_msg) if uic_in_msg is not None else None
                ),
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
                    f"[ws] START_SIM type={state.sim_type} T={state.params.temperature} "
                    f"dt={state.params.timestep_fs} friction={state.params.friction}",
                    flush=True,
                )
                print(
                    f"[ws:state] counters uic={state.user_interaction_count} sim_step={state.sim_step}",
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
            async def read_next_message_from_socket() -> Optional[pb.ClientAction]:
                """
                Return the next message to process.
                - If the head is USER_INTERACTION, coalesce all consecutive UI frames.
                - If the head is non-UI, return it as-is.
                - Block briefly if empty until data arrives or reader finishes.
                - Return None once reader is done and the buffer is empty.
                """
                while True:
                    if msg_buf:
                        head = msg_buf[0]
                        if _is_user_interaction(head):
                            msg = _coalesce_user_interactions(msg_buf)
                        else:
                            msg = msg_buf.popleft()
                        if not msg_buf:
                            buf_has_data.clear()
                        return msg

                    if reader_done:
                        return None

                    try:
                        await asyncio.wait_for(buf_has_data.wait(), timeout=0.05)
                    except asyncio.TimeoutError:
                        continue
                    finally:
                        if not msg_buf:
                            buf_has_data.clear()

            try:
                while True:
                    msg = await read_next_message_from_socket()
                    if msg is None:
                        break

                    # shared header handling
                    _handle_ack_and_client_seq(msg)
                    uic_in_msg = _extract_correlation_fields(msg)

                    # dispatch
                    tname = _type_name(msg)
                    if tname == "USER_INTERACTION":
                        await _handle_user_interaction(msg, uic_in_msg)
                    elif tname == "START_SIMULATION":
                        _handle_start_simulation(msg)
                    elif tname == "STOP_SIMULATION":
                        await _handle_stop_simulation()
                    elif tname == "PING":
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
            return max(1, torch.cuda.device_count())
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

    # If a specific HTTP port is requested, start Serve with that port
    if http_port is not None:
        try:
            serve.start(
                detached=False,
                http_options={"host": "0.0.0.0", "port": int(http_port)},
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
        logging_config=LoggingConfig(enable_access_log=False, log_level="ERROR"),
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
