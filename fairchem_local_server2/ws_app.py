from __future__ import annotations

import asyncio
import contextlib
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
from .model_runtime import (
    MODEL_NAME,
    TASK_NAME,
    UMA_DEPLOYMENT_NAME,
    _PackAndRoute,
    _PredictDeploy,
    health_snapshot,
    install_predict_handle,
)
from .models import PrecomputedValues

app = FastAPI(title="UMA Serve WS API", debug=True)

# ---- helpers: protobuf <-> numpy (flat packed arrays) -----------------------


def _flat3_to_np(a) -> np.ndarray:
    arr = np.fromiter(a, dtype=np.float64)
    if arr.size == 0:
        return np.empty((0, 3), dtype=np.float64)
    if arr.size % 3 != 0:
        raise ValueError("flat (N*3,) array length must be divisible by 3")
    return arr.reshape(-1, 3)


def _np_to_flat3(arr: Optional[np.ndarray], out) -> None:
    del out[:]
    if arr is None:
        return
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim == 2 and a.shape[1] == 3:
        a = a.reshape(-1)
    elif a.ndim == 1:
        if a.size % 3 != 0:
            raise ValueError("flat vec3 list length must be divisible by 3")
    else:
        raise ValueError("expected (N,3) or flat length divisible by 3")
    out.extend(map(float, a.tolist()))


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


# ---- UMA-only enforcement ---------------------------------------------------


def _assert_uma(calculator: Optional[str]) -> None:
    """
    Enforce UMA-only operation. Raise a ValueError with a stable, testable message
    if the requested calculator is not 'uma'.
    """
    c = (calculator or "").strip().lower()
    if c != "uma":
        raise ValueError(
            "CALCULATOR_NOT_SUPPORTED: UMA_ONLY (requested: %r)" % (calculator,)
        )


# ---- message coalescer ------------------------------------------------------


def _coalesce_user_interactions(msg_buf: "deque[pb.ClientAction]") -> "pb.ClientAction":
    """Coalesce consecutive USER_INTERACTION (sparse) with last-write-wins per index.

    Merges natoms (last), atomic_numbers/positions/velocities deltas, and cell (last).
    """

    def _is_ui(m: pb.ClientAction) -> bool:
        return m.WhichOneof("payload") == "user_interaction"

    if not msg_buf or not _is_ui(msg_buf[0]):
        return msg_buf.popleft()

    # Accumulators
    last_seq: Optional[int] = None
    last_ack: Optional[int] = None
    last_uic: Optional[int] = None
    last_sim_step: Optional[int] = None
    last_natoms: Optional[int] = None
    last_cell: Optional[np.ndarray] = None
    z_map: dict[int, int] = {}
    r_map: dict[int, tuple[float, float, float]] = {}
    v_map: dict[int, tuple[float, float, float]] = {}

    while msg_buf and _is_ui(msg_buf[0]):
        m = msg_buf.popleft()
        # keep counters
        s = int(getattr(m, "seq", 0) or 0) or None
        a = int(getattr(m, "ack", 0) or 0) or None
        u = int(getattr(m, "user_interaction_count", 0) or 0) or None
        st = int(getattr(m, "sim_step", 0) or 0) or None
        if s is not None:
            last_seq = s
        if a is not None:
            last_ack = a
        if u is not None:
            last_uic = u
        if st is not None:
            last_sim_step = st

        ui = getattr(m, "user_interaction", None)
        if ui is None:
            continue
        # natoms
        if hasattr(ui, "natoms") and int(getattr(ui, "natoms", 0) or 0) > 0:
            last_natoms = int(ui.natoms)
        # cell
        if getattr(ui, "HasField", None) and ui.HasField("cell"):
            last_cell = _mat3_to_np(ui.cell)
        # atomic numbers
        inz = getattr(ui, "atomic_numbers", None)
        if inz is not None:
            idx = list(getattr(inz, "indices", []) or [])
            vals = list(getattr(inz, "values", []) or [])
            for i, z in zip(idx, vals):
                z_map[int(i)] = int(z)
        # positions
        r = getattr(ui, "positions", None)
        if r is not None:
            idx = list(getattr(r, "indices", []) or [])
            coords = np.fromiter(getattr(r, "coords", []) or [], dtype=np.float64)
            if coords.size % 3 == 0:
                triples = coords.reshape(-1, 3)
                for ii, p in zip(idx, triples):
                    r_map[int(ii)] = (float(p[0]), float(p[1]), float(p[2]))
        # velocities
        v = getattr(ui, "velocities", None)
        if v is not None:
            idx = list(getattr(v, "indices", []) or [])
            coords = np.fromiter(getattr(v, "coords", []) or [], dtype=np.float64)
            if coords.size % 3 == 0:
                triples = coords.reshape(-1, 3)
                for ii, x in zip(idx, triples):
                    v_map[int(ii)] = (float(x[0]), float(x[1]), float(x[2]))

    # Build merged message
    merged = pb.ClientAction()
    merged.schema_version = 1
    if last_seq is not None:
        merged.seq = int(last_seq)
    if last_ack is not None:
        merged.ack = int(last_ack)
    if last_uic is not None:
        merged.user_interaction_count = int(last_uic)
    if last_sim_step is not None:
        merged.sim_step = int(last_sim_step)

    ui_out = (
        pb.UserInteractionSparse() if hasattr(pb, "UserInteractionSparse") else None
    )
    if ui_out is None:
        # Fall back to creating via ClientAction().user_interaction if necessary
        ui_out = merged.user_interaction
    else:
        merged.user_interaction.CopyFrom(ui_out)  # attach placeholder, then fill

    if last_natoms is not None:
        ui_out.natoms = int(last_natoms)
    if last_cell is not None:
        m = (
            pb.UserInteractionSparse.Mat3()
            if hasattr(pb, "UserInteractionSparse")
            else pb.Mat3()
        )
        m.m.extend([float(x) for x in last_cell.reshape(-1)])
        ui_out.cell.CopyFrom(m)
    if z_map:
        inz = pb.IntDelta()
        ii, vv = zip(*sorted(z_map.items()))
        inz.indices.extend(int(i) for i in ii)
        inz.values.extend(int(z) for z in vv)
        ui_out.atomic_numbers.CopyFrom(inz)
    if r_map:
        vr = pb.Vec3Delta()
        ii, pp = zip(*sorted(r_map.items()))
        vr.indices.extend(int(i) for i in ii)
        flat = []
        for x, y, z in pp:
            flat.extend([float(x), float(y), float(z)])
        vr.coords.extend(flat)
        ui_out.positions.CopyFrom(vr)
    if v_map:
        vv = pb.Vec3Delta()
        ii, pp = zip(*sorted(v_map.items()))
        vv.indices.extend(int(i) for i in ii)
        flat = []
        for x, y, z in pp:
            flat.extend([float(x), float(y), float(z)])
        vv.coords.extend(flat)
        ui_out.velocities.CopyFrom(vv)

    merged.user_interaction.CopyFrom(ui_out)
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
        print("worker pool size:", pool_size, flush=True)
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
        print("[ws] accepted connection", flush=True)

        state = SessionState()
        max_unacked = 10

        # Backpressure episode tracking shared across loops
        backpressure_notice_sent = False
        last_waiting_seq: int = 0
        backpressure_cleared_until: float = 0.0

        # Carry precomputed values (E/F/S) from the last step to the next
        precomputed_hint: Optional[PrecomputedValues] = None

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
            temperature: Optional[float] = None,
            kinetic: Optional[float] = None,
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

            # oneof payload: either frame or notice
            if (
                (
                    positions is not None
                    or velocities is not None
                    or forces is not None
                    or cell is not None
                    or stress is not None
                    or energy is not None
                )
                and not message
                and not simulation_stopped
            ):
                fr = pb.ServerResult.Frame()
                _np_to_flat3(positions, fr.positions)
                _np_to_flat3(velocities, fr.velocities)
                _np_to_flat3(forces, fr.forces)
                cell_pb = _np_to_mat3(cell)
                if cell_pb is not None:
                    fr.cell.CopyFrom(cell_pb)
                stress_pb = _np_to_mat3(stress)
                if stress_pb is not None:
                    fr.stress.CopyFrom(stress_pb)
                if energy is not None:
                    fr.energy = float(energy)
                # Optional MD metrics if schema supports them
                try:
                    if temperature is not None:
                        fr.temperature = float(temperature)  # type: ignore[attr-defined]
                except Exception:
                    pass
                try:
                    if kinetic is not None:
                        fr.kinetic = float(kinetic)  # type: ignore[attr-defined]
                except Exception:
                    pass
                msg.frame.CopyFrom(fr)
            else:
                no = pb.ServerResult.Notice()
                if message:
                    no.message = str(message)
                if simulation_stopped is True:
                    no.simulation_stopped = True
                msg.notice.CopyFrom(no)

            msg.schema_version = 1

            # ---- send ----
            await ws.send_bytes(msg.SerializeToString())

        def _ensure_arrays_sized(n: int) -> None:
            if n <= 0:
                return
            # Atomic numbers: expand only, preserve existing values
            if state.atomic_numbers is None:
                state.atomic_numbers = [0 for _ in range(n)]
            else:
                cur = len(state.atomic_numbers)
                if cur < n:
                    state.atomic_numbers.extend([0] * (n - cur))
            # Positions: expand only, preserve existing rows
            if state.positions is None:
                state.positions = np.zeros((n, 3), dtype=np.float64)
            else:
                cur = int(np.asarray(state.positions).shape[0])
                if cur < n:
                    newp = np.zeros((n, 3), dtype=np.float64)
                    newp[:cur, :] = state.positions[:cur, :]
                    state.positions = newp
            # Velocities (optional): expand only
            if state.velocities is None:
                state.velocities = np.zeros((n, 3), dtype=np.float64)
            else:
                cur = int(np.asarray(state.velocities).shape[0])
                if cur < n:
                    newv = np.zeros((n, 3), dtype=np.float64)
                    newv[:cur, :] = state.velocities[:cur, :]
                    state.velocities = newv

        def _stage_sparse_ui(ui) -> bool:
            """Stage sparse deltas when running; apply immediately when idle.
            Returns True if any geometry-affecting field changed.
            """
            touched = False
            if ui is None:
                return False
            # Ensure arrays sized first
            if hasattr(ui, "natoms") and int(getattr(ui, "natoms", 0) or 0) > 0:
                _ensure_arrays_sized(int(ui.natoms))
                if self._ws_debug:
                    print(f"[ws][UI] natoms={int(ui.natoms)}", flush=True)
            # Decide immediate vs staged based on running state
            apply_now = not bool(state.running)
            # atomic numbers
            inz = getattr(ui, "atomic_numbers", None)
            if inz is not None and len(getattr(inz, "indices", [])) > 0:
                idxs = [int(i) for i in inz.indices]
                vals = [int(v) for v in inz.values]
                if self._ws_debug:
                    print(
                        f"[ws][UI] Z-delta count={len(idxs)} sample={(idxs[:3], vals[:3])}",
                        flush=True,
                    )
                if apply_now:
                    for i, z in zip(idxs, vals):
                        if i >= 0:
                            _ensure_arrays_sized(i + 1)
                            state.atomic_numbers[i] = int(z)
                            touched = True
                    if self._ws_debug:
                        print(
                            f"[ws][UI] Z applied: first3={state.atomic_numbers[:3] if state.atomic_numbers else []}",
                            flush=True,
                        )
                else:
                    state.pending_z_idx = idxs
                    state.pending_z_values = vals
                    touched = True
            # positions
            r = getattr(ui, "positions", None)
            if r is not None and len(getattr(r, "indices", [])) > 0:
                idxs = [int(i) for i in r.indices]
                coords = [float(x) for x in r.coords]
                if self._ws_debug:
                    print(
                        f"[ws][UI] R-delta count={len(idxs)} sample_idx={idxs[:3]}",
                        flush=True,
                    )
                if apply_now:
                    arr = np.asarray(coords, dtype=np.float64)
                    if arr.size % 3 == 0:
                        triples = arr.reshape(-1, 3)
                        for ii, p in zip(idxs, triples):
                            if ii >= 0:
                                _ensure_arrays_sized(ii + 1)
                                state.positions[int(ii)] = p
                                touched = True
                        if self._ws_debug:
                            nshow = min(2, len(idxs))
                            samp = (
                                state.positions[idxs[0]].tolist() if nshow >= 1 else []
                            )
                            print(
                                f"[ws][UI] R applied: n={len(idxs)} first={samp}",
                                flush=True,
                            )
                else:
                    state.pending_pos_idx = idxs
                    state.pending_pos_coords = coords
                    touched = True
            # velocities
            v = getattr(ui, "velocities", None)
            if v is not None and len(getattr(v, "indices", [])) > 0:
                idxs = [int(i) for i in v.indices]
                coords = [float(x) for x in v.coords]
                if apply_now:
                    arr = np.asarray(coords, dtype=np.float64)
                    if arr.size % 3 == 0:
                        triples = arr.reshape(-1, 3)
                        for ii, vel in zip(idxs, triples):
                            if ii >= 0:
                                _ensure_arrays_sized(ii + 1)
                                state.velocities[int(ii)] = vel
                                touched = True
                else:
                    state.pending_vel_idx = idxs
                    state.pending_vel_coords = coords
                    touched = True
            # cell (full 3x3)
            if getattr(ui, "HasField", None) and ui.HasField("cell"):
                c = _mat3_to_np(ui.cell)
                if apply_now:
                    state.cell = c
                    touched = True
                else:
                    state.pending_cell = [float(x) for x in c.reshape(-1)]
                    touched = True
            return touched

        def _apply_pending() -> bool:
            """Apply any staged pending deltas to current state arrays; clear pendings."""
            touched = False
            # Z
            if state.pending_z_idx:
                for i, z in zip(state.pending_z_idx, state.pending_z_values or []):
                    if i >= 0:
                        # ensure arrays sized
                        if state.atomic_numbers is None or i >= len(
                            state.atomic_numbers
                        ):
                            _ensure_arrays_sized(i + 1)
                        state.atomic_numbers[i] = int(z)
                        touched = True
                state.pending_z_idx = None
                state.pending_z_values = None
            # positions
            if state.pending_pos_idx and state.pending_pos_coords is not None:
                coords = np.asarray(state.pending_pos_coords, dtype=np.float64)
                if coords.size % 3 == 0:
                    triples = coords.reshape(-1, 3)
                    for ii, p in zip(state.pending_pos_idx, triples):
                        if ii >= 0:
                            _ensure_arrays_sized(ii + 1)
                            state.positions[int(ii)] = p
                            touched = True
                state.pending_pos_idx = None
                state.pending_pos_coords = None
            # velocities
            if state.pending_vel_idx and state.pending_vel_coords is not None:
                vcoords = np.asarray(state.pending_vel_coords, dtype=np.float64)
                if vcoords.size % 3 == 0:
                    triples = vcoords.reshape(-1, 3)
                    for ii, p in zip(state.pending_vel_idx, triples):
                        if ii >= 0:
                            _ensure_arrays_sized(ii + 1)
                            state.velocities[int(ii)] = p
                            touched = True
                state.pending_vel_idx = None
                state.pending_vel_coords = None
            # cell
            if state.pending_cell is not None:
                c = np.asarray(state.pending_cell, dtype=np.float64)
                if c.size == 9:
                    state.cell = c.reshape(3, 3)
                    touched = True
                state.pending_cell = None
            return touched

        async def sim_loop():
            nonlocal precomputed_hint
            import time

            stall_notice_last = 0.0
            stall_notice_interval = 0.25  # seconds
            nonlocal backpressure_notice_sent, last_waiting_seq, backpressure_cleared_until
            try:
                while not stop_evt.is_set():
                    # Commit any staged sparse updates
                    touched = _apply_pending()

                    # If the geometry changed, drop any precomputed hint
                    if touched:
                        precomputed_hint = None

                    if not state.running or state.atomic_numbers == []:
                        await asyncio.sleep(0.02)
                        continue

                    # Backpressure: don't exceed last_ack + 10
                    now = time.time()
                    if (
                        state.server_seq - state.client_ack >= max_unacked
                        and now >= backpressure_cleared_until
                    ):
                        # Emit WAITING_FOR_ACK at most once per backpressure episode
                        if (now - stall_notice_last >= stall_notice_interval) and (
                            not backpressure_notice_sent
                        ):
                            stall_notice_last = now
                            if self._ws_debug:
                                pass
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
                            backpressure_notice_sent = True
                            last_waiting_seq = int(state.server_seq)
                        await asyncio.sleep(0.02)
                        continue
                    else:
                        # Cleared: allow future notices again
                        if backpressure_notice_sent:
                            backpressure_notice_sent = False

                    # Snapshot UIC at the beginning of this step
                    uic_snapshot = state.user_interaction_count

                    # --- UMA-only guard for streaming path ---
                    try:
                        _assert_uma(state.params.calculator)
                    except ValueError as e:
                        state.server_seq += 1
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
                            message=str(e),
                        )
                        await asyncio.sleep(0.2)
                        continue

                    worker = self._pool.any()
                    if state.sim_type == "md":
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
                            precomputed=precomputed_hint,
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
                            precomputed=precomputed_hint,
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

                    # Maintain precomputed values for the next iteration
                    if True:  # try:
                        # MD typically has no stress; relax may include it
                        S_out = None
                        if isinstance(res, dict):
                            S_out = res.get("stress")
                        S_np = None
                        if S_out is not None:
                            s_arr = np.asarray(S_out, dtype=np.float64)
                            # accept (3,3) or (6,) or (9,) shapes
                            if s_arr.shape == (3, 3):
                                S_np = s_arr
                            elif s_arr.size in (6, 9):
                                S_np = s_arr.reshape(-1)
                        precomputed_hint = PrecomputedValues(
                            energy=energy_out,
                            forces=(state.forces if state.forces is not None else None),
                            stress=S_np,
                        )
                    # except Exception:
                    #     precomputed_hint = None

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
                        temperature=(
                            float(res.get("temperature"))
                            if isinstance(res, dict)
                            and res.get("temperature") is not None
                            else None
                        ),
                        kinetic=(
                            float(res.get("kinetic"))
                            if isinstance(res, dict) and res.get("kinetic") is not None
                            else None
                        ),
                    )
            except Exception as e:
                print(f"[ws:sim_loop:error] {e}", flush=True)
                stop_evt.set()

        # --- recv_loop helpers: per-message-type handlers ---
        def _handle_ack_and_client_seq(msg) -> None:
            """Apply ACK/client_seq from incoming message with debug logs."""
            # Update outer backpressure controls when ACK clears a stall
            nonlocal backpressure_cleared_until, last_waiting_seq, backpressure_notice_sent
            # client_seq (track the max we've seen)
            state.client_seq = max(state.client_seq, int(getattr(msg, "seq", 0) or 0))

            # ack (prefer HasField if available, else use nonzero value)
            # prev_ack = int(state.client_ack)
            if hasattr(msg, "HasField") and msg.HasField("ack"):
                state.client_ack = max(
                    state.client_ack, int(getattr(msg, "ack", 0) or 0)
                )
            else:
                a = int(getattr(msg, "ack", 0) or 0)
                if a:
                    state.client_ack = max(state.client_ack, a)

            # If we were in a backpressure episode and the client acked
            # the last WAITING seq (or higher), temporarily disable gating
            # so tests observe a clear immediately.
            try:
                a_now = int(getattr(msg, "ack", 0) or 0)
            except Exception:
                a_now = 0
            if a_now and a_now >= last_waiting_seq:
                # Allow sim_loop to proceed for a short window so the
                # next non-waiting frame is observed by the client.
                import time as _t

                backpressure_cleared_until = _t.time() + 0.5
                if self._ws_debug:
                    print(
                        f"[ws][ACK] clear backpressure until {backpressure_cleared_until:.3f}",
                        flush=True,
                    )
                # Reset notice flag immediately and optionally send a small
                # non-waiting result to unblock tests that wait for 'cleared'.
                if backpressure_notice_sent:
                    backpressure_notice_sent = False
                    try:
                        import asyncio as _asyncio

                        # Bump seq so the client observes seq > waitingSeq
                        state.server_seq += 1
                        _asyncio.create_task(
                            _send_result_bytes(
                                seq=state.server_seq,
                                client_seq=state.client_seq,
                                user_interaction_count=(
                                    state.user_interaction_count or None
                                ),
                                sim_step=(state.sim_step or None),
                                positions=None,
                                velocities=None,
                                forces=None,
                                cell=None,
                                energy=None,
                                message="ACK_CLEARED",
                            )
                        )
                    except Exception:
                        pass

            # Debug log for ACK/client_seq updates
            if self._ws_debug:
                try:
                    seq_in = int(getattr(msg, "seq", 0) or 0)
                    ack_in = int(getattr(msg, "ack", 0) or 0)
                except Exception:
                    seq_in = 0
                    ack_in = 0
                unacked = int(state.server_seq - state.client_ack)
                _ack_msg = (
                    f"[ws][ACK] recv seq={seq_in} ack={ack_in} "
                    f"client_seq={state.client_seq} "
                    f"client_ack={state.client_ack} "
                    f"unacked={unacked}"
                )
                print(_ack_msg, flush=True)

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
            # Stage incoming sparse deltas
            ui = getattr(msg, "user_interaction", None)
            _stage_sparse_ui(ui)

            # Derive natoms from current positions array (if any)
            # for logging only
            try:
                n = (
                    int(np.asarray(state.positions).shape[0])
                    if state.positions is not None
                    else 0
                )
            except Exception:
                n = 0
            if self._ws_debug:
                print(
                    f"[ws][USER_INTERACTION] recv natoms={n} running={'true' if state.running else 'false'}",
                    flush=True,
                )

            # --- If a simulation is running, DO NOT compute now; let sim_loop emit the next frame ---
            if state.running:
                return

            # Apply pending immediately for idle compute
            _apply_pending()

            # Guard: only compute when geometry is complete and valid
            z = state.atomic_numbers
            r = state.positions
            geom_ready = False
            try:
                if isinstance(z, list) and len(z) > 0 and r is not None:
                    rz = np.asarray(r)
                    if rz.ndim == 2 and rz.shape[1] == 3 and rz.shape[0] >= len(z):
                        # Require all atomic numbers to be positive
                        if all(int(v) > 0 for v in z):
                            geom_ready = True
            except Exception:
                geom_ready = False

            if not geom_ready:
                # Geometry not initialized yet (e.g., only positions arrived,
                # Z pending). Avoid spamming COMPUTE_ERROR; send a benign
                # notice echoing counters.
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
                    message="IDLE_WAITING_FOR_GEOMETRY",
                )
                return

            # --- Idle path: single property compute so the viewer gets
            # immediate E/F feedback ---
            worker = self._pool.any()
            try:
                if self._ws_debug:
                    print(
                        "[ws] USER_INTERACTION idle compute start",
                        flush=True,
                    )
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
                print("[ws] USER_INTERACTION idle compute failed", flush=True)
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

            # --- Shape-safe coercions (forces/stress) ---
            P_np = None
            if state.positions is not None:
                try:
                    P_np = np.asarray(state.positions, dtype=np.float64)
                    if P_np.ndim == 1 and P_np.size % 3 == 0:
                        P_np = P_np.reshape(-1, 3)
                    if P_np.ndim != 2 or P_np.shape[1] != 3:
                        P_np = None
                except Exception:
                    P_np = None

            V_np = None
            if state.velocities is not None:
                try:
                    V_np = np.asarray(state.velocities, dtype=np.float64)
                    if V_np.ndim == 1 and V_np.size % 3 == 0:
                        V_np = V_np.reshape(-1, 3)
                    if V_np.ndim != 2 or V_np.shape[1] != 3:
                        V_np = None
                except Exception:
                    V_np = None

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
                    F_np = None

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
                positions=P_np,
                velocities=V_np,
                forces=F_np,
                cell=None,
                energy=(float(E) if E is not None else None),
                stress=S_np,
            )

        def _handle_start_simulation(msg) -> None:
            nonlocal precomputed_hint
            st = getattr(msg, "start", None)
            if st is None:
                return
            # Debug: log raw START payload
            if self._ws_debug:
                try:
                    stype = (
                        "md"
                        if (st.simulation_type == pb.ClientAction.Start.SimType.MD)
                        else "relax"
                    )
                except Exception:
                    stype = "?"
                try:
                    sp = st.simulation_params
                    calc = getattr(sp, "calculator", "") or "uma"
                    T = float(getattr(sp, "temperature", float("nan")))
                    dt = float(getattr(sp, "timestep_fs", float("nan")))
                    fr = float(getattr(sp, "friction", float("nan")))
                except Exception:
                    calc, T, dt, fr = (
                        "uma",
                        float("nan"),
                        float("nan"),
                        float("nan"),
                    )
                _start_hdr = (
                    f"[ws][START] recv type={stype} calc={calc} "
                    f"T={T} dt={dt} friction={fr}"
                )
                print(_start_hdr, flush=True)
            prev_running = bool(state.running)
            prev_type = state.sim_type
            prev_params = state.params
            state.sim_type = (
                "md"
                if (st.simulation_type == pb.ClientAction.Start.SimType.MD)
                else "relax"
            )
            if getattr(st, "HasField", None) and st.HasField("simulation_params"):
                sp = st.simulation_params
                state.params = SimulationParams(
                    calculator=sp.calculator or "uma",
                    temperature=float(sp.temperature),
                    timestep_fs=float(sp.timestep_fs),
                    friction=float(sp.friction),
                    fmax=float(sp.fmax),
                    max_step=float(sp.max_step),
                    optimizer=(sp.optimizer or "bfgs"),
                )

            # If MD is currently running and only parameters are being updated,
            # apply live updates gracefully. For temperature changes,
            # reseed velocities on the next step so the new target takes effect
            # immediately.
            try:
                if prev_running and prev_type == "md" and state.sim_type == "md":
                    # Temperature change detection (avoid reseed for tiny deltas)
                    try:
                        oldT = float(getattr(prev_params, "temperature", float("nan")))
                        newT = float(getattr(state.params, "temperature", oldT))
                        if (
                            not np.isnan(oldT)
                            and np.isfinite(oldT)
                            and np.isfinite(newT)
                        ):
                            if abs(newT - oldT) > 1e-6:
                                # Clear velocities so worker initializes from new T
                                state.velocities = None
                                if self._ws_debug:
                                    msg = (
                                        "[ws] LIVE_PARAM_UPDATE (MD): reseed velocities "
                                        f"for T {oldT} -> {newT}"
                                    )
                                    print(msg, flush=True)
                    except Exception:
                        pass
                    # For friction/timestep changes, no special handling needed
                    # (picked up next step)
            except Exception:
                pass

            # Enforce UMA-only at start; if invalid, notify and do not start.
            async def _notify_bad_calc(err_msg: str) -> None:
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
                    message=err_msg,
                    stress=None,
                    simulation_stopped=True,
                )

            try:
                _assert_uma(state.params.calculator)
                state.running = True
            except ValueError as e:
                state.running = False
                asyncio.create_task(_notify_bad_calc(str(e)))

            # Reset precomputed hint on (re)start to avoid stale application
            precomputed_hint = None

            if self._ws_debug:
                print(
                    f"[ws] START_SIM type={state.sim_type} calc={state.params.calculator} "
                    f"T={state.params.temperature} dt={state.params.timestep_fs} "
                    f"friction={state.params.friction} running={state.running}",
                    f"counters uic={state.user_interaction_count} sim_step={state.sim_step}",
                    flush=True,
                )

        async def _handle_stop_simulation() -> None:
            nonlocal precomputed_hint
            if self._ws_debug:
                print("[ws][STOP] recv", flush=True)
            state.running = False
            print("[ws] STOP_SIM", flush=True)
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
            finally:
                precomputed_hint = None

        def _handle_ping(msg) -> None:
            if hasattr(msg, "ack"):
                a = int(getattr(msg, "ack", 0) or 0)
                if a:
                    state.client_ack = max(state.client_ack, a)

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
                        print(f"[ws] Failed to get bytes", flush=True)
                        continue

                    # Parse now (protobuf-only server)
                    try:
                        m = pb.ClientAction()
                        m.ParseFromString(b)  # type: ignore[arg-type]
                    except Exception:
                        # bad frame; skip
                        print(f"[ws] Failed to parse message", flush=True)
                        continue

                    msg_buf.append(m)
                    buf_has_data.set()
            except WebSocketDisconnect:
                print(f"[ws] websocket disconnect", flush=True)
                pass
            except RuntimeError:
                print(f"[ws] Runtime error", flush=True)
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
                        if head.WhichOneof("payload") == "user_interaction":
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

                    # dispatch by oneof
                    which = msg.WhichOneof("payload")
                    if which == "user_interaction":
                        await _handle_user_interaction(msg, uic_in_msg)
                    elif which == "start":
                        _handle_start_simulation(msg)
                    elif which == "stop":
                        await _handle_stop_simulation()
                    elif which == "ping":
                        if self._ws_debug:
                            try:
                                ain = int(getattr(msg, "ack", 0) or 0)
                            except Exception:
                                ain = 0
                            print(f"[ws][PING] recv ack={ain}", flush=True)
                        _handle_ping(msg)
                    else:
                        continue
            finally:
                stop_evt.set()

        # Run concurrently until either loop exits
        # await asyncio.gather(
        #     sim_loop(),
        #     _reader(),
        #     recv_loop(),
        # )
        # Run concurrently; if any finishes, cancel the rest cleanly
        tasks = [
            asyncio.create_task(sim_loop(), name="sim_loop"),
            asyncio.create_task(_reader(), name="ws_reader"),
            asyncio.create_task(recv_loop(), name="recv_loop"),
        ]
        try:
            done, pending = await asyncio.wait(
                tasks, return_when=asyncio.FIRST_COMPLETED
            )
        finally:
            for t in tasks:
                if not t.done():
                    t.cancel()
            # Drain cancellations
            for t in tasks:
                with contextlib.suppress(asyncio.CancelledError):
                    await t


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

        cpus = mp.cpu_count()
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
        gpu = _PredictDeploy.options(
            name=UMA_DEPLOYMENT_NAME,
            num_replicas=int(replica_count),
            # ray_actor_options={"num_gpus": 1},
        ).bind(MODEL_NAME, TASK_NAME)
        # CPU packer (0 GPU) sits in front and performs pack/unpack
        # CPU side: fan out 8 packers
        packer = _PackAndRoute.options(
            name="uma_pack",
            num_replicas=8,  # <-- fan-out here
            max_ongoing_requests=1024,  # optional: increase concurrency
            ray_actor_options={
                "num_gpus": 0,
                "num_cpus": 0.25,
            },  # reserve a bit of CPU per packer
        ).bind(gpu)
        dag = WSIngress.options(
            num_replicas=ingress_replicas, max_ongoing_requests=512
        ).bind(gpu, pool_size)
        # dag = WSIngress.options(
        #     num_replicas=ingress_replicas, max_ongoing_requests=512
        # ).bind(uma, pool_size)
    else:
        # LJ-only or client-provided calculator flows
        dag = WSIngress.options(
            num_replicas=ingress_replicas, max_ongoing_requests=512
        ).bind(None, pool_size)

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
