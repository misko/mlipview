# Backend Architecture (current)

## Runtime topology
- **Entrypoint**: `fairchem_local_server2/ws_app.py` exposes a FastAPI app and is wrapped by a Ray Serve deployment (`WSIngress`). `serve_ws_app.py` boots the deployment graph.
- **Serve DAG**: when GPUs are available, `deploy()` binds
  - `_PredictDeploy` (GPU inference) and `_PackAndRoute` (CPU pack/unpack fan-out) for UMA predictions,
  - `WSIngress` replicas that host the FastAPI/WebSocket server and a per-replica `WorkerPool`.
- **Worker pool**: `WorkerPool` spins up Ray actors of class `ASEWorker` (remote CPU processes). Each exposes synchronous `run_simple`, `run_md`, and `run_relax` methods that wrap ASE calculators with UMA.
- **Calculator handle**: `_PredictDeploy` installs its Ray handle via `install_predict_handle()`. `model_runtime.get_calculator()` builds lightweight `FAIRChemCalculator` instances on demand so each worker can attach UMA to ASE `Atoms`.

## Session state and concurrency
- **Per-connection state**: `SessionState` tracks current geometry arrays (numpy), velocities, forces, simulation parameters, counters, and a simple `need_idle_emit` flag used to schedule idle recomputes.
- **Concurrency model**: each WebSocket spawns three async tasks:
  1. `_reader()` to pull protobuf frames into a deque.
  2. `recv_loop()` to coalesce and dispatch messages.
  3. `sim_loop()` to drive MD/Relax when `state.running` is true.
  Tasks share an `asyncio.Event` (`stop_evt`) so any failure cancels the session cleanly.
- **Message queue**: incoming frames are stored in `msg_buf`. Consecutive `USER_INTERACTION` (sparse) messages are merged by `_coalesce_user_interactions()` (last-write wins per index) before processing.

## Protocol handling
- **Transport**: protobuf-only WebSocket (`session.proto`). Text frames are ignored. ServerResult frames contain either `frame` (positions/forces/etc.) or `notice` messages like `WAITING_FOR_ACK`.
- **Sequencing**: every outbound frame increments `state.server_seq`. Client-provided `seq`/`ack` update `state.client_seq` and `state.client_ack`. Backpressure triggers when `server_seq - client_ack >= max_unacked` (default 10).
- **Correlation counters**: frontend `user_interaction_count` and `sim_step` are echoed back to support per-atom gating on the client. Simulation frames capture a UIC snapshot at step start.
- **Full snapshot flag**: `UserInteractionSparse.full_update` is a required field (explicit presence via proto oneof). Idle or running paths reject payloads that omit the flag with `PROTO_VIOLATION`.
- **UMA-only guard**: `_assert_uma()` rejects non-UMA calculators for streaming paths; MD/Relax calls fail fast with `CALCULATOR_NOT_SUPPORTED`.

## User interaction & idle compute flow
- `_apply_user_interaction()` applies sparse deltas immediately. Geometry edits while a simulation is running are rejected with `STRUCTURE_CHANGED`; idle sessions apply updates in place and flag `need_idle_emit`.
- `natoms` ships on every `UserInteractionSparse` so resize intent is explicit; `_ensure_arrays_sized()` grows/shrinks arrays accordingly.
- Full snapshots (`full_update=true`) resize to the requested atom count, zero velocities if none were provided, reset cached forces/sim-step counters, and leave the session paused until a new `START_SIMULATION` arrives. An idle recompute is queued right away.
- When the viewer is idle (`state.running == False`) and geometry is valid, `_handle_user_interaction()` runs `ASEWorker.run_simple` to compute forces/energy. Replies include forces (shape-checked), positions/velocities echo, and optional stress.
- If geometry is incomplete (e.g., only Z or only positions arrived), the server responds with `IDLE_WAITING_FOR_GEOMETRY` notices instead of errors.

## Simulation loop
- `sim_loop()` acquires a worker via `WorkerPool.any()` and calls `run_md` or `run_relax` for each iteration:
  - Requests include current atomic numbers, positions, velocities, cell, simulation parameters, and optional precomputed hints.
  - Worker responses are coerced to numpy arrays; velocities/forces shapes are validated.
  - Energies, optional kinetic/temperature, and stress (when provided) are forwarded to the client.
- Precomputed caching: after each frame, the loop captures `energy`, `forces`, and `stress` into a `PrecomputedValues` object so the next step can reuse them if the geometry is unchanged.
- Backpressure: when unacked frames exceed the threshold, the loop emits a `WAITING_FOR_ACK` notice (echoing counters) and stalls until a higher ack is observed. Upon ack, an `ACK_CLEARED` notice is optionally sent.

## Worker behaviours
- `ASEWorker.run_simple` → `SimpleIn` → `_simple_calculate()` to compute energy/forces (and stress when available).
- `run_relax` performs BFGS steps (`_relax_run`), applies precomputed values before the first energy evaluation, and returns updated positions/forces/stress.
- `run_md` sets or initializes velocities (respecting provided arrays), runs a Langevin integrator one step at a time, enforces displacement sanity checks, and reports kinetic/temperature metrics.
- All worker paths validate atomic numbers, positions, atom counts, and velocities; violations raise `HTTPException` propagated to the ingress.

## Deployment utilities & observability
- `deploy()` auto-detects GPU/CPU availability (override with `MLIPVIEW_FORCE_CPU`). It configures Serve logging to suppress access logs.
- `model_runtime.health_snapshot()` inspects UMA handle stats and reports device/capacity info for `/serve/health`.
- Optional debug flags (`WS_DEBUG`, `WS_LOG_CALLS`, UMA geom debug) add verbose logging for tracing sessions.

## Error handling & cleanup
- Compute failures in idle path send `COMPUTE_ERROR` notices with the exception string.
- Simulation loop exceptions log to stdout and tear down the session (`stop_evt.set()`).
- While a simulation is running, arrival of a `full_update` is rejected with a `STRUCTURE_CHANGED` notice (`simulation_stopped=True`); clients must pause before sending a full reset.
- `STOP_SIMULATION` triggers a `SIMULATION_STOPPED` notice and resets cached precomputed hints.
- WebSocket disconnects cancel all running tasks and release the worker loop.

## Observability updates
- `install_predict_handle()` now logs the device (`predict_unit_device=cuda|cpu`) so test runs confirm GPU usage when available.
- Proto generation is consolidated to a single `session_pb2.py`; stale nested copies were removed to avoid drift between server and client stubs.
