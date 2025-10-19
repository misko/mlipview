# mlipview v1 WebSocket Architecture — Current Implementation

Date: 2025-10-19
Status: Draft (reflects code on branch `ws` at time of writing)

## Overview

mlipview renders real-time output from UMA Fairchem MLIP across desktop/VR/AR.
The legacy REST backend (`fairchem_local_server`) is being replaced by a
WebSocket + protobuf backend (`fairchem_local_server2`). This document captures
what the current code actually does end-to-end (frontend and backend).

## Actors and processes

- Frontend (public/)
  - `index.js`: orchestrates viewer init, molecule state, interactions, energy
    plotting, continuous MD/relax runs, and WS usage.
  - `fairchem_ws_client.js`: thin protobuf-over-WS client with helpers for
    sending actions and subscribing to results.
- Backend (fairchem_local_server2/)
  - `ws_app.py`: FastAPI endpoint `/ws` using Ray Serve; hosts a worker pool to
    run ASE + UMA/LJ computations. Only protobuf WS frames are supported.
  - `worker_pool.py`: executes `run_simple`, `run_md`, `run_relax` using UMA or
    LJ fallback.
  - `session_models.py`: in-memory session state for one WS connection.

## Protocol (protobuf)

Message types are defined in `fairchem_local_server2/session.proto`.

- ClientAction (current codebase)
  - Fields: seq, optional ack, atomic_numbers, positions, velocities, cell,
    optional simulation_type/params, optional counters (user_interaction_count,
    sim_step)
  - Type enum (current repo state): INIT_SYSTEM, USER_INTERACTION,
    START_SIMULATION, STOP_SIMULATION, PING, SIMPLE_CALCULATE
- ServerResult
  - Fields: seq, client_seq echo, positions, forces, velocities, cell,
    optional energy, optional message; counters echo (user_interaction_count,
    sim_step)

Notes:
- The frontend already treats WS as protobuf-only. The backend ignores incoming
  text messages, and sends binary protobuf frames.
- INIT_SYSTEM is still present in code paths but we are removing it (see changes
  below).
- SIMPLE_CALCULATE exists in the current code and is used both explicitly by the
  client and implicitly by the server on idle USER_INTERACTION. This duplication
  is slated for removal in the proposed design; only USER_INTERACTION and
  simulation frames will be used to deliver forces/energy.

## Backend behavior (ws_app.py)

- On connection: accept, create a new `SessionState`, spawn two async tasks:
  - sim_loop(): produces MD/relax frames when `state.running` is true.
  - recv_loop(): processes client actions, updates state, may trigger simple
    calculations when idle.

- Backpressure:
  - The server tracks `server_seq` vs `client_ack`. If `server_seq - client_ack`
    exceeds `max_unacked` (10), it stops producing simulation frames and instead
    occasionally sends a "WAITING_FOR_ACK" ServerResult (with last known data)
    until the client `ACK`s.

- Actions handled by recv_loop():
  - INIT_SYSTEM (legacy): sets atomic_numbers, positions, optional velocities,
    optional cell, clears forces, and sends a single frame marked `initialized`.
  - SIMPLE_CALCULATE: runs `run_simple` with current state (UMA, LJ fallback),
    returns forces and optional energy; does not increment `server_seq`.
    Note: slated for removal; idle USER_INTERACTION responses already cover
    this functionality.
  - USER_INTERACTION: Partial update handler. When `state.running` is false,
    runs a one-shot simple calculation and emits forces/energy update. When MD
    or relax is running, updates state and allows sim_loop to continue.
    - If atomic_numbers and positions are provided together, this is treated
      as an initialization (sets atoms/positions/cell/velocities).
  - START_SIMULATION: sets simulation type and params, flips `running=true`.
  - STOP_SIMULATION: flips `running=false`.
  - PING: may carry `ack`; updates client_ack.

- sim_loop():
  - While running, each step:
    - Calls `run_md` or `run_relax` for 1 step and updates positions/vel/forces.
    - Increments `server_seq` and server-side `sim_step` counter.
    - Sends ServerResult with current positions/vel/forces and optional energy.

- Correlation fields:
  - Both Action and Result carry optional `user_interaction_count` and `sim_step`.
    Today the server typically echoes the most recent value in state. In the
    proposed design, the server must echo the exact counter associated with the
    specific triggering request (even if newer counters arrive while computing
    a response).

## Frontend behavior

- On viewer init (`index.js:initNewViewer`):
  - Sets up molecule state, bonds/selection/manipulation services, scene & view.
  - Creates singleton WS client via `getWS()`.
  - Ensures connection and calls `__wsEnsureInit()` which currently sends
    INIT_SYSTEM; in our patch we switch to USER_INTERACTION for init.
  - Requests a baseline energy/forces by sending USER_INTERACTION and then an
    explicit SIMPLE_CALCULATE one-shot. Proposed design removes the explicit
    one-shot and relies only on the USER_INTERACTION idle response.

- WS client (`public/fairchem_ws_client.js`):
  - Handles connect/open/close/error.
  - Serializes ClientAction via Bufbuild schemas to protobuf binary.
  - Decodes ServerResult binary. Listener API `onResult(cb)` plugs into
    `index.js` streaming loops.
  - API surface:
    - initSystem({ atomic_numbers, positions, velocities?, cell? })
    - userInteraction({ positions?, velocities?, cell?, dragLockIndex? })
    - setCounters({ userInteractionCount, simStep })
    - startSimulation({ type: 'md'|'relax', params })
    - stopSimulation()
    - requestSimpleCalculate(): one-shot forces/energy
    - requestSingleStep({ type, params }): start -> first frame -> stop
    - ack(seq): send PING with ack for backpressure

- Interaction gating (frontend):
  - The viewer tracks `userInteractionVersion` and which atoms are being dragged
    or were modified during an operation. When a server result arrives, the
    viewer applies new positions partially: atoms being interacted with remain
    unchanged. Forces for those atoms are also merged from current state.
  - Discards stale WS frames when `userInteractionCount` on the server result is
    lower than the frontend’s current interaction counter.
  - The proposed design formalizes per-atom state client-side: each atom keeps
    a `last_user_interaction` counter. Server remains unaware of these per-atom
    values due to unknown client latency.

- Continuous streams:
  - Relax and MD `start*Continuous` use WS streaming. They set counters,
    subscribe to `onResult`, apply frames with partial updates (respecting drag
    and latching), update forces/energy, and auto-stop after N steps or when
    stopped by UI.

## Known gaps/quirks (current code)

- Duplicate init verbs: both INIT_SYSTEM and USER_INTERACTION exist; code
  supports both. This can be simplified.
- Server does not gate per-atom updates by interaction count; instead
  front-end handles partial applies and staleness discard. This matches desired
  product intent but leaves some policy on the client.
- Some long lines and debug prints were reflowed; linting should be clean now.
- Several REST code paths in `index.js` remain but are guarded to throw (REST
  removed). They can be deleted to reduce confusion.
- SIMPLE_CALCULATE duplication (client + server implicit) adds complexity and
  will be removed in the proposed design.
- Tests and scripts still reference INIT_SYSTEM; these will need updating after
  the protobuf change.

## Sequence diagrams (simplified)

1) Baseline forces on idle (current):
- Frontend connects WS -> sends USER_INTERACTION init -> may also call
  SIMPLE_CALCULATE
- Backend echoes initialized -> computes and returns forces+energy either
  implicitly (idle USER_INTERACTION) or via the one-shot
- Frontend updates forces, draws energy point

2) USER_INTERACTION while idle:
- Frontend sends USER_INTERACTION positions (debounced) -> may also call
  SIMPLE_CALCULATE -> Backend runs run_simple -> returns forces/energy

3) Start MD:
- Frontend startSimulation(MD) -> Backend running=true
- Backend sim_loop produces frames -> Frontend applies partial positions,
  merges forces, updates visuals and energy
- Frontend ACKs periodically to maintain flow

## Stability and performance

- Backpressure: server limits to 10 unacked frames; clients call ack() on each
  frame to continue streaming.
- RPS counter/metrics on frontend for runtime pacing and UI feedback.
- UMA errors fall back to LJ for simple_calculate; MD/relax use UMA when
  available.

## Current risks

- INIT_SYSTEM removal touches session.proto, generated Python/JS code, backend
  handler, frontend client and tests. Regeneration steps must be executed.
- Per-atom interaction count policy is only on the frontend; any divergence in
  future may require server awareness to prune updates earlier.

