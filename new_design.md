# mlipview v1 WebSocket Architecture — Proposed Design (v1)

**Date:** 2025-10-19  
**Status:** Proposal ready for implementation

---

## 1. Goals

- Single, protobuf-based WebSocket API
- Remove legacy actions (`INIT_SYSTEM`, `SIMPLE_CALCULATE`)
- Use `USER_INTERACTION` for both init and updates
- Use `START_SIMULATION` / `STOP_SIMULATION` for active MD/Relax streams
- Deterministic, low-latency force/energy updates while idle  
  (positions advance only during simulations)
- Client-side per-atom gating only (server unaware of per-atom state)
- Cleaner backend/frontend loops with fewer special cases

---

## 2. Protocol Overview

### 2.1 ClientAction Types

| Type               | Purpose                               | Notes                                                             |
| ------------------ | ------------------------------------- | ----------------------------------------------------------------- |
| `USER_INTERACTION` | Initialization and subsequent updates | Sends atomic_numbers, positions, velocities?, cell?, and counters |
| `START_SIMULATION` | Begin MD or Relax loop                | Includes params and counters                                      |
| `STOP_SIMULATION`  | End active simulation                 | —                                                                 |
| `PING`             | Acknowledge frames (ack=seq)          | Maintains backpressure                                            |

### 2.2 ServerResult Payload

| Context              | Contents                                                           |
| -------------------- | ------------------------------------------------------------------ |
| Idle (no simulation) | `forces`, `energy`, `stress?`; omit positions                      |
| Simulation frames    | `positions`, `velocities?`, `forces`, `energy`, `cell?`, `stress?` |
| Informational        | `message="WAITING_FOR_ACK"`                                        |

Note: `stress` is included when a `cell` is present and is normally a 3x3 tensor.

**Counters**

- `user_interaction_count` (UIC): Monotonic client counter attached to every action.
  - Idle/init responses: server echoes the exact UIC from the triggering request.
  - Simulation frames: server echoes the UIC snapshot captured at the start of the step.
- `sim_step`: Server-authoritative step index incremented per simulation frame; the client may also send a local sim_step for UI correlation

---

## 3. Message Contract

### 3.1 `ClientAction.USER_INTERACTION`

**Input**

```protobuf
{
  atomic_numbers?,
  positions?,
  velocities?,
  cell?,
  user_interaction_count,
  sim_step
}
```

**Behavior**

- If `atomic_numbers + positions` provided → initialize new session
- Otherwise → update existing session state (omit unchanged fields)

**Errors**

- Invalid array lengths
- Non-finite values
- Missing atoms on init

---

### 3.2 `ClientAction.START_SIMULATION`

**Input**

```protobuf
{
  simulation_type: MD | RELAX,
  simulation_params,
  counters
}
```

**Errors**

- Missing or invalid parameters

---

### 3.3 `ServerResult`

**Success**

```protobuf
{
  positions? (sim only),
  velocities? (sim only),
  forces?,
  energy?,
  cell?,
  stress?,
  message?,
  user_interaction_count,
  sim_step
}
```

**Informational**

```protobuf
{
  message = "WAITING_FOR_ACK"
}
```

---

## 4. Backend Behavior

- Remove `INIT_SYSTEM` and `SIMPLE_CALCULATE`; use protobuf-only WS
- **On `USER_INTERACTION`:**
  - **Init and idle:** if `running == false` → ensure session is initialized (when `atomic_numbers + positions` are provided), then compute one-shot `forces`/`energy` (include `stress` if `cell` present); positions are omitted.
  - **Running:** if `running == true` → update state and let `sim_loop` continue (no immediate frame to avoid duplicates)
- Maintain backpressure via ACK (`WAITING_FOR_ACK` informational frames)
- Emit lightweight `WAITING_FOR_ACK` frames (message + counters only) and pause stepping until ack advances
- Validate atom counts, matrix sizes; report errors in `ServerResult.message`
- No per-atom masking (client responsibility)

---

## 5. Frontend Behavior

1. **Initialization**
   - Open WS connection
   - Load default molecule (e.g. `roy.xyz`)

- Send initial `USER_INTERACTION` with `atomic_numbers`, `positions`, (and optional `cell`); expect an idle result with `forces`/`energy` (and `stress` when `cell` present); positions are not included.

2. **User Interaction (drag/rotate, etc.)**
   - Increment `user_interaction_count`
   - Send `USER_INTERACTION` with updated positions (omit unchanged fields)

3. **Simulation Control**
   - Send `START_SIMULATION` (`MD` or `RELAX`)
   - Subscribe to streaming `ServerResult` frames
   - Acknowledge each frame (`PING(ack=seq)`)

4. **Frame Handling**

- Never discard an entire frame solely due to a lower `userInteractionCount`.
- Apply updates per-atom:
  - For each atom i, apply incoming positions/forces only if `last_user_interaction[i] <= frame.userInteractionCount`.
  - While dragging, atoms with newer local interaction counters are excluded from position application on simulation frames and from force overrides on idle frames.
- Scalars (e.g., `energy`, and `stress` when present) may be updated from the frame regardless, since they’re global and not per-atom gated.

---

## 6. Counters and Gating

- **Global `user_interaction_count` (UIC):**  
  Attached to every outgoing action. Idle/init responses echo the exact request UIC; simulation frames echo the step-start UIC snapshot.
- **Per-atom `last_user_interaction`:**  
  Client-only mechanism preventing overwriting of actively edited atoms
- Server remains stateless regarding per-atom counts
- Client discards stale frames (based on UIC comparison)

---

## 7. Run Loops and State Machines

### 7.1 Backend WebSocket Loop

**States**

`idle` | `running-md` | `running-relax` | `backpressured` | `closed`

**Events**

`USER_INTERACTION(init/update)` | `START_SIMULATION` | `STOP_SIMULATION` | `PING(ack)` | `WebSocketDisconnect`

**Transitions**

| From    | Event                           | To                 | Action                                                                                     |
| ------- | ------------------------------- | ------------------ | ------------------------------------------------------------------------------------------ |
| idle    | `USER_INTERACTION(init)`        | idle               | Initialize session; compute one-shot forces/energy (and stress if cell); positions omitted |
| idle    | `USER_INTERACTION(update)`      | idle               | Compute one-shot forces/energy (and stress if cell); positions omitted                     |
| idle    | `START_SIMULATION`              | running-{md,relax} | Begin sim loop                                                                             |
| running | `USER_INTERACTION(update)`      | running            | Update state; next frame reflects change                                                   |
| running | `STOP_SIMULATION`               | idle               | Stop sim                                                                                   |
| any     | Backpressure threshold exceeded | backpressured      | Emit `WAITING_FOR_ACK`, pause until ack                                                    |
| any     | `WebSocketDisconnect`           | closed             | Clean up                                                                                   |

**Emitted Frames**

- _Idle compute (includes init):_ forces, energy (and stress if cell); positions omitted
- _Simulation frame:_ positions, velocities?, forces, energy, stress?
- _Informational:_ message = `WAITING_FOR_ACK`

---

### 7.2 Frontend WebSocket Loop

**States**

- `disconnected`, `connecting`, `connected-idle`, `connected-running-md`, `connected-running-relax`, `closed`
- `dragging` is a transient UI sub-state (orthogonal to the WS state)

Allowed transitions among connected states

- From any of: `connected-idle`, `connected-running-md`, `connected-running-relax`
  - To any of: `connected-idle`, `connected-running-md`, `connected-running-relax`
  - Typical flows:
    - `connected-idle` → `connected-running-md` on START_SIMULATION(MD)
    - `connected-idle` → `connected-running-relax` on START_SIMULATION(RELAX)
    - `connected-running-md` → `connected-idle` on STOP_SIMULATION
    - `connected-running-relax` → `connected-idle` on STOP_SIMULATION
    - `connected-running-md` → `connected-running-relax`: implement as STOP then START(RELAX)
    - `connected-running-relax` → `connected-running-md`: implement as STOP then START(MD)

**Actions**

- `USER_INTERACTION` (init/update)
- `START_SIMULATION` / `STOP_SIMULATION`
- `PING(ack)`

**Incoming Frames**

| Frame Type        | Behavior                                                                                     |
| ----------------- | -------------------------------------------------------------------------------------------- |
| `idle compute`    | Apply per-atom gating for forces; update energy (and stress when present); positions omitted |
| `simulation`      | Apply positions (gated), forces, energy                                                      |
| `WAITING_FOR_ACK` | Log only; ensure ack loop active                                                             |

**Counter Semantics (client)**

- Maintain monotonic `UIC`; discard any frame with smaller value
- Track per-atom last interaction markers for gating
- Maintain optional local `sim_step` for UI; prefer server’s authoritative value

---

## 8. Migration Plan

1. **Protocol**
   - Remove `INIT_SYSTEM` and `SIMPLE_CALCULATE` from `session.proto`
   - Regenerate Python/JS protobuf stubs

2. **Backend**
   - Unify init/update logic under `USER_INTERACTION`
   - Add `stress` field to `ServerResult`
   - Implement input validation helpers

3. **Frontend**

- Switch all flows to WS (`USER_INTERACTION`, `START/STOP_SIMULATION`)
- Remove REST and `SIMPLE_CALCULATE` references
- Update initialization handling to consume the first idle compute result instead of a special init echo

4. **Tests**
   - Update integration tests and scripts
   - Validate UIC echo and frame gating

---

## 9. Known Issues / Open Items

- **Stress field:** Add optional `stress` (normally 3x3) to `ServerResult`
- **Energy semantics:** Currently labeled `energy`; can be clarified later
- **Backpressure frames:** `WAITING_FOR_ACK` frames should be lightweight (message + counters); revisit further reductions if needed
- **Mid-step updates:** Reflected only in next frame (by design)
- **Cell updates:** Allowed; verify UMA/ASE stability under changing cell
- **Legacy JSON:** Some dev paths still support JSON; prune post-migration

---

## 10. Acceptance Criteria

- Default molecule loads and WS connection established
- `USER_INTERACTION` init returns an idle compute result (`forces`/`energy`/`stress?`); positions are omitted
- Idle updates produce deterministic force/energy results
- Simulations advance positions only when running
- Frames during drags do not overwrite active atoms (per-atom gating verified)
- UIC gating verified: frames are applied per-atom only when `last_user_interaction[i] <= frame.userInteractionCount`
- No `INIT_SYSTEM` or `SIMPLE_CALCULATE` remain in protocol or code
- All tests pass

---

## 11. Implementation TODOs and Task Breakdown

This section enumerates concrete engineering tasks to deliver the proposed design. Grouped by area with suggested file paths and acceptance checks.

### 11.1 Protocol and Protobuf

- Remove legacy action from protocol
  - [ ] Remove `SIMPLE_CALCULATE` from `ClientAction.Type` in `fairchem_local_server2/session.proto`
  - [ ] Ensure comment states: “INIT_SYSTEM and SIMPLE_CALCULATE removed; use USER_INTERACTION for init and idle compute”
- Add stress output
  - [ ] Add `optional Mat3 stress = 15;` to `ServerResult` (3x3, row-major)
  - [ ] Document omission rules: positions omitted on idle responses; stress included when `cell` present or when available
- Regenerate stubs
  - [ ] Regenerate Python stubs used by backend
  - [ ] Regenerate JS/TS stubs used by frontend (Buf.build or protoc pipeline used by `public/`)
  - [ ] Update any import paths in `public/fairchem_ws_client.js` and `fairchem_local_server2/*.py`
- Contract validation
  - [ ] Verify binary WS frames only; text frames rejected with clear error
  - [ ] Versioning note: bump minor schema version in comments

Acceptance: repo builds, unit tests compile with new stubs, no references to `SIMPLE_CALCULATE` remain.

### 11.2 Backend (Python, `fairchem_local_server2/`)

- Unify init/update under `USER_INTERACTION`
  - [ ] In `ws_app.py` recv loop, treat `atomic_numbers + positions` as init
  - [ ] When `running == false`, compute one-shot forces/energy (and stress) and emit an idle frame with positions omitted
  - [ ] When `running == true`, update state only; no immediate extra frame
- Remove SIMPLE_CALCULATE
  - [ ] Delete/disable any `SIMPLE_CALCULATE` handling paths
  - [ ] Ensure idle `USER_INTERACTION` covers the previous one-shot behavior
- Counters and echoing
  - [ ] Capture and echo the exact `user_interaction_count` from the triggering action in idle responses
  - [ ] For simulation frames, echo the snapshot UIC captured at step start
  - [ ] Maintain and emit server-side `sim_step` (monotonic per stream)
- Backpressure
  - [ ] Keep `max_unacked` threshold; on exceed, pause stepping and emit light `WAITING_FOR_ACK` frames (message + counters only)
  - [ ] Resume when `ack` advances
- Validation and error handling
  - [ ] Validate array sizes (N atoms across all vectors), matrix size (3x3), finite values
  - [ ] Return error via `ServerResult.message` and do not crash the connection
- UMA/LJ compute
  - [ ] Ensure `run_simple` returns forces, energy, and stress (when feasible) consistently for idle responses
  - [ ] Ensure `run_md`/`run_relax` return positions, velocities?, forces, energy, stress?
- Logging/metrics
  - [ ] Structured logs for connect/disconnect, start/stop, backpressure, errors

Files: `fairchem_local_server2/ws_app.py`, `worker_pool.py`, `session_models.py`, `services.py`, `models.py`.

Acceptance: manual smoke via local run; unit tests (Python and JS) pass; Playwright e2e init/sim flows green.

### 11.3 Frontend (JS, `public/`)

- WS client API cleanup (`public/fairchem_ws_client.js`)
  - [ ] Remove `requestSimpleCalculate()` and references
  - [ ] Ensure `userInteraction()` handles both init and updates
  - [ ] ACK loop: send `PING(ack=seq)` for every frame; ensure `WAITING_FOR_ACK` frames trigger ACK loop if needed
  - [ ] Decode idle frames: positions omitted; forces/energy (and stress) present
  - [ ] Expose `onResult(cb)` semantics unchanged
- Viewer orchestration (`public/index.js` and related services)
  - [ ] Switch initialization to a single `USER_INTERACTION` and wait for idle compute
  - [ ] Remove any REST remnants and guards
  - [ ] Per-atom gating only: do not discard entire frames solely due to lower global UIC; instead gate per-atom by `last_user_interaction[i]`
  - [ ] While dragging, suppress position application for edited atoms on sim frames; merge forces appropriately on idle frames
  - [ ] Handle optional `stress` (store/update; UI may ignore or display in dev panel)
- Cleanups
  - [ ] Remove dead code for `INIT_SYSTEM`/`SIMPLE_CALCULATE` and REST
  - [ ] Update energy plotting and any consumers for idle frame semantics

Files: `public/index.js`, `public/fairchem_ws_client.js`, `public/fairchem_provider.js`, `public/core/*` (where state apply lives), selection/interaction services.

Acceptance: existing UI interactions work; energy updates appear after init and edits; MD/Relax run and stop as expected.

### 11.4 Tests and Tooling

- Lint/type/test harness
  - [ ] Ensure ESLint/Prettier pass
  - [ ] Verify Jest config runs in jsdom and Node where needed
  - [ ] Playwright config runs headed/headless
- CI
  - [ ] Add jobs for unit + e2e (optional if using local scripts initially)

Acceptance: CI green or local equivalent run green.

### 11.5 Migration and Cleanup

- [ ] Delete REST-only code paths and deprecated flags
- [ ] Update README and internal docs where protocol is described
- [ ] Remove any remaining references to `SIMPLE_CALCULATE`/`INIT_SYSTEM`

Acceptance: workspace grep returns zero matches for removed verbs.

---

## 12. Test Plan (Jest mock browser + Playwright)

The following specifies concrete tests to implement. File paths are suggestions aligned with repo structure.

### 12.1 Unit Tests (Jest, jsdom + Node)

Mocking strategy

- Use jsdom for DOM-dependent modules (viewer integration)
- Mock WebSocket at the transport layer to inject binary protobuf frames
- Use lightweight fake clock where timing matters

Proposed test files (create under `tests/`)

- `tests/wsProtocol.initIdleCompute.spec.js`
  - Should send `USER_INTERACTION` with atoms+positions, receive idle frame without positions, with forces/energy (stress when cell present)
  - Should increment and echo `user_interaction_count`
- `tests/wsProtocol.noSimpleCalculate.spec.js`
  - Ensure client does not call any removed `requestSimpleCalculate()` and that init/update flows still yield forces/energy
- `tests/wsFrameGating.perAtom.spec.js`
  - Given per-atom `last_user_interaction` map, apply simulation frames and verify only non-active atoms update positions
  - Ensure scalar `energy` updates regardless of per-atom gating
- `tests/wsCounters.echoSemantics.spec.js`
  - Idle responses echo exact UIC of triggering action
  - Simulation frames echo UIC snapshot captured at step start
- `tests/wsBackpressure.ackLoop.spec.js`
  - Simulate server sending frames; verify client sends ACK promptly
  - Simulate `WAITING_FOR_ACK`; verify client resumes ACK and frames apply
- `tests/wsIdleUpdate.positionsOmitted.spec.js`
  - Verify that idle responses do not contain positions and client does not try to apply any
- `tests/wsCellAndStress.spec.js`
  - When cell is present, verify stress parsed and stored; absence tolerated
- `tests/wsValidation.errors.spec.js`
  - Malformed incoming frames (missing required fields) do not crash; error message surfaced/logged

Notes

- Provide binary payloads via small helper encoders using generated JS protobufs so tests remain schema-stable.

### 12.2 Node Integration Tests (optional, fast WS loop)

If practical, add a thin fake WS server:

- `tests/wsFakeServer.integration.spec.js`
  - Spins up a minimal WS endpoint that echoes protobuf `ServerResult`
  - Verifies end-to-end encoding/decoding and ACK loop without launching Python backend

### 12.3 Playwright E2E Scenarios (under `tests-e2e/`)

- `tests-e2e/ws.init-and-idle.spec.ts`
  - Load app, connect to WS (mock or local server), send initial `USER_INTERACTION`
  - Expect energy number appears; ensure no positions jump until simulation starts
- `tests-e2e/ws.start-md-stream.spec.ts`
  - Start MD; verify positions update over multiple frames; ACKs sent (instrumented via test hook or server logs)
- `tests-e2e/ws.drag-gating.spec.ts`
  - Begin MD; start dragging one atom; verify dragged atom position is not overwritten by incoming frames while others update
- `tests-e2e/ws.stop-simulation.spec.ts`
  - Start then stop; positions stop advancing; idle user interaction still returns forces without positions
- `tests-e2e/ws.backpressure.spec.ts`
  - Simulate client not ACKing, server emits `WAITING_FOR_ACK`; after resuming ACKs, frames continue
- `tests-e2e/ws.cell-stress.spec.ts`
  - Load system with cell; verify stress read and optionally displayed/logged
- `tests-e2e/ws.reconnect.spec.ts`
  - Force server disconnect; client reconnects and re-initializes via `USER_INTERACTION`

Implementation hints

- Prefer Playwright’s request interception or a test-dedicated backend profile to simulate backpressure and inject messages deterministically.

### 12.4 Coverage Matrix

- Protocol
  - No references to removed verbs; idle frames omit positions
  - Counters echo semantics covered (idle and sim)
- Frontend
  - Per-atom gating validated; no frame-wide discards
  - ACK loop resiliency
- Backend
  - Backpressure path and error validation (covered via Playwright where feasible)

---

## 13. Milestones and Ownership

1. Protocol update and stub regen (Day 0-1)

- Owner: Backend + Frontend shared
- Output: Updated `session.proto`, regenerated stubs, zero references to removed verbs

2. Backend refactor (Day 2-4)

- Owner: Backend
- Output: Unified `USER_INTERACTION`, idle compute semantics, backpressure tidy

3. Frontend refactor (Day 2-5)

- Owner: Frontend
- Output: API cleanup, init flow, per-atom gating finalized, UI uses idle compute

4. Tests (Day 3-6)

- Owner: QA/Dev
- Output: Jest unit tests + Playwright scenarios implemented and passing

5. Polish and cleanup (Day 6-7)

- Owner: Shared
- Output: Docs updated, CI green, release candidate tag

---

## 14. Risks and Mitigations

- Schema regen drift across FE/BE
  - Mitigation: commit generated artifacts or lock versions; add a script and document steps
- Per-atom gating inconsistencies
  - Mitigation: centralize gating in a single function with unit tests; use golden-frame fixtures
- Backpressure deadlocks
  - Mitigation: watchdog ACK loop; Playwright test to simulate and recover
- UMA stress availability
  - Mitigation: treat stress as optional; provide LJ fallback or omit with a clear message

---

## 15. Progress Log (2025-10-19)

- Protocol
  - Updated `fairchem_local_server2/session.proto`:
    - Removed `SIMPLE_CALCULATE` from `ClientAction.Type`
    - Added zero-valued `TYPE_UNSPECIFIED = 0` (proto3 requirement)
    - Added optional `stress` to `ServerResult`
    - Clarified idle compute semantics in comments
  - Regenerated JS ESM stubs via npm script. Python `session_pb2.py` regenerated via `protoc --python_out`.
  - Issue: `protoc` initially failed due to non-zero first enum; fixed by adding `TYPE_UNSPECIFIED`.

- Frontend
  - `public/fairchem_ws_client.js`:
    - Removed `requestSimpleCalculate()` API and all SIMPLE_CALCULATE usage
    - Decode optional `stress` from ServerResult
    - `waitForEnergy()` now resolves on energy-bearing frames only
  - `public/index.js`:
    - Replaced explicit one-shot SIMPLE_CALCULATE calls with `USER_INTERACTION` + `waitForEnergy()`
  - `public/fairchem_provider.js`: error guidance updated to recommend `userInteraction + waitForEnergy`

- Backend (in-progress)
  - Updating `ws_app.py` to remove SIMPLE_CALCULATE handling and unify idle compute under USER_INTERACTION
  - Added optional `stress` plumb in `_send_result_bytes`
  - Issues encountered:
    - Syntax errors during large refactor (stray else/elif, malformed f-strings)
    - Fixed major structural issues; remaining lint warnings about long lines will be cleaned in follow-up pass

Next steps

- Finish ws_app.py cleanup: remove legacy comments, shorten lines, ensure idle compute emits positions omitted and message unset
- Run unit tests and quick smoke to validate WS path

- Tests (Python)
  - Migrated remaining WS tests to USER_INTERACTION for initialization:
    - Updated `test_ws_md_0k.py`, `test_ws_md_parity.py`, `test_ws_roy_benchmark.py`
    - Replaced `INIT_SYSTEM` with `USER_INTERACTION` and adjusted comments
  - Fixed style/lint issues: long lines, slice spacing, and unused variables
  - Verified no syntax errors remain in the updated tests

- Additional updates
  - Stubs
    - Regenerated Python protobuf stubs via `protoc --python_out` (successful)
    - JS ESM stubs already regenerated and integrated in the frontend
  - Backend
    - `ws_app.py`: idle compute frames increment `server_seq` and include optional `stress` when available; backpressure/ACK semantics unchanged
  - Documentation
    - Added `testing.md` with updated WS testing patterns (waitForEnergy, requestSingleStep, setTestHook, injectTestResult)
    - Created `test-results/test_inventory.md` with a one-line-per-test inventory grouped by folder
  - Tests (Python)
    - Migrated `test_ws_simple_calc.py` and `test_ws_counters_optional.py` to USER_INTERACTION init
    - Confirmed lint/style pass after edits

Next steps (delta)

- Add a short “Test inventory” reference in `testing.md` pointing to `test-results/test_inventory.md`
- Continue migrating any JS/Playwright specs that still reference INIT_SYSTEM/SIMPLE_CALCULATE
- Optional cleanup: remove legacy generated artifacts (e.g., old JS proto `clientaction.js`, nested outdated `session_pb2.py`), and wrap remaining long lines in `ws_app.py`
