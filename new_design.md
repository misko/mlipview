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

| Type | Purpose | Notes |
|------|----------|-------|
| `USER_INTERACTION` | Initialization and subsequent updates | Sends atomic_numbers, positions, velocities?, cell?, and counters |
| `START_SIMULATION` | Begin MD or Relax loop | Includes params and counters |
| `STOP_SIMULATION` | End active simulation | — |
| `PING` | Acknowledge frames (ack=seq) | Maintains backpressure |

### 2.2 ServerResult Payload

| Context | Contents |
|----------|-----------|
| Idle (no simulation) | `forces`, `energy`, `stress?`; omit positions |
| Simulation frames | `positions`, `velocities?`, `forces`, `energy`, `cell?`, `stress?` |
| Informational | `message="WAITING_FOR_ACK"` |

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

| From | Event | To | Action |
|------|--------|----|--------|
| idle | `USER_INTERACTION(init)` | idle | Initialize session; compute one-shot forces/energy (and stress if cell); positions omitted |
| idle | `USER_INTERACTION(update)` | idle | Compute one-shot forces/energy (and stress if cell); positions omitted |
| idle | `START_SIMULATION` | running-{md,relax} | Begin sim loop |
| running | `USER_INTERACTION(update)` | running | Update state; next frame reflects change |
| running | `STOP_SIMULATION` | idle | Stop sim |
| any | Backpressure threshold exceeded | backpressured | Emit `WAITING_FOR_ACK`, pause until ack |
| any | `WebSocketDisconnect` | closed | Clean up |

**Emitted Frames**

- *Idle compute (includes init):* forces, energy (and stress if cell); positions omitted  
- *Simulation frame:* positions, velocities?, forces, energy, stress?  
- *Informational:* message = `WAITING_FOR_ACK`

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

| Frame Type | Behavior |
|-------------|-----------|
| `idle compute` | Apply per-atom gating for forces; update energy (and stress when present); positions omitted |
| `simulation` | Apply positions (gated), forces, energy |
| `WAITING_FOR_ACK` | Log only; ensure ack loop active |

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
