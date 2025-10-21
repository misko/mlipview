# Protobuf migration guide (flat arrays + oneof payloads)

This document explains the protocol changes introduced in the new `fairchem_local_server2/session.proto`, how to migrate client and server code, and how to reason about simulation lifecycle and counters under coalescing.

No legacy interface is kept; migrate fully to the new schema.

## What changed (at a glance)

- Replaced Vec3 messages with flat, packed `repeated double` arrays:
  - Positions, velocities, forces are now flattened `(N,3)` arrays: `[x0,y0,z0, x1,y1,z1, ...]`.
  - Fewer allocations, simpler parsing, faster transfer.
- Introduced `oneof` payloads for both directions:
  - ClientAction.oneof payload: `user_interaction`, `start`, `stop`, `ping`.
  - ServerResult.oneof payload: `frame`, `notice`.
- All numeric fields are `double` for consistency.
- Added `schema_version` (field id = 100) to both ClientAction and ServerResult.
- Kept `Mat3` for 3x3 matrices (row‑major: 9 doubles).
- Added correlation counters: `user_interaction_count` and `sim_step`.
- Idle results are delivered as `ServerResult.frame` with `energy` present and `positions` omitted.

## New schema (summary)

File: `fairchem_local_server2/session.proto`

- `ClientAction` (from client to server)
  - `seq` (uint64): monotonic client sequence
  - `ack` (uint64, optional): ACK a `ServerResult.seq`
  - `user_interaction_count` (uint64, optional): client-side monotonic counter of UI changes
  - `sim_step` (uint64, optional): optional client hint/correlation
  - `oneof payload`:
    - `user_interaction`:
      - `atomic_numbers: repeated int32`
      - `positions: repeated double` (flat packed, length % 3 == 0)
      - `velocities: repeated double` (flat packed, optional)
      - `cell: Mat3` (row-major, optional)
    - `start`:
      - `simulation_type: enum { MD, RELAX }`
      - `simulation_params: SimulationParams` (double fields)
    - `stop`: signal to end currently running simulation
    - `ping`: used mostly for ACK/backpressure
  - `schema_version: uint32 = 100`

- `ServerResult` (from server to client)
  - `seq` (uint64): server sequence number for this result
  - `client_seq` (uint64, optional): last seen client seq
  - `user_interaction_count` (uint64, optional): echo of latest UI counter
  - `sim_step` (uint64, optional): increments per produced frame during sim
  - `oneof payload`:
    - `frame`:
      - `positions, velocities, forces: repeated double` (flat)
      - `cell, stress: Mat3` (row‑major, optional)
      - `energy: double` (optional)
    - `notice`:
      - `message: string` (e.g., "SIMULATION_STOPPED")
      - `simulation_stopped: bool` (optional)
  - `schema_version: uint32 = 100`

## Migration steps

### 1) Regenerate code from the new schema

- Replace any older session.proto in your build with `fairchem_local_server2/session.proto`.
- Re-run `protoc` for your languages with `-I` pointing to repo root or the folder with the proto file.

Example for Python:
- Input proto: `fairchem_local_server2/session.proto`
- Output: `fairchem_local_server2/session_pb2.py`

### 2) Replace Vec3 with flat arrays

Old style (pseudocode):
- positions: repeated Vec3 { v: [x, y, z] }

New style:
- positions: repeated double (packed), in row-major triples: `[x0,y0,z0, x1,y1,z1, ...]`

Client-side construction:
- When sending `ClientAction.user_interaction`, fill `atomic_numbers`, then do:
  - for each atom `p = [x,y,z]`, push `x,y,z` into `positions`.
  - optionally do the same for `velocities` if needed.

Server-side/Client-side parsing:
- Convert to `(N,3)` by reshaping or chunking triples.
- Ensure `len(array) % 3 == 0`.

### 3) Switch to oneof payloads

- ClientAction: choose exactly one of `user_interaction`, `start`, `stop`, `ping`.
- ServerResult: check `WhichOneof("payload")` (or language equivalent):
  - `frame` for numeric results.
  - `notice` for text/flags events (including stop confirmations).

Remove any code using legacy type enums (e.g., TYPE_USER_INTERACTION, etc.).

### 4) Set and check schema_version

- Set `ClientAction.schema_version = 1` (or your negotiated version) on every message.
- Optionally validate `ServerResult.schema_version` from server.

### 5) Counters and lifecycle with coalescing

- `user_interaction_count` (UIC):
  - Send a monotonically increasing integer with each UI update.
  - The server echoes the latest seen value in `ServerResult.user_interaction_count`.
  - Due to coalescing, you may receive fewer frames than UI updates.
  - To know whether the latest UI was processed, compare `echoed_uic == last_sent_uic`.

- `sim_step`:
  - In `ServerResult`, increments per produced frame during a running simulation.
  - May start at 0 depending on implementation. Don’t assume 1-based.

- Detect simulation start:
  - Send `ClientAction.start` with `Start.SimType` and params.
  - Simulation is considered started when you receive the first `ServerResult.frame` that includes positions (non-empty positions array) for that run.
  - You can optionally gate on `sim_step` >= 0 and/or updated `client_seq` echo.

- Detect simulation stop:
  - Send `ClientAction.stop`.
  - Wait for a `ServerResult.notice` with either `notice.message == "SIMULATION_STOPPED"` or `notice.simulation_stopped == true`.
  - After stop, you may still see any in-flight frames due to buffering; only treat the `notice` as authoritative.

### 6) Idle results vs simulation frames

- Idle results (when no simulation is running):
  - Delivered as `ServerResult.frame` where `positions` are omitted (length 0) but `energy` and sometimes `forces` are present.
  - Use this for quick feedback on “dragging atoms” updates.

- Simulation frames:
  - Delivered as `ServerResult.frame` with positions, velocities, forces (non-empty positions array).

### 7) ACK/backpressure (ping oneof)

- The server applies backpressure and may not send more frames until receiving client ACKs.
- To ACK a frame:
  - Send a `ClientAction` with `ack = res.seq` and set the oneof payload to `ping`.
- ACK after each accepted frame to sustain throughput.

### 8) Matrices (Mat3)

- Kept as `Mat3 { repeated double m }` with length 9.
- Row-major order: `m = [r00, r01, r02, r10, r11, r12, r20, r21, r22]`.

## Code migration examples

### Client: initialize system (UI) and drain idle result

Pseudocode (Python protobuf):

```python
init = pb.ClientAction()
init.seq = next_seq()
init.schema_version = 1
ui = pb.ClientAction.UserInteraction()
ui.atomic_numbers.extend(Z)
for p in R:
    ui.positions.extend([float(p[0]), float(p[1]), float(p[2])])
init.user_interaction.CopyFrom(ui)
ws.send(init.SerializeToString())

# Drain one idle energy frame and ACK
while True:
    data = ws.recv()
    res = pb.ServerResult(); res.ParseFromString(data)
    if res.WhichOneof("payload") != "frame":
        continue
    fr = res.frame
    if getattr(fr, "energy", None) is not None and len(fr.positions) == 0:
        ack = pb.ClientAction(); ack.seq = next_seq(); ack.ack = int(res.seq)
        ack.ping.CopyFrom(pb.ClientAction.Ping())
        ws.send(ack.SerializeToString())
        break
```

### Client: start simulation, stream frames, and stop

```python
# Start MD
start = pb.ClientAction(); start.seq = next_seq(); start.schema_version = 1
st = pb.ClientAction.Start(); st.simulation_type = pb.ClientAction.Start.SimType.MD
sp = pb.SimulationParams(); sp.calculator = "uma"; sp.temperature = 300.0; sp.timestep_fs = 1.0; sp.friction = 0.02
st.simulation_params.CopyFrom(sp)
start.start.CopyFrom(st)
ws.send(start.SerializeToString())

# Read frames
while True:
    data = ws.recv()
    res = pb.ServerResult(); res.ParseFromString(data)
    if res.WhichOneof("payload") != "frame":
        continue
    fr = res.frame
    if len(fr.positions) == 0:
        # ignore idle frames during start-up
        continue
    # process positions/velocities/forces
    # ACK
    ack = pb.ClientAction(); ack.seq = next_seq(); ack.ack = int(res.seq)
    ack.ping.CopyFrom(pb.ClientAction.Ping())
    ws.send(ack.SerializeToString())
    if done:
        break

# Stop
stop = pb.ClientAction(); stop.seq = next_seq()
stop.stop.CopyFrom(pb.ClientAction.Stop())
ws.send(stop.SerializeToString())

# Await stop notice
while True:
    data = ws.recv(); res = pb.ServerResult(); res.ParseFromString(data)
    if res.WhichOneof("payload") == "notice":
        n = res.notice
        if (getattr(n, "simulation_stopped", False) is True) or (n.message == "SIMULATION_STOPPED"):
            break
```

### Counting UIC (robust to coalescing)

```python
# send UI updates
for i, R in enumerate(moves, start=1):
    ui = pb.ClientAction.UserInteraction();
    for p in R: ui.positions.extend([p[0], p[1], p[2]])
    msg = pb.ClientAction(); msg.seq = next_seq(); msg.schema_version = 1
    msg.user_interaction.CopyFrom(ui)
    msg.user_interaction_count = i
    ws.send(msg.SerializeToString())

# Later, drain frames and check echoed UIC
last_sent = i
while True:
    data = ws.recv(); res = pb.ServerResult(); res.ParseFromString(data)
    if res.WhichOneof("payload") != "frame": continue
    echoed = int(getattr(res, "user_interaction_count", 0))
    if echoed == last_sent: break
    # ACK to advance
    ack = pb.ClientAction(); ack.seq = next_seq(); ack.ack = int(res.seq)
    ack.ping.CopyFrom(pb.ClientAction.Ping()); ws.send(ack.SerializeToString())
```

## Testing and validation

- Unit tests updated under `fairchem_local_server2/tests_py` demonstrate usage patterns:
  - Idle compute: `test_ws_simple_calc.py`.
  - MD parity and 0K operation: `test_ws_md_parity.py`, `test_ws_md_0k.py`.
  - Start/Stop & idle behavior: `test_ws_start_stop_idle_relax.py`.
  - ROY frames benchmark: `test_ws_roy_benchmark.py`.
  - Coalescing behavior: `test_ws_coalescing.py` (rapid UIs coalesced into fewer idle frames; last UIC echoed).

## Do’s and don’ts

- Do: always set `schema_version` on client messages; optionally validate it on server results.
- Do: ACK frames using `ping` with `ack=res.seq` to avoid backpressure stalls.
- Do: treat `notice` as the authoritative simulation-stop signal.
- Don’t: rely on 1:1 mapping between UI updates and frames; use UIC echo to detect progress.
- Don’t: assume `sim_step` starts at 1; it may start at 0.
- Don’t: use the legacy Vec3 or type enum fields; they are removed.

## Troubleshooting

- No frames arriving after start:
  - Ensure you’re ACKing with `ping`.
  - Check that you’re ignoring idle frames (positions empty) and waiting for sim frames.
- Latest UI not reflected:
  - Compare echoed `user_interaction_count` to your last sent value.
  - If behind, continue ACKing frames to advance.
- Protobuf mismatches:
  - Confirm you regenerated stubs against the updated `session.proto` and your code imports the correct file.
