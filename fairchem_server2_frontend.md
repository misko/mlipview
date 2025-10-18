# Frontend migration to stateful WS backend (fairchem_local_server2)

This doc inventories current frontend API usage, reviews the new websocket server (v2), and proposes an incremental plan to switch from stateless REST to a stateful WS stream with interaction and step counters.

## Current frontend API calls (stateless)

Primary locations in `public/` and `scripts/` that call the backend via HTTP:

- public/index.js
	- simple forces/energy: POST `${base}/serve/simple` (resolved by `getEndpointSync('simple')`)
	- relax step: POST `${base}/serve/relax`
	- md step: POST `${base}/serve/md`
	- Notes: carries optional `cell`, `pbc`, `precomputed`, and for MD includes optional `velocities` from local cache. Includes rich client-side versioning and partial-apply safeguards.
- public/fairchem_provider.js
	- fairchemCalculate: POST `/simple_calculate` (standalone provider used by createFairChemForcefield). In the main app, simple/relax/md endpoints above are used; this file is a compatibility path and test harness.
- scripts/relax_via_api.js, scripts/md_via_api.js
	- CLI helpers that perform REST calls to check server behavior.

No direct WebSocket calls exist in the frontend yet. A placeholder `public/fairchem_ws_client.js` file exists (empty in repo snapshot).

## New WS server v2 overview

`fairchem_local_server2/ws_app.py` exposes:
- GET /serve/health
- WebSocket /ws carrying protobuf messages defined in `session.proto` and `session_pb2.py`.
- The server maintains per-connection session state and streams simulation frames continuously when running (MD/Relax, 1 step per tick).

Messages:
- ClientAction (INIT_SYSTEM, UPDATE_POSITIONS, START_SIMULATION, STOP_SIMULATION, PING) with optional `simulation_params`.
- ServerResult carrying positions, velocities, forces, optional cell and message.

Extension added in this change:
- ClientAction optional fields: `user_interaction_count`, `sim_step` (client-provided counters for correlation).
- ServerResult optional fields: `user_interaction_count`, `sim_step` (echoed and/or server-managed step index).
- Server now increments `sim_step` per produced frame and echoes the last seen `user_interaction_count`.

## Migration plan (frontend)

Goal: shift relax/MD from request/response (REST) to continuous streaming over WS, preserving existing user interaction semantics and partial-apply rules.

Key steps:
1) Introduce a WS client in `public/fairchem_ws_client.js` that:
	 - opens ws to `${base}/ws` (protocol chosen from page origin),
	 - encodes ClientAction protobuf messages,
	 - decodes ServerResult frames and emits them through an event API,
	 - carries and echoes `userInteractionVersion` and a `simStep` counter.
2) Integrate in `public/index.js`:
	 - map existing relax/md loops to start/stop the WS stream instead of polling REST,
	 - feed incoming frames into the same partial-apply path (respecting modified/dragged/latch exclusions),
	 - maintain RPS/latency metrics and energy plotting as before.
3) Keep REST simple/relax/md as fallback behind a feature flag for A/B and tests during transition.

Data contract (WS):
- Inputs: ClientAction with seq, optional ack, INIT with atomic_numbers, positions, velocities, cell. START_SIMULATION with params (md or relax). Optional user_interaction_count and sim_step.
- Outputs: ServerResult with seq, client_seq, positions, velocities, forces, optional cell, and the echoed `user_interaction_count` and server `sim_step`.
- Error modes: disconnect, timeouts, malformed frames. Client should auto-reconnect and re-INIT when appropriate.

Edge cases to handle:
- User editing during in-flight frames: use existing partial-apply exclusion sets (dragging, latched, modifiedSince).
- Reset/abort: bump `resetEpoch`, stop reading frames, clear caches, re-INIT on next start.
- Focus gating: pause consumption/start requests when window is unfocused (desktop) similar to current behavior.

## Test impact and plan

Existing browser unit tests rely on REST pacing and single-step semantics. Changes:
- Update or wrap tests that expect `fetch` calls for md/relax to either enable the REST fallback flag or run with WS using a small mock server.
- Add tests:
	- WS connection lifecycle (connect, INIT_SYSTEM, START_SIMULATION MD then STOP, frames received > 0).
	- Counter propagation: ensure `user_interaction_count` and `sim_step` are echoed in ServerResult.
	- Partial-apply with exclusions: when a user drag happens, frames do not overwrite dragged atoms.

Minimal concrete changes proposed in this PR:
- Extend protobuf and server to accept/echo counters (done).
- Add a basic `fairchem_ws_client.js` that can connect and expose a small API (next step).
- Keep `public/index.js` logic unchanged for now (still REST), then stage a follow-up to switch orchestration to use the WS client when a feature flag is on.

## Phased rollout
- Phase 1: server + proto counters; provide WS client library; add a couple of WS integration tests under `fairchem_local_server2/tests_py`.
- Phase 2: behind feature flag, route `startMDContinuous`/`startRelaxContinuous` to WS; REST remains fallback.
- Phase 3: remove REST paths and update scripts/ docs.

## Protobuf JS generation (browser)

- Install deps (already in package.json): google-protobuf, protoc-gen-js.
- Generate JS stubs for browser usage:
	- Run: `npm run gen:proto:js`
	- Output will be placed under `public/proto/session_pb.js` (CommonJS with binary wire format).
	- `public/fairchem_ws_client.js` imports this module and encodes/decodes ClientAction/ServerResult.
