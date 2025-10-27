# Frontend Architecture (current)

## Entry point & initialization
- `public/index.js` exports `initNewViewer(canvas, initialData)`, which bootstraps the Babylon.js scene, Molecule state, services, and WebSocket wiring.
- On load it normalizes runtime configuration (drag throttles, min step interval, temperature defaults) via `window.__MLIP_CONFIG`.
- `createMoleculeState()` (domain/moleculeState.js) builds the authoritative state object holding elements, position vectors, cell metadata, event bus, and version counters.
- Initial geometry is clamped into a "safe sphere" if configured, baseline energy/forces are fetched, and optional initial cell snapshots are stored for reset.

## Domain services & event flow
- **Event bus**: the Molecule state exposes `bus` with `positionsChanged`, `forcesChanged`, `bondsChanged`, `selectionChanged`, etc. All render/UI modules subscribe to these events.
- **Bond/selection/manipulation**:
  - `createBondService()` recomputes bonds based on current positions, handling periodic wrapping and opacity heuristics.
  - `createSelectionService()` maintains `state.selection` and emits selection events.
  - `createManipulationService()` handles drag and bond-rotation transforms, updating positions and emitting bus events.
- **Picking & UI**: `createPickingService()` orchestrates Babylon picking, forwards drag gestures to the manipulation wrapper, and invokes energy recomputes. Desktop and VR interactions layer on top via `ui/touchControls.js`, `vr/setup.js`, and `vr/vr-picker.js`.
- **Rendering**: `createScene()` initializes Babylon engine/camera while `createMoleculeView()` builds thin-instance meshes (atoms, bonds, forces, ghosts) that react to bus events.

## WebSocket client
- `public/fairchem_ws_client.js` provides a singleton (via `getWS()`) that manages protobuf messages (`session_pb.js` schema), reconnect logic, and frame listeners.
- Functions exposed: `ensureConnected`, `userInteraction` (sparse or dense payloads), `startSimulation`, `stopSimulation`, `requestSingleStep`, `waitForEnergy`, `ack(seq)`, `setCounters`, etc.
- `userInteraction()` always attaches a `fullUpdate` boolean flag so the backend can enforce dense snapshot semantics; callers pass `{ full_update: true }` whenever natoms/atomic numbers/positions are sent as a full frame.
- Incoming frames are decoded to plain objects with `positions`, `forces`, `velocities`, `energy`, `stress`, `userInteractionCount`, and `simStep`. Every frame triggers auto-ACK (batched ping) and fan-out to registered listeners.
- Reconnect logic emits state events (`connecting`, `open`, `close`, `error`, `reconnect-scheduled`). UI hooks in `index.js` display a banner with countdown and “Reconnect now” button.

## Viewer orchestrator (index.js)
- Maintains counters: `userInteractionVersion`, `totalInteractionVersion`, `structureVersion`, and `lastAppliedUIC`. These gate outgoing edits and incoming frames to prevent stale data during drags.
- `modifiedByVersion`, `draggingAtoms`, `rotatingAtoms`, and `latchedUntil` track per-atom edits so the viewer can exclude those indices when applying server frames (per-atom gating).
- `buildExcludeSet()` composes the active drag/latch sets; `applyTriples()` applies coordinates unless the atom is excluded.
- `bumpUser()` and `bumpSim()` mutate version counters and invalidate the cached forces (`state.forceCache`) so idle recompute requests know when to refresh.
- `applyFullSnapshot()` is the single entry point for dense geometry changes (file load, add/remove atoms). It rewrites state, resets interaction caches, and reuses `ensureWsInit()` so the next upload is a dense `full_update` snapshot.
- `emitDuringDrag()` throttles sparse `USER_INTERACTION` updates to the backend while dragging (indices + triple positions only). Bond rotations queue sparse updates for the affected atoms with optional latching windows.
- Energy plot (`energyPlot` closure) records energies per step for UI display; `noteReqDone()` tracks loop rate (RPS label).

## Idle vs simulation flows
- **Initialization / idle**:
  - `ensureWsInit()` ensures the WS connection is open and, when geometry changed, sends a dense `userInteraction` snapshot (natoms + atomic numbers + positions [+ velocities/cell]) with `full_update: true`. It waits for `waitForEnergy()` to resolve so force cache is primed.
  - `ff.computeForces()` delegates to the WS client’s single-step helper for one-shot `run_simple` computations on the backend, updating energy/force cache.
- **Manual steps**: `doOneStepViaWS(kind)` ensures a position upload occurs (including velocities when fresh), calls `requestSingleStep`, applies positions from the response, and updates dynamics metadata.
- **Continuous streaming**:
- `startContinuous()` sets `running.kind` (`md` or `relax`), sends `startSimulation`, and subscribes to WS frames.
- `handleStreamFrame(kind, ws, frame)` applies positions (respecting per-atom gating), updates forces/energy/temperature, pushes energy to the plot, acks the sequence, and updates `lastAppliedUIC`.
- Timeline playback reuses the same path via `applyTimelineFrame()`, which honours the stored `cell.enabled` flag so scrubbing through DVR history never re-enables periodic rendering unless the snapshot captured it on. When the viewer returns to live, `ensureWsInit()` re-uploads the current cell vectors plus the enabled/disabled bit so the backend stays aligned.
- `stopSimulation()` sends STOP, updates UI state, clears bond latches, and reinstates idle compute listeners.

## Interaction gating & latching
- Dragging atoms adds their indices to `draggingAtoms`; on release the viewer sends a full `positions` update and records the indices in `modifiedByVersion`.
- Bond rotations build a `bondLatch` set (side atoms) to continue sending sparse updates until rotation completes; latching ensures server frames don’t overwrite those atoms immediately afterward.
- `latchedUntil` and `rotatingAtoms` map indices to expiry timestamps; `buildExcludeSet()` purges expired items before returning the active exclusion set.
- `userInteractionVersion` increments on every user edit. Incoming frames with `userInteractionCount` lower than `lastAppliedUIC` are ignored to avoid reverting freshly edited atoms.

## UI & state synchronization
- Desktop panel controls (temperature slider, MD/Relax buttons, toggles) call into exported orchestration helpers that adjust `cfg` values, call `startContinuous`, `mdStep`, `relaxStep`, or `stopSimulation`.
- Temperature display (`instTemp`) and forces toggle respond to bus events and WS frames.
- Optional test hooks (`window.__WS_TEST_HOOK__`, `setTestHook`) allow integration and E2E tests to inspect outgoing/incoming WS messages without modifying runtime logic.
- Reconnect banner uses `showErrorBanner` / `hideErrorBanner` to present loss-of-connection status; `reconnectNow` triggers an immediate `ensureConnected`.

## Timeline DVR overlay
- **Hidden rail UI**: `public/ui/timelineOverlay.js` injects a hover-revealed rail anchored to the canvas bottom edge. It contains an index label, discrete slider (tick labels in offsets from “live”), and `Play/Pause/Live` buttons. Pointer enter/leave animates between a 6 px hint and a 74 px expanded panel.
- **Frame store**: `public/core/timelineStore.js` keeps a ring buffer (default capacity 500) of recent frames. Each frame stores flattened `Float32Array` payloads (positions/velocities/forces), metadata (kind, sim step, counters, timestamps), and a `cell` record that includes the `enabled` flag so periodic visibility survives playback.
- **Controller**: `public/core/timelineController.js` exposes modes (`live`, `scrub`, `paused`, `playback`), tracks the playback index, and drives a 20 fps timeline via `requestAnimationFrame`. It emits events that the overlay and orchestrator consume to keep UI state, button disabled states, and slider bounds in sync.
- **Playback semantics**:
  - Scrubbing a frame transitions to `paused`, applies the stored positions/cell snapshot locally, bumps the local interaction counter, issues a `STOP_SIMULATION`, and sets `ignoreStreamFrames` so live frames are buffered but not rendered.
  - `Play` resumes from the selected index, iterating frames until the live head. While playing the Play button disables, Pause and Live remain available, and incoming WS frames continue to queue so ack/backpressure stays healthy.
  - Upon reaching the live frame the controller emits `playback-complete`; the orchestrator re-runs `ensureWsInit()` and restarts the remembered MD/relax mode so streaming continues seamlessly.
  - `Live` jumps directly to the newest buffered frame, reapplies it, clears playback state, and restarts streaming immediately.
  - `Pause` simply keeps the viewer detached while still buffering new frames.
- **Testing**: The Playwright suite `ws-timeline-dvr.spec.js` covers scrubbing, bond retention, playback cadence, slider alignment, rapid scrub performance, and cell visibility (both default “periodic off” and lattice-on scenarios).

## Supporting modules (public/)
- `domain` directory contains other services (eventBus, selection, manipulation).
- `render` holds Babylon-specific scene/camera/molecule view logic.
- `core/pickingService.js` centralizes selection and drag pipelines.
- `ui/*` modules build desktop panel widgets, toggles, sliders, and touch controls.
- `util/` includes molecule loading, periodic boundary helpers (`pbc.js`), constants, and element metadata. These utilities are consumed by the orchestrator, domain services, and UI.
- `proto/` hosts the compiled protobuf schemas consumed by `fairchem_ws_client.js`.

## Error handling & diagnostics
- Viewer logs to the console when `window.__MLIPVIEW_DEBUG_API` or `?debug=1` is set (positions, velocities, counter updates).
- WS notices (`WAITING_FOR_ACK`, `SIMULATION_STOPPED`, `COMPUTE_ERROR`) propagate through `handleStreamFrame` so the UI can pause loops and show banners when needed.
- Safety checks (`SAFE_SPHERE_RAD`, displacement sanity, per-atom gating) reduce the risk of runaway positions when the backend sends surprising frames.
