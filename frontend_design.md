# Frontend Architecture (current)

## Entry point & initialization
- `public/index.js` exports `initNewViewer(canvas, initialData)`, which bootstraps the Babylon.js scene, Molecule state, services, the timeline dock, and WebSocket wiring.
- On load it normalizes runtime configuration (drag throttles, min step interval, temperature defaults) via `window.__MLIP_CONFIG`.
- `createMoleculeState()` (domain/moleculeState.js) builds the authoritative state object holding elements, position vectors, cell metadata, event bus, and version counters.
- Initial geometry is clamped into a "safe sphere" if configured, baseline energy/forces are fetched, timeline history is seeded, and optional initial cell snapshots are stored for reset.

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
- `handleStreamFrame(kind, ws, frame)` records the frame into the timeline buffer, then applies positions (respecting per-atom gating), updates forces/energy/temperature, pushes energy to the plot, acks the sequence, and updates `lastAppliedUIC`.
- `stopSimulation()` sends STOP, updates UI state, clears bond latches, and reinstates idle compute listeners.

### Timeline controller
- `createFrameBuffer()` powers the timeline history (bonded to a 500 frame capacity by default).
- `installTimeline()` builds the UI dock (`timeline.js`), exposes `viewerApi.timeline`, and wires button/slider callbacks to `handleTimelineOffsetRequest`, `handleTimelinePlayRequest`, and `handleTimelineLiveRequest`.
- Timeline mode introduces a modal state (`Mode.Timeline`) that disables molecule editing while keeping camera controls active. The overlay and manipulation facade enforce the read-only policy until `resumeLiveFromTimeline()` runs.
- Pointer-driven scrubbing maps screen coordinates to stored offsets, ensuring one click selects the intended frame; tests cover the mapping via `ws-timeline-slider-select.spec.js`.

## Interaction gating & latching
- Dragging atoms adds their indices to `draggingAtoms`; on release the viewer sends a full `positions` update and records the indices in `modifiedByVersion`.
- Bond rotations build a `bondLatch` set (side atoms) to continue sending sparse updates until rotation completes; latching ensures server frames don’t overwrite those atoms immediately afterward.
- `latchedUntil` and `rotatingAtoms` map indices to expiry timestamps; `buildExcludeSet()` purges expired items before returning the active exclusion set.
- `userInteractionVersion` increments on every user edit. Incoming frames with `userInteractionCount` lower than `lastAppliedUIC` are ignored to avoid reverting freshly edited atoms.

## UI & state synchronization
- Desktop panel controls (temperature slider, MD/Relax buttons, toggles) call into exported orchestration helpers that adjust `cfg` values, call `startContinuous`, `mdStep`, `relaxStep`, or `stopSimulation`.
- Temperature display (`instTemp`) and forces toggle respond to bus events and WS frames.
- Timeline controls expose `viewerApi.timeline` helpers for Playwright and debugging, including `getState()` and `getSignature(offset)`.
- Optional test hooks (`window.__WS_TEST_HOOK__`, `setTestHook`) allow integration and E2E tests to inspect outgoing/incoming WS messages without modifying runtime logic.
- Reconnect banner uses `showErrorBanner` / `hideErrorBanner` to present loss-of-connection status; `reconnectNow` triggers an immediate `ensureConnected`.

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
