# Timeline Replay Design

## Goals
- Surface a hover-revealed timeline along the bottom edge of the viewer that can scrub through the most recent 500 authoritative frames.
- Persist the last 500 protobuf `ServerResult` frames (with positions) entirely on the client so a user can inspect prior states without a server round-trip.
- Provide `Play`, `Pause`, and `Live` controls. `Play` replays stored frames at 20 fps, `Pause` freezes on the selected frame, and `Live` resumes real-time streaming from the backend.
- When a historical frame is selected, immediately stop consuming and emitting WebSocket traffic, inject the stored frame as the active viewer state, and block further geometry edits until the user returns to live mode.
- Guarantee frame identity: selecting the same offset twice should yield byte-identical positions (needed for the Playwright assertion).
- Leave the system in a safe live-streaming state once `Live` is tapped, with new frames continuing to fill the 500-frame buffer.

## Current Behaviour
- `public/index.js` owns the streaming state machine (`Mode.Idle|MD|Relax`), the per-frame handler (`handleStreamFrame`), and the idle/continuous subscription wiring (`startContinuous`, `attachIdleWSListener`).
- Frames are applied in-place; no history is retained once `applyTriples` mutates the molecule state.
- WebSocket sequencing expects timely ACKs; `handleStreamFrame` immediately calls `ws.ack(seq)` after accepting a frame.
- User interactions (dragging, bond rotation) are still enabled whenever the viewer is idle, so simply rewinding state today would race against in-flight edits.
- There is no UI real estate reserved along the bottom of the canvas; CSS is injected on demand via `desktopPanel.js`.

## Proposed Changes

### Frame Buffer Service
- Add `public/core/frameBuffer.js` exporting `createFrameBuffer({ capacity = 500 })`.
  - Stores frames in a circular buffer of `capacity` entries keyed by monotonic `frameId`.
  - Each stored frame clones:
    - `kind` (`'md' | 'relax' | 'idle'`), `timestamp`, `simStep`, `userInteractionCount`, `seq`.
    - `positions` as an array of `Float32Array` triples (ensures no sharing with live state).
    - Optional `velocities`, `forces`, `stress`, `energy`.
    - A `signature` (e.g., `Uint32Array` created by hashing the first N triples) used by tests to prove equality.
  - Provides helpers:
    - `record(kind, payload)` returns `{ id, offset: -1 }`.
    - `getByOffset(offset)` with offsets `-1` (latest) to `-capacity`.
    - `latestOffset()` and `size()`.
    - `listOffsets()` for UI tick generation.
  - Ignores frames that do not include `positions`.

### Timeline Controller
- Add `public/ui/timeline.js` exporting `installTimeline({ host, frameBuffer, applyFrame, requestLiveResume, blockInteractions, unblockInteractions })`.
  - Renders a docked `<div id="timelineDock">` anchored to `bottom: 0`.
  - Hidden by default via CSS transform; `:hover` on a 20 px-high hit strip reveals it (also expose explicit `data-visible` for programmatic control).
  - Child structure:
    - Control bar with `button[data-action="play"]`, `button[data-action="pause"]`, `button[data-action="live"]`.
    - A discrete slider implemented as `input[type=range]` with `step=1`, `min=-capacity`, `max=-1`. CSS overlays tick marks every 10 frames and labels `-1`, `-250`, `-500`.
    - `div.timeline-ticks` showing small dividers for accessibility. Each tick element carries `data-offset`.
  - Event flow:
    - On hover reveal, slider syncs to the latest offset.
    - On slider change or tick click, fires `onFrameRequested(offset)`.
    - Buttons dispatch `play`, `pause`, `live` events.
  - Expose a testing bridge via `window.viewerApi.timeline = { getOffsets(), getActiveOffset(), getFrameSignature(offset), isLive(), play(), pause(), goLive(), setOffset(offset) }`.
  - Add `data-testid` attributes (`timeline-dock`, `timeline-play`, `timeline-slider`, etc.) to keep Playwright selectors stable.

### Viewer State Machine Extensions
- Extend `Mode` in `public/index.js` with `Timeline`.
  - Add `setMode(Mode.Timeline)` branch to:
    - Detach idle/MD/relax WS listeners.
    - Send `ws.stopSimulation()` if streaming was active.
    - Prevent auto-resume logic from re-arming (skip `rememberResume`).
  - When switching into `Timeline`, call `blockInteractions()` to:
    - Attach a full-screen overlay element (`#timelineOverlay`) intercepting pointer events.
    - Expose a boolean to picking/manipulation services so programmatic drags (e.g. tests) bail early.
  - `Live` button transitions back to `Mode.Idle` (or restarts the last streaming mode if it was active when the user scrubbed) and replays any pending resume info.

### Frame Application Unification
- Refactor the body of `handleStreamFrame` that mutates state into `applyServerFrame({ kind, frame, acknowledge = true })`.
  - `handleStreamFrame` becomes a thin wrapper that stores the frame (`frameBuffer.record(kind, frame)`) before calling the shared applicator.
  - `applyServerFrame` accepts an option to skip ACK for replayed frames.
  - Replay path (timeline) calls `applyServerFrame({ kind: stored.kind, frame: stored, acknowledge: false })`. This ensures energy/force caches, temperature labels, and energy plots update consistently.
  - While the viewer is in timeline mode, suppress `energyPlot.push` to avoid polluting the live trend (tracked by a flag).

### Playback Loop
- Timeline controller owns a `setInterval` ticking every 50 ms (20 fps) while playing.
  - Starting from the selected offset, it increments toward `-1`. When reaching `-1`, it auto-pauses (remains in timeline mode until the user presses `Live`).
  - Each tick obtains the next offset, invokes `applyFrame`, updates the active slider value, and records the new active offset.
  - `Pause` clears the interval, maintaining the current offset selection.
  - Playback guards against frames disappearing (e.g., buffer underflow) by clamping to available offsets.

### WebSocket Coordination
- When a historical frame is requested:
  1. Immediately set mode to `Timeline`.
  2. Call `ws.stopSimulation()` inside a try/catch. Do not send further `userInteraction` payloads while in this mode.
  3. Stop scheduling ACK flushes by detaching the active `onResult` subscription. `FrameBuffer` storage already acknowledged the last live frame before the pause, so the backend sees the stream as cleanly stopped.
  4. Apply the stored frame locally without emitting a new ACK (frames carry stale seq numbers by design).
  - Backend confirmation (`fairchem_local_server2/ws_app.py`): the server only enforces `state.server_seq - state.client_ack < 10` before producing another frame. It happily tolerates a missing ACK for the final live frame because the next ACK we send after resuming will advance `state.client_ack`. No extra bookkeeping required during rewind.
- `Live` button flow:
  1. Clear timeline playback, unblock interactions, switch mode to `Mode.Idle`.
  2. Re-run `ensureWsInit` to re-upload the current geometry snapshot (now equal to the historical frame).
  3. If a simulation was running before rewinding, restart it with the previous options (`pendingResume` infrastructure already caches `lastContinuousOpts`).
  4. Resume idle listener so new frames refill the ring buffer.

### Interaction Lock
- `blockInteractions()` inserts an absolutely positioned `div` covering the canvas with a subtle backdrop to signal read-only mode; the overlay keeps `pointer-events: none` so camera input, touch gestures, and wheel zoom continue to hit the canvas.
- `unblockInteractions()` removes the overlay when returning to live streaming.
- Picking/manipulation services expose `setInteractionEnabled(bool)` to gracefully reject programmatic drags; `public/index.js` disables only manipulation while timeline mode is active so camera orbit/zoom stays responsive.

### Energy Plot & Metrics
- While in timeline mode, freeze `energyPlot.push`. Instead update HUD labels (`instEnergy`, `instTemp`, `rps`) directly via `applyServerFrame` so the UI reflects the stored frame but the historical playback does not mutate the running energy history.
- On returning to live mode, rebuild the plot using the buffer if necessary (e.g., seed it with the most recent live energy value) and unfreeze pushes.

## Playwright Coverage

### Mandatory Spec
- `tests-e2e/ws-timeline-replay.spec.js`
  1. Launch viewer, ensure MD run for ≥30 frames (use existing helper).
  2. Reveal timeline (`page.hover('[data-testid="timeline-hitbox"]')`).
  3. Click the tick for frame `-15`, wait for `viewerApi.timeline.getActiveOffset()` to equal `-15`.
  4. Record `viewerApi.timeline.getActiveSignature()`.
  5. Click tick `-16`, then `-15` again; confirm the signature matches.
  6. Hit `Play`, let it advance to offset `-1`.
  7. Press `Live`, rerun MD for ≥30 additional frames, asserting frame buffer growth and `viewerApi.getMetrics().running === 'md'` after resume.

### Additional Tests
1. **Buffer Capacity** (`tests-browser` or Playwright hybrid):
   - Stream > 520 frames, ensure `viewerApi.timeline.getBufferStats().size === 500` and the oldest accessible offset is `-500`.
2. **Interaction Lock**:
   - Enter timeline mode, attempt `viewerApi.manipulation.beginDrag`, expect `false` and no outgoing `USER_INTERACTION` frames via test hook.
3. **Live Resume Idle**:
   - From idle mode (no active simulation), scrub backwards, then press `Live` and verify idle frames resume arriving (buffer newest offset updates).
4. **Play From Middle**:
   - Select offset `-100`, press `Play`, ensure it lands on `-1` without throwing and remains in `Mode.Timeline` until `Live` pressed.
5. **UI Visibility & Buttons** (DOM-focused test):
   - Ensure timeline dock is hidden without hover, becomes visible on hover, `Play`/`Pause` buttons toggle `data-state` attributes.
6. **Control Policy Enforcement**:
   - **Live-State Disablement**: While live streaming, assert Play and Live buttons carry `disabled` and Pause is enabled; clicking Pause transitions to timeline paused with overlay visible.
   - **Timeline Paused Buttons**: After selecting an older frame, confirm Play and Live are enabled while Pause is disabled; clicking Live returns to live streaming and re-disables Live/Play.
   - **Timeline Playing Buttons**: Trigger Play from a historical frame, verify Play renders disabled, Pause & Live enabled; click Pause to stop replay, then Play to resume.
   - **Pause During Live (Regression for Pause requirement)**: With MD streaming active, click Pause and assert we enter timeline paused at `-1` with WS stopped.
   - **Auto Resume on Replay Completion**: Start playback from a mid-buffer frame and wait for offset `-1`; verify timeline auto-resumes live (overlay hidden, all buttons in live configuration).

## Risks & Mitigations
- **Memory Pressure**: storing full `forces` arrays multiplies footprint. Use `Float32Array` and reuse typed buffers where possible; clamp capacity to 500.
- **WS Backpressure**: ensure the last live frame is acknowledged before stopping (store/ack first, then detach) so the server does not throttle future runs.
- **Energy Plot Drift**: timeline playback should not distort live metrics. Suppressing `energyPlot.push` while timeline mode is active keeps live trends intact.
- **Concurrency Edge Cases**: guard all transitions with `mode` checks so repeated button presses do not leave the viewer in an inconsistent state. Expose internal state for tests to assert invariants.

## Open Items
- Confirm whether forces are required for timeline playback. Current plan stores them; if UI becomes sluggish we can drop forces from the stored payload later.
- Determine whether we should auto-restart the exact MD/relax loop the user was running upon `Live`. Proposal mirrors existing `pendingResume` behaviour; feedback welcome before implementation.

## Timeline Control Policy

- **Live Streaming**: `timelineState.active=false`, `playing=false`. Buttons: Play disabled, Pause enabled (enters timeline mode at -1), Live disabled.
- **Timeline Paused**: `active=true`, `playing=false`. WebSocket paused, overlay active. Buttons: Play enabled (resumes replay), Pause disabled (no-op), Live enabled (resume streaming).
- **Timeline Playing**: `active=true`, `playing=true`. Buttons: Play disabled (already playing), Pause enabled (stops replay, stay paused), Live enabled (resume streaming).
- **Auto Resume**: When replay advances to offset -1, automatically transition to Live Streaming (re-enable WS, disable overlay).

Buttons should reflect these policies: disable Play in Live/Playing states, disable Pause in Live/Paused states, disable Live only while already live.
