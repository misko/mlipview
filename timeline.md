# Timeline Replay Guide

The timeline system lets users pause a live UMA Fairchem stream, inspect any of the last 500 frames, and resume without losing state. Playback is fully client-side so scrubbing never blocks the viewer on new backend requests.

## Goals & User Experience
- Surface a hover-revealed dock at the bottom of the viewer with play, pause, live, and loopable controls.
- Keep the most recent 500 authoritative frames (positions, forces, energy metadata) in memory so playback is byte-identical to the live stream.
- Enter a read-only mode when a historical frame is selected: geometry edits are blocked while camera rotation/zoom gestures stay enabled.
- Respect scripted playback: optional JSON control messages can auto-play, loop, display callouts, adjust playback speed, and highlight subsets of atoms/bonds during a walkthrough.
- Maintain the existing stream state so clicking **Live** or reaching offset −1 resumes the prior MD/relax loop without reloading the molecule.

## Implementation Overview
- **Frame buffer** (`public/core/frameBuffer.js`)
  - Circular buffer keyed by negative offsets (`-1` latest → `-capacity` oldest).
  - Frames carry stable string IDs (`frame-00001`, …) alongside numeric metadata so session snapshots and control ranges can address frames deterministically.
  - Provides helpers to map IDs ↔ offsets/indices (`resolveOffset`, `resolveFrameIndex`) for the control-message engine.
- **Timeline UI** (`public/ui/timeline.js`)
  - Builds the dock, hover hitbox, and the control row (`timeline-play`, `timeline-pause`, `timeline-live`).
  - Slider uses pointer events so single-clicks map directly to the nearest stored offset; scrubbing with pointer move updates continuously.
  - Exposes `refresh()`, `setMode()`, `setActiveOffset()`, and `getState()` for orchestration and Playwright hooks.
- **Playback controller** (`public/core/timelinePlaybackController.js`)
  - Drift-aware scheduler that advances frames according to the current FPS.
  - Supports auto-play on load, looping over an optional sub-range, and temporary speed overrides injected by control messages.
- **Control message engine** (`public/core/controlMessageEngine.js`)
  - Preprocesses the session’s `timeline.controlMessages[]` and determines, per frame, which actions should fire (speed override, callout, opacity mask).
  - Supplies the playback controller, callout layer, and molecule opacity mask with the active directives.
- **Callout layer** (`public/render/calloutLayer.js`)
  - Renders floating annotations as Babylon billboards that always face the viewer and anchor to world/atom/bond positions.
- **Viewer orchestration** (`public/index.js`)
  - `handleStreamFrame` records live frames into the buffer before applying them.
  - `enterTimelineMode()` transitions into timeline mode, halts active simulations, enables the read-only overlay, and applies the requested frame.
- `applyTimelineFrame()` sets the active offset, applies stored positions, triggers control message side effects (callouts/opacity/speed), and updates the energy marker without emitting user interactions.
- `resumeLiveFromTimeline()` restores the previous continuous mode (idle/md/relax), clears overlays, and resumes incoming frames with defensive logging (`[Timeline][resumeLiveFromTimeline] …`).

## Authoring Mode (`?edit=1`)
- Append `?edit=1` to the viewer URL to unlock the timeline authoring toolkit. A right-rail editor lists all control messages, exposes their playback ranges, and renders dedicated editors for playback speed, callouts, and opacity masks.
- The editor persists changes immediately on save: messages are written back to the control-message engine, the active frame is reapplied so callouts/opacity updates are visible straight away, and the session baseline is refreshed so JSON exports pick up new metadata.
- Playback presets (default FPS, auto-play, loop range, start frame) are editable at the top of the panel. The editor calls into `timelinePlayback` so loop bounds and cadence update in real time without reloading.
- A bottom status bar mirrors the currently selected frame (`Live` vs `frame-0012 · offset -4`) and highlights the active atom/bond selection. Utility buttons in the editor allow authors to pull the current frame offset or selection into range/action fields.
- Viewer API surface: `viewerApi.timelineEditor` exposes `refresh`, `getState`, `select(id)`, `getDraft()`, and `status()` helpers for tests and automation. `viewerApi.timeline.getFrameMeta(offset)` returns `{ frameId, offset, frameIndex }` to make assertions easier.

## Mode & Control Flow
State is tracked in `timelineState` (active, playing, offset) plus the playback controller’s runtime (effective FPS, loop bounds, start frame).

| Mode | Description | Controls Enabled |
| --- | --- | --- |
| **Live** | Default streaming state. | Pause only (`Play`/`Live` disabled) |
| **Paused** | Timeline active, playback halted at a specific offset. | Play + Live |
| **Playing** | Timeline auto-advances according to the effective FPS (default 20 fps, control messages may override). | Pause + Live |

Key transitions:
- Pointer down on the slider or pause button invokes `handleTimelineOffsetRequest()` → `enterTimelineMode()` → `applyTimelineFrame()`.
- Playback (`Play`) calls `timelinePlayback.start()`, which schedules `timelineStepForward()`; looping ranges restart automatically, otherwise playback returns to live mode on offset −1.
- `Live` triggers `resumeLiveFromTimeline()` which re-enables interactions and restores the last continuous run (idle/MD/relax) if one was active.

## Interaction Policy
- **Geometry edits:** Disabled while timeline mode is active. Manipulation services stay latched off so atom drags and bond rotations are rejected.
- **Camera controls:** Remain enabled so viewers can rotate/zoom historical frames. The picking layer explicitly keeps navigation responsive during timeline playback.
- **Energy plot:** Historical frames set a yellow marker at the recorded energy index; returning to live mode clears the marker.
- **Overlay:** A translucent overlay communicates the read-only state and suppresses accidental mesh interaction.
- **Control messages:** Optional JSON directives can override playback speed, display callouts, and fade non-target atoms/bonds while their range is active. Outside the range the viewer restores the baseline opacity/callout state automatically.

## Pointer & Slider Behaviour
- The slider maps clientX to offsets using the current offset list; scrubbing always clamps to the nearest recorded frame.
- Pointer-up events call `handleOffsetRequest` once more and suppress the subsequent `change` event so single clicks do not require double interaction.
- The Playwright regression `tests-e2e/ws-timeline-slider-select.spec.js` covers this mapping to guard against the historical “double-click to select” regression.

## Diagnostics & Test Hooks
- `window.viewerApi.timeline` exposes:
  - `select(offset)`, `play(offset)`, `pause()`, `live()` delegating to timeline handlers.
  - `getState()` returning `{ mode, offset, active, playing }`.
  - `getSignature(offset)` for byte-level equality checks during playback.
  - `getOffsets()` for assertions on the backing buffer.
  - `getPlaybackConfig()` and `getControlState()` so tests can assert control-message effects.
- Additional coverage:
  - `ws-timeline-controls.spec.js`, `ws-timeline-interaction-lock.spec.js`, `ws-timeline-visibility.spec.js`, `ws-timeline-replay.spec.js`, `ws-timeline-camera.spec.js`, `ws-timeline-energy-marker.spec.js`.
  - `ws-session-playback-resume.spec.js` saves a JSON snapshot, reloads it, scrubs five frames back, plays forward, and verifies the viewer receives a bounded burst of fresh MD frames before handing control back to the live stream.
  - Unit suites `tests/controlMessageEngine.spec.js` and `tests/timelinePlaybackController.spec.js` guard the control/message plumbing.

## Operational Notes
- When entering timeline mode, the viewer acknowledges the latest live frame **before** stopping the simulation to avoid triggering backend backpressure.
- Frame buffer capacity is configurable (default 500) and can be tuned via `installTimeline({ capacity })`.
- Play/pause/loop state can be pre-seeded by JSON snapshots (`timeline.playback`). When `autoPlay` is true the viewer enters timeline mode and begins playback immediately after loading.
- Control messages are optional; when absent the playback controller falls back to fixed 20 fps and no callouts/opacity masks are applied.
- Tests cap resume bursts to 170 frames (via `viewerApi.startMDContinuous` interception) so Playwright coverage remains fast while still proving the stream resumes correctly.
