# Timeline Replay Guide

The timeline system lets users pause a live UMA Fairchem stream, inspect any of the last 500 frames, and resume without losing state. It is implemented entirely on the client so scrubbing never blocks on new backend requests.

## Goals & User Experience
- Surface a hover-revealed dock at the bottom of the viewer with play, pause, live, and range controls.
- Keep the most recent 500 authoritative frames (positions, forces, energy metadata) in memory so playback is byte-identical to the live stream.
- Enter a read-only mode when a historical frame is selected: geometry edits are blocked, but camera rotation/zoom gestures stay enabled for inspection.
- Maintain the existing stream state so clicking **Live** or reaching offset −1 resumes the prior MD/relax loop without reloading the molecule.

## Implementation Overview
- **Frame buffer** (`public/core/frameBuffer.js`)
  - Circular buffer keyed by negative offsets (`-1` latest → `-capacity` oldest).
  - Stores positions/velocities/forces as `Float32Array`s plus metadata (energy, stress, counters, derived energy index).
  - `listOffsets()` powers the slider bounds; `getByOffset()` materialises a frame for replay.
- **Timeline UI** (`public/ui/timeline.js`)
  - Creates the dock, a hover hitbox, and the control row (`timeline-play`, `timeline-pause`, `timeline-live`).
  - Slider uses pointer events so single-clicks map directly to the nearest stored offset; scrubbing with pointer move updates continuously.
  - Exposes `refresh()`, `setMode()`, `setActiveOffset()`, and `getState()` for orchestration and Playwright hooks.
- **Viewer orchestration** (`public/index.js`)
  - `handleStreamFrame` records live frames into the buffer before applying them.
  - `enterTimelineMode()` transitions into timeline mode, halts active simulations, enables the read-only overlay, and applies the requested frame.
  - `resumeLiveFromTimeline()` restores the previous continuous mode (idle/md/relax), clears overlays, and resumes incoming frames.
  - `applyTimelineFrame()` sets the active offset, applies stored positions, and updates the energy marker without emitting user interactions.

## Mode & Control Flow
State is tracked in `timelineState` (active, playing, offset, interval timer):

| Mode | Description | Controls Enabled |
| --- | --- | --- |
| **Live** | Default streaming state. | Pause only (`Play`/`Live` disabled) |
| **Paused** | Timeline active, playback halted at a specific offset. | Play + Live |
| **Playing** | Timeline auto-advances toward offset −1 at 20 fps. | Pause + Live |

Key transitions:
- Pointer down on the slider or pause button invokes `handleTimelineOffsetRequest()` → `enterTimelineMode()` → `applyTimelineFrame()`.
- Playback (`Play`) sets `timelineState.playing=true` and schedules `timelineStepForward()` to walk offsets until it hits `-1`, then resumes live mode automatically.
- `Live` triggers `resumeLiveFromTimeline()` which re-enables interactions and restarts the last continuous loop if one was running.

## Interaction Policy
- **Geometry edits:** Disabled while timeline mode is active. Manipulation services stay latched off so atom drags and bond rotations are rejected.
- **Camera controls:** Remain enabled so viewers can rotate/zoom historical frames. The picking layer explicitly keeps navigation enabled during timeline playback.
- **Energy plot:** Historical frames set a yellow marker at the recorded energy index; returning to live mode clears the marker.
- **Overlay:** A translucent overlay communicates the read-only state and suppresses accidental mesh interaction.

## Pointer & Slider Behaviour
- The slider maps clientX to offsets using the current offset list; scrubbing always clamps to the nearest recorded frame.
- Pointer-up events call `handleOffsetRequest` once more and suppress the subsequent `change` event so single clicks do not require double interaction.
- The Playwright regression `tests-e2e/ws-timeline-slider-select.spec.js` covers this mapping to protect against regressions where the first click landed on offset −1.

## Diagnostics & Test Hooks
- `window.viewerApi.timeline` exposes:
  - `select(offset)`, `play(offset)`, `pause()`, `live()` delegating to timeline handlers.
  - `getState()` returning `{mode, offset, active, playing}`.
  - `getSignature(offset)` for byte-level equality checks during playback.
  - `getOffsets()` for assertions on the backing buffer.
- Additional coverage:
  - `ws-timeline-controls.spec.js`, `ws-timeline-interaction-lock.spec.js`, `ws-timeline-visibility.spec.js`, `ws-timeline-replay.spec.js`, `ws-timeline-camera.spec.js`, `ws-timeline-energy-marker.spec.js`.

## Operational Notes
- When entering timeline mode, the viewer acknowledges the latest live frame **before** stopping the simulation to avoid triggering backend backpressure.
- Frame buffer capacity is configurable (default 500) and can be tuned via `installTimeline({ capacity })`.
- The system tolerates sparse frame lists; when the buffer is empty the slider disables itself and the dock stays hidden.
