# State System Update (Release V1)

## Findings From Current Review
- Timeline playback in `public/index.js` is hard-coded to a 50 ms interval (20 fps) with no hook for dynamic speed, looping, or external orchestration; control of the dock lives entirely inside `index.js`.
- `SessionStateManager` exports `schemaVersion: 1`, zeroes counters on restore, and calls `seedSequencing` with hard-coded zeros; the saved JSON fixtures (`fixtures/sessionSnapshots/*.json`) already use schema v2, so the mismatch drops metadata during load.
- Frontend state is still split across ad-hoc singletons (`stateStore`, timeline globals, viewer closures). Control surfaces (`viewerApi`) have no way to expose additional playback state or UI annotations.
- Backend session state (`fairchem_local_server2/session_models.py`, `ws_app.py`) remains isolated from new timeline metadata, but the JSON loader always pushes a dense `full_update`; resuming loops after loading still depends on the client remembering `lastContinuousOpts`.
- There is no library selector for curated sessions in `public/examples/`; sample manifest `public/examples/library.json` is unused.

## JSON Session Schema (v3 Proposal)

### Top-level shape

```json
{
  "schemaVersion": 3,
  "savedAt": "2025-02-15T18:42:33.410Z",
  "source": { "kind": "json", "label": "SN2 demo" },
  "viewer": { "...": "unchanged geometry/cell payload" },
  "energyPlot": { "series": [...], "markerIndex": 14 },
  "timeline": {
    "capacity": 500,
    "frames": [
      {
        "id": "frame-0001",
        "seq": 401,
        "kind": "md",
        "simStep": 1,
        "userInteractionCount": 92,
        "timestamp": 1761703692101,
        "energy": -31987.19,
        "positions": [[...]]
      }
      /* oldest → newest */
    ],
    "playback": {
      "startFrame": { "frameId": "frame-0120", "offset": -42 },
      "autoPlay": true,
      "loop": true,
      "loopRange": { "startFrameId": "frame-0080", "endFrameId": "frame-0145" },
      "defaultFps": 20
    },
    "controlMessages": [
      {
        "id": "focus-approach",
        "priority": 0,
        "range": {
          "start": { "frameId": "frame-0075" },
          "end": { "frameId": "frame-0105", "inclusive": true }
        },
        "actions": [
          { "type": "timeline.playbackSpeed", "fps": 10, "transitionMs": 150 },
          {
            "type": "overlay.callout",
            "text": "Approach phase\nNucleophile forming bond",
            "textSize": 0.5,
            "panelSize": { "width": 1.2, "height": 0.6 },
            "anchor": { "mode": "bond", "atoms": [2, 15], "offset": [0, 0.8, 0] },
            "style": { "background": "rgba(12,18,32,0.85)", "color": "#f0f7ff" }
          },
          {
            "type": "visual.opacityFocus",
            "focus": { "atoms": [2, 15, 18], "includeBonds": "connected" },
            "focusOpacity": { "atoms": 1.0, "bonds": 1.0 },
            "backgroundOpacity": { "atoms": 0.12, "bonds": 0.05 },
            "transitionMs": 200
          }
        ]
      }
    ]
  },
  "websocket": {
    "nextSeq": 413,
    "lastAck": 409,
    "userInteractionCount": 37,
    "totalInteractionCount": 51,
    "simStep": 98
  }
}
```

### Boundary semantics
- `timeline.frames` continue to serialize oldest → newest. Each frame gains a stable string `id` (`frame-0001`, `frame-0002`, …) so ranges stay readable and resilient; the numeric `id` is preserved for backward compatibility.
- `timeline.playback.startFrame` accepts any combination of `frameId`, `frameIndex` (0 = oldest), or `offset` (-1 = newest). The loader resolves in that order.
- `loopRange` defaults to the full frame list when omitted. If `loop` is false, playback stops at the newest frame and resumes live mode as today.
- `controlMessages[].range.start|end` accept `frameId`, `frameIndex`, or `offset`. `inclusive` defaults to `true` for `end`, `false` for `start`.
- `priority` controls action overrides when multiple messages overlap (higher wins per action type). Absent priority defaults to `0`.

### Action catalog

| Action type                 | Fields (required → optional)                                                                                                             | Effect |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| `timeline.playbackSpeed`    | `fps` **or** `speedMultiplier`; optional `transitionMs`, `easing`                                                                         | Overrides the scheduler cadence while active; transitions are eased when provided. |
| `overlay.callout`          | `text`, `anchor` (`mode: 'world' | 'atom' | 'bond'`, `position` or `atoms`, optional `offset`), optional `textSize`, `panelSize`, `style` | Creates or updates a billboarded mesh that always faces the active camera and renders above atoms/bonds. |
| `visual.opacityFocus`      | `focus` (`atoms`, optional `bonds`, `includeBonds`), `focusOpacity` (`atoms`, `bonds`), `backgroundOpacity`, optional `transitionMs`      | Applies per-instance alpha buffers to highlight selected atoms/bonds while fading everything else. |

Future actions (camera bookmarks, HUD toggles) can append to `actions` without altering the schema.

## Frontend Design

### Session state refactor
- Bump `SessionStateManager.SCHEMA_VERSION` to 3 and ensure save/load round-trips include `timeline.playback` and `timeline.controlMessages`.
- Replace the hard-coded counter reset with the persisted values (`websocket.userInteractionCount`, `totalInteractionCount`, `simStep`). Expose these via `viewerApi.session.getLastSnapshot()` for diagnostics.
- Store `controlMessages` and `playback` in the baseline snapshot so `resetToLastSource()` restores annotations and playback policy alongside geometry.
- Normalize inputs when loading: validate `fps`, filter control messages with invalid ranges, and emit console diagnostics tagged `[SessionManager][controlMessages]`.

### Timeline orchestration
- Introduce `createTimelinePlaybackController` (new module in `public/core/timelinePlayback.js`). Responsibilities:
  - Resolve playback start frame, initiate `enterTimelineMode`, and schedule the first frame render.
  - Run a dynamic scheduler using `requestAnimationFrame` plus drift-aware timers so per-frame spacing adheres to the current `fps`.
  - Apply `loop` and `loopRange`, jumping back to the loop start instead of leaving timeline mode when enabled.
  - Surface hooks `setBaseFps`, `setSpeedOverride`, `setLoopConfig`, `setAutoPlay`, and `stop`.
- Replace the existing `timelineState.interval` logic with the new controller. Retain legacy behaviour (resume live mode) when neither loop nor auto-play is requested.

### Control message engine
- Add `public/core/controlMessageEngine.js` to ingest the session’s `controlMessages`, map numeric/string frame references to buffer indices, and compute per-frame active actions.
- Expose methods:
  - `prime({ frames, geometryAccessor })` – preprocess frame ranges and precompute any static anchors.
  - `evaluate(frameId)` – return the set of actions active for the requested frame.
  - `onFrameApplied(frameId, framePayload)` – invoke action handlers (speed, callouts, opacity) with the latest geometry.
  - `onTimelineCleared()` – tear down callouts and restore baseline opacities.
- Install the engine inside `public/index.js` alongside the timeline buffer. Whenever `applyTimelineFrame` succeeds, forward the current frame metadata to `controlEngine.onFrameApplied`.
- Speed overrides feed into the playback controller; the controller selects the highest-priority active `timeline.playbackSpeed` action each step.

### Callout overlay layer
- Create `public/render/calloutLayer.js` that builds/maintains billboard planes:
  - Uses Babylon `MeshBuilder.CreatePlane` with `billboardMode = BABYLON.Mesh.BILLBOARDMODE_ALL`.
  - Text rendered via `AdvancedDynamicTexture` (or HTML fallback in headless tests) so multi-line strings keep newline formatting.
  - Assign `renderingGroupId` > bonds/atoms so callouts always draw last.
  - Updates anchor coordinates each frame by reading current atom/bond positions (`moleculeView.getAtomPosition(i)` / midpoint helper).
- Provide `controlEngine` with a thin wrapper for create/update/destroy to keep UI logic separated from geometry state.

### Opacity transitions
- Extend `moleculeView` with `setOpacityMask({ focusAtoms, focusOpacity, backgroundOpacity })` so we can apply alpha changes without mutating the underlying bond/atom collections.
- Cache previous opacities; when no control message applies, restore the original transparency (respecting the existing "atoms render last" rule).
- Ensure transitions interpolate alpha in the shader buffer rather than re-creating meshes; this keeps render order untouched.

### System → Library dropdown
- Read `public/examples/library.json` at startup (lazy via `fetch` when the System panel first expands). Store results in `viewerState.library`.
- Add a `select` element labelled `Library` beneath the Load/Save buttons. Each option shows the manifest `label`; hover tooltip displays `description`.
- Selecting an entry loads the JSON via `fetch(entry.path)` and passes the blob to `session.loadFromFile`. Handle errors with `showErrorBanner`.
- For smoke-test compatibility, expose `viewerApi.session.loadFromLibrary(id)` so Playwright can bypass the DOM.

## Backend Coordination
- No protobuf changes required. Keep the backend oblivious to UI control metadata.
- During snapshot load, continue to push a dense `full_update` but reuse the stored counters instead of zeroing them; this avoids unnecessary `WAITING_FOR_ACK` states.
- Add a debug log in `ws_app.py` gated by `MLIPVIEW_RESUME_DEBUG` when a restored client re-sends a `full_update`, so we can correlate with playback-led resumes.

## Regression Testing

### Jest / unit
- `tests/controlMessageEngine.spec.js`: feed synthetic frame sets and control messages, assert active actions for interior/start/end frames, overlap resolution via priority, and graceful handling of unknown frame IDs.
- `tests/timelinePlaybackController.spec.js`: simulate time progression, verify fps overrides apply, loop restarts at the configured start, and auto-play toggles call the frame callback.
- Extend `tests/sessionStateManager.spec.js` to assert schema v3 round-trips playback metadata, counters persist, and invalid control messages are dropped with warnings.
- `tests/moleculeView.opacityMask.spec.js`: verify focus/background opacities update buffers without breaking render order.

### Playwright
- `ws-library-load.spec.js`: select the SN2 library entry, wait for the callout to appear, confirm auto-play kicks in, then ensure the loop repeats without returning to live mode.
- `ws-control-messages.spec.js`: scrub into the middle of a control range, assert `viewerApi.debug.getTimelineState()` reports the overridden fps, verify opacity changes via `viewerApi.debug.getRenderState()`, and ensure callout text matches multiline expectations.
- `ws-session-start-frame.spec.js`: load a session with `autoPlay` disabled, confirm the specified `startFrame` is displayed on load, trigger `viewerApi.timeline.play`, and ensure playback resumes from that offset.

### Backend (pytest)
- Add a smoke test that loads a v3 snapshot via the existing JSON loader helper, confirms the returned `SessionState` counters match the snapshot, and that the next MD frame is accepted without triggering backpressure.

## Options & Trade-offs
1. **Dedicated control-message engine (recommended)**  
   - *Pros*: Encapsulates preprocessing, prioritization, and lifecycle; keeps `index.js` manageable; straightforward to unit test; future actions (camera cues) slot in cleanly.  
   - *Cons*: Adds another moving part; requires new dependency wiring in `index.js`.
2. **Inline control handling inside `index.js`**  
   - *Pros*: Fewer files touched up front.  
   - *Cons*: Exacerbates the existing mega-module, harder to compose with tests, and risks regressions when more control types arrive.

For callout rendering we evaluated HTML overlays vs. Babylon GUI meshes. Mesh-based overlays keep depth ordering consistent with VR/AR paths; HTML overlays would require duplicate logic for XR, so we stick with the Babylon approach.

## Implementation Checklist
1. Adjust `SessionStateManager` (schema v3, counter restore, playback/control serialization) and update fixtures.
2. Introduce `timelinePlaybackController` and replace the interval-based loop in `public/index.js`.
3. Add `controlMessageEngine` plus the callout/opacity helpers; wire them into timeline playback and viewer initialization.
4. Teach the System panel about the Library manifest and add `viewerApi.session.loadFromLibrary`.
5. Update `public/examples/sn2/sn2.json` (and other fixtures) to schema v3 with sample control messages.
6. Extend Jest + Playwright + pytest coverage as outlined.
7. Document the new schema and hooks in `testing.md` / `test_hooks.md` where relevant.

Upon approval of this plan, we can proceed with implementation (Step 6) followed by converting `sn2.json` to the new format (Step 7).
