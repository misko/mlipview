# State System Update (Release V1)

## Findings From Review
- Timeline playback state is orchestrated inside `public/index.js#L1700` with no extension points for editing; control-message data lives solely in `createControlMessageEngine` but has no authoring UI or metadata (labels, descriptions). Editing today requires manual JSON edits.
- `SessionStateManager` (`public/core/sessionStateManager.js`) persists timeline playback + control messages, yet it enforces `schemaVersion: 3` and drops unknown fields. Any authoring features must update this module so saves/loads remain authoritative.
- Selection/camera/timeline HUDs are scattered: `state.selection` events are only used to manage bond latches; there is no consolidated place to surface selection metadata or the currently previewed frame, which the editing workflow needs.
- Frontend state is still split across `stateStore`, `timelineState`, and implicit singletons; however, we have sufficient hooks (`controlEngine`, `timelinePlayback`, `frameBuffer`, `sessionStateManager`) to expose an edit surface without a wholesale rewrite.
- Documentation (`timeline.md`, `frontend_design.md`, `backend_design.md`, `testing.md`) already assumes control messages exist but does not describe authoring; introducing an editor requires updating docs and structured regression coverage.

## Proposed UI & State Flow

### Edit Mode Activation & Layout
- Detect edit mode via `const EDIT_MODE = qbool('edit');` in `public/index.js`. When true:
  - Mount a new fixed right-rail panel (`public/ui/timelineEditorPanel.js`) inside `#app`.
  - Mount a bottom status HUD (`public/ui/timelineEditorStatus.js`) that mirrors the selected frame and selection metadata.
- Expose editor helpers on `viewerApi.timelineEditor` (for tests) with methods like `getMessages()`, `setMessage(id, draft)`, `commit()`, `getPlaybackConfig()`.

### Timeline Control Editor (Right Rail)
- Panel layout:
  1. **Playback Preset** section (loop, auto-play, default FPS).
  2. **Control Messages** list with add/remove buttons.
  3. **Message Editor** detail pane.
- Messages list:
  - Displays `label || id` plus active range summary (`frame-0010 → frame-0040`).
  - Selecting an item loads its draft into the detail form.
  - `Add` creates a draft with generated id (`cm-${timestamp}`), default priority 0, range spanning entire buffer.
  - `Remove` prompts confirmation and deletes from engine.
- Draft form fields:
  - **Metadata:** id (read-only), label (new optional string), priority (number).
  - **Range:** start/end inputs offering three entry modes:
    - `Frame ID` (validated against `frameBuffer.resolveOffset`).
    - `Offset` (numeric, e.g., `-42`).
    - `Index` (0-based).
    - Each input row includes "Use current frame" button to copy `timelineState.offset`.
  - **Actions:** toggles for each supported action. Toggling off removes the action from the message.

### Action-Specific Editors
- **Playback Speed (`timeline.playbackSpeed`):**
  - Radio buttons for `FPS` or `Speed multiplier`.
  - Numeric input with validation (≥1 fps, >0 multiplier).
  - Optional easing dropdown (`linear`, `ease-in`, `ease-out`) and transition duration slider (0–1000 ms).
  - Preview chip shows resulting effective FPS.
- **Callout (`overlay.callout`):**
  - Multiline text area (supports `\n` line breaks).
  - Anchor selector (`world`, `atom`, `bond`):
    - `world`: XYZ inputs (Angstrom), "Use camera target" button.
    - `atom`: "Use current atom selection" button; displays selected index/element fallback when no selection.
    - `bond`: "Use current bond selection" button plus orientation radio.
  - Panel size inputs (width/height in Å) and text size slider (0.1–2.0).
  - Style editor (background/color inputs with alpha, optional border toggle).
- **Opacity Focus (`visual.opacityFocus`):**
  - Focus atoms input (comma-separated indices) with "Use current selection" button (atom adds single index, bond adds both).
  - `includeBonds` dropdown (`none`, `connected`, `exact`).
  - Sliders (0–1) for focus/background atom and bond opacities.
  - Transition duration slider plus "Preview mask" button (applies without closing form).

All forms validate inputs as the user types; invalid fields show inline errors and disable `Save`.

### Playback Configuration Editor
- Located at top of the panel; binds directly to `timelinePlayback`.
- Controls:
  - `Default FPS` number input (10–120).
  - `Auto play` toggle.
  - `Loop` toggle; when enabled, range inputs mirror message range widgets (with "Use current frame" convenience).
  - `Start frame` ref; "Use oldest" and "Use current" buttons.
- On change: call `timelinePlayback.setBaseConfig` and `recomputePlaybackRuntime()` followed by `timelinePlayback.shouldAutoPlay()` check; if auto-play newly enabled, prompt to preview (optional) but do not auto-start until user clicks play.

### Selection & Frame HUD (Bottom Bar)
- Status bar anchored bottom-center:
  - Left chunk: `Frame` display showing `Live` or `frame-0123 (offset -5)`.
  - Right chunk: `Selection` display:
    - Atom: `Atom 14 (O)` plus XYZ (rounded).
    - Bond: `Bond 12–18 (C–O)` plus orientation label.
    - None: `No selection`.
- HUD updates on:
  - `applyTimelineFrame` (frame text).
  - `state.bus.on('selectionChanged')` (selection text).
- Provide `viewerApi.timelineEditor.getStatus()` for Playwright asserts.

### Immediate State Synchronisation
- `Save` button serialises the draft list, updates `controlEngine.setMessages`, calls `controlEngine.refresh()`, and re-applies `applyTimelineFrame(timelineState.offset)` to show effects.
- After updates, call `sessionStateManager.setBaselineFromState({ kind: 'json', label: 'edit' }, { includeTimeline: true })` so subsequent exports include the edits.
- Removing a message also re-applies the current frame; if the deleted message owned the active callout/opacity, `applyControlActions` clears the visuals.

### SessionStateManager & Engine Updates
- Bump schema to **v4** (`SCHEMA_VERSION = 4`).
- Extend control-message snapshots with optional `label`, `notes`, and action-specific metadata pass-through (already supported, but ensure `sanitizeAction` retains new fields such as style colors).
- When loading v3 snapshots (fixtures), upgrade in place:
  - Inject default labels (id) and convert existing playback config.
  - Mark timeline playback defaults if missing.
- Add `timeline.playback.editor` metadata to record defaults (used by UI but optional).

### Supporting Changes
- `controlMessageEngine`:
  - Preserve `label`/`notes` in `setMessages`/`getSnapshot`.
  - Expose `getMessages()` to ease editor sync.
- `timelinePlaybackController`:
  - Add `setAutoPlay(boolean)` & `setLoopRange(range)` helpers invoked by editor.
  - Emit `onPlaybackStateChange` callbacks when base config changes (for UI toggles).
- `viewerApi` additions:
  - `viewerApi.timelineEditor = { list(), select(id), setDraft(draft), commit(), remove(id), setPlayback(cfg) }`.
  - `viewerApi.timeline.getFrameMeta(offset)` to aid tests.

## JSON Session Format (schemaVersion 4)

### Top-Level Structure
```jsonc
{
  "schemaVersion": 4,
  "savedAt": "2025-03-01T18:42:33.410Z",
  "source": { "kind": "json", "label": "SN2 control walkthrough" },
  "viewer": { ... },
  "energyPlot": { ... },
  "timeline": { ... },
  "websocket": { ... }
}
```

### `viewer`
- `elements`: array of element symbols (strings).
- `positions`: `[[x,y,z], ...]` in Å.
- `velocities`: optional; same shape as positions.
- `cell`: `{ a:[3], b:[3], c:[3], originOffset:[3], enabled:bool }`.
- `showCell`, `showGhostCells`: booleans.

### `energyPlot`
- `series`: array of `{ energy: number, kind: 'idle'|'md'|'relax', label?: string }`.
- `markerIndex`: nullable integer identifying highlighted point.

### `timeline`
- `capacity`: integer (buffer size).
- `frames`: oldest → newest array with entries:
  - `id`: string (`frame-00042`).
  - `numericId`: integer (optional mirror of `id`).
  - `kind`: `'idle'|'md'|'relax'`.
  - `seq`, `simStep`, `userInteractionCount`, `timestamp`.
  - `energy`, `temperature`, `energyIndex`.
  - `positions`, `velocities`, `forces` (triple arrays).
- `playback`:
  - `defaultFps`: integer.
  - `autoPlay`: boolean.
  - `loop`: boolean.
  - `startFrame`: `{ frameId?, frameIndex?, offset? }`.
  - `loopRange`: `{ startFrameId?, endFrameId?, start?: {...}, end?: {...} }`.
- `controlMessages`: array of messages (see below).
- `lastLiveMode`: `'idle'|'md'|'relax'`.
- `wasRunning`: boolean.
- `pendingSimParams`: MD/relax parameter snapshot when `wasRunning=true`.

### Control Message Format
```jsonc
{
  "id": "focus-approach",
  "label": "Approach phase",
  "priority": 10,
  "notes": "Slow motion over nucleophilic attack",
  "range": {
    "start": { "frameId": "frame-0080" },
    "end": { "frameId": "frame-0120", "inclusive": true }
  },
  "actions": [
    {
      "type": "timeline.playbackSpeed",
      "fps": 12,
      "transitionMs": 150,
      "easing": "ease-in-out"
    },
    {
      "type": "overlay.callout",
      "text": "Approach phase\nNucleophile forming bond",
      "anchor": { "mode": "bond", "atoms": [2, 15], "orientation": "midpoint" },
      "panelSize": { "width": 1.2, "height": 0.6 },
      "textSize": 0.5,
      "style": { "background": "rgba(12,18,32,0.85)", "color": "#f0f7ff" },
      "offset": [0, 0.8, 0]
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
```

### `websocket`
- `nextSeq`, `lastAck`, `userInteractionCount`, `totalInteractionCount`, `simStep`.
- Legacy aliases (`seq`, `clientAck`) accepted on load.

## Implementation Plan
1. **Editor scaffolding**
   - Create `timelineEditorPanel.js` and `timelineEditorStatus.js` with mount/unmount helpers.
   - Detect `?edit=1` in `public/index.js` and initialise the editor with injected dependencies.
2. **Control message data layer**
   - Extend `controlMessageEngine` to store labels/notes and expose `list()`/`upsert()`.
   - Add draft handling + validation utilities (`timelineEditorStore.js`).
3. **Playback config bridging**
   - Add setter/getter wrappers on `timelinePlayback` and integrate `recomputePlaybackRuntime` hooks.
4. **HUD integration**
   - Update `applyTimelineFrame` and selection listeners to feed the status HUD.
   - Ensure timeline resume clears HUD (set to `Live`).
5. **Session manager updates**
   - Bump schema to 4; add migration from v3 snapshots.
   - Ensure `captureSnapshot` persists new metadata and editor commits flush baseline.
6. **Viewer API + tests hooks**
   - Expose editor helpers under `viewerApi.timelineEditor`.
   - Document new hooks in `test_hooks.md`.
7. **Documentation refresh**
   - Update `timeline.md`, `testing.md`, and `frontend_design.md` sections describing authoring and schema v4.
8. **Fixtures**
   - Prepare plan to migrate `public/examples/sn2/sn2.json` (executed after code approval per Step 7).

## Regression Testing Strategy

### Jest (unit)
- `tests/timelineEditor.store.spec.js`: draft validation, range conversions, action toggles.
- `tests/controlMessageEngine.spec.js`: extend to cover label preservation and overlapping priority updates triggered by editor saves.
- `tests/sessionStateManager.spec.js`: verify schema v4 round-trip, editor metadata retained, v3 upgrade path.
- `tests/timelinePlaybackController.spec.js`: ensure new setters update snapshot + runtime.

### Playwright (E2E)
- `ws-timeline-editor-basic.spec.js`: launch with `?edit=1`, add message, set range via "Use current frame", verify callout/opacity apply immediately.
- `ws-timeline-editor-playback.spec.js`: tweak playback to loop subset, play through range, confirm loop stops at end.
- `ws-timeline-editor-selection.spec.js`: select atom/bond, assert HUD shows expected text and "Use selection" buttons populate form fields.
- Extend `ws-session-save-load.spec.js` to ensure edited messages + playback persist after JSON export/import.

### Backend (pytest)
- `test_snapshot_schema_v4.py`: load v4 session, confirm counters + control messages are ignored by backend (no proto changes) and streaming resumes.
- Upgrade helper test ensuring v3 → v4 conversion leaves control messages intact.

## Options & Trade-offs
1. **Dedicated editor module (recommended)**
   - **Pros:** Keeps `index.js` manageable; encapsulates UI/validation; reusable in VR/desktop contexts; clear test surface.
   - **Cons:** Introduces new modules and cross-component wiring to maintain.
2. **Extend existing desktop panel**
   - **Pros:** Fewer files; leverages existing panel styles.
   - **Cons:** Overloads left panel, complicates non-edit flows, and mixing authoring controls with run controls risks accidental edits.

## Risks & Follow-Ups
- Need careful validation to prevent references to non-existent frame IDs after trimming buffer; plan to disable edits when buffer underflows and prompt to reselect range.
- Looping + auto-play toggles may surprise users if saved sessions auto-start; include confirmation modal before enabling auto-play in editor.
- VR/XR parity: editor is desktop-only, but control messages/callouts must still render in VR—tests should include a VR smoke run post-implementation.
- After implementation approval, convert `public/examples/sn2/sn2.json` to schema v4 with curated control messages (Step 7).

Pending your approval, we will proceed with the implementation (Step 6) and subsequently migrate `sn2.json`.
