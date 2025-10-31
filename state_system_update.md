# State System Update Proposal

## Implementation Status (current branch)
- `SessionStateManager` now owns snapshot capture/restore, JSON persistence, and reset plumbing. It resets interaction counters and reseeds WS sequencing before uploading a fresh `full_update`, eliminating stalled resumes after loading JSON timelines.
- Timeline resume flow logs `[Timeline][resumeLiveFromTimeline] …` events and replays at most the requested burst of frames before restarting continuous MD/relax loops.
- Playwright coverage (`ws-session-playback-resume.spec.js`) runs the full JSON → timeline playback → live resume cycle and bounds live-frame replay (≤170) to prevent runaway simulations during CI.
- Backend ingress exposes `MLIPVIEW_RESUME_DEBUG` for selective ACK/seq tracing while diagnosing resume issues.

## Current Architecture Findings
- The viewer is orchestrated from `public/index.js`, which owns the Babylon scene, energy plot, WebSocket client, timeline buffer, and reset baseline (`seedResetBaseline`) without an external coordinator (`public/index.js:356-455`, `public/index.js:520-581`).  
- Timeline replay is purely client-side: frames are cached in a circular buffer as `Float32Array`s and replayed by the viewer while interactions are locked (`public/core/frameBuffer.js:82-157`, `timeline.md:11-49`).  
- Interaction counters and latch state previously lived in an ad-hoc singleton, while the WebSocket client owned sequencing/ACK bookkeeping with no way to seed values from persisted state (`public/fairchem_ws_client.js:170-242`).  
- Backend session state mirrors client counters and tracks `server_seq`, `client_seq`, pending full updates, and simulation parameters (`fairchem_local_server2/session_models.py:27-118`).  
- XYZ/SMILES imports feed through `applyParsedToViewer`, which refreshes the reset baseline and toggles periodic-cell flags (`public/util/moleculeLoader.js:90-178`). Existing docs emphasise a 500-frame timeline and resume semantics (`README.md:3-9`, `frontend_design.md:25-48`, `backend_design.md:11-33`, `timeline.md:5-24`).  
- Regression coverage already exercises timeline controls, energy markers, and molecule loaders (`testing.md:40-75`), providing scaffolding for new persistence tests.

## Observed Pain Points
- State is fragmented: reset baselines, timeline cache, energy series, and WS counters live in private closures with no serialization hooks (`public/index.js:356-455`, `public/index.js:520-581`).  
- `frameBuffer` stores typed arrays but exposes no export/import helpers, making timeline persistence ad-hoc (`public/core/frameBuffer.js:82-195`).  
- The energy plot is write-only; markers and series cannot be reconstructed after reload, preventing faithful timeline playback.  
- The reset button depends on `state.__resetBaseline`, which is only updated by loaders and manual calls—JSON timeline loads would need to touch multiple internals to stay consistent.  
- WebSocket sequencing cannot be reseeded; loading a saved session risks desynchronised `seq`/`ack` values and stalled streams (`public/fairchem_ws_client.js:170-176`, `fairchem_local_server2/session_models.py:99-104`).  
- Backend has no notion of persisted timelines, so full snapshots must preserve `server_seq`/`client_ack` alignment to avoid triggering the server’s backpressure guard (`backend_design.md:21-33`).

## Proposed Unified Session State Layer
Create a `SessionStateManager` module under `public/core/sessionStateManager.ts` that centralises “authoritative” client session data and exposes a stable API consumed by `index.js`, loaders, and UI.

### Responsibilities
1. **Snapshot management** – capture, serialize, and restore molecule geometry, dynamics, energy history, and timeline frames.  
2. **Back-end sync** – coordinate `full_update` uploads, seeding WS sequence/counter state before replaying or resuming simulations.  
3. **Reset binding** – surface `resetToLastLoad()` so the UI reset button simply calls the manager, regardless of whether the last load came from XYZ, SMILES, or JSON.  
4. **Timeline orchestration** – provide replay iterators and export/import helpers that wrap the existing `frameBuffer`.  
5. **Energy trace persistence** – store the energy plot series + marker in plain arrays suitable for JSON.

### API Surface (draft)
```ts
type SourceKind = 'xyz' | 'smiles' | 'json';
type SessionSnapshot = {
  schemaVersion: 1;
  savedAt: string;
  source: { kind: SourceKind; label?: string; params?: Record<string, unknown> };
  viewer: {
    elements: string[];
    positions: number[][];
    velocities?: number[][];
    cell?: number[][];
    showCell: boolean;
    showGhostCells: boolean;
  };
  energyPlot: { series: Array<{ energy: number; kind: 'idle'|'md'|'relax'; label?: string }>; markerIndex?: number|null };
  timeline: {
    capacity: number;
    frames: Array<{
      kind: 'idle'|'md'|'relax';
      seq: number;
      simStep?: number|null;
      userInteractionCount?: number|null;
      energy?: number|null;
      temperature?: number|null;
      timestamp: number;
      energyIndex?: number|null;
      positions: number[][];
      velocities?: number[][];
      forces?: number[][];
      stress?: number[][];
    }>;
    lastLiveMode: 'idle'|'md'|'relax';
    wasRunning: boolean;
    pendingSimParams?: Record<string, unknown>;
  };
  websocket: { nextSeq: number; lastAck: number; userInteractionCount: number; simStep: number };
};
```

Key helpers:
- `SessionStateManager.capture({ source })` – harvest current viewer + timeline + counters from `index.js` and return `SessionSnapshot`.  
- `SessionStateManager.hydrate(snapshot)` – populate state, reset energy plot/timeline UI, seed WS client, and emit the final frame as a full snapshot.  
- `SessionStateManager.loadXYZ(text|parsed)` – existing loaders delegate here; manager updates `lastSource`, reset baseline, and triggers `ensureWsInit()`.  
- `SessionStateManager.saveToFile()` / `loadFromFile()` – use the above helpers while handling file dialogs.  
- `SessionStateManager.resetToLastSource()` – reapply the manager’s cached snapshot (driven by reset button).

### Front-end Refactors
- Factor the energy plot closure into `energyPlotStore` with `getSeries()`, `setSeries(series, marker?)`, and `reset()` so it can be serialized and restored (`public/index.js:356-455`).  
- Wrap `createFrameBuffer` with an adapter exposing `exportFrames()` and `importFrames(frames)` that internally writes to the existing buffer (`public/core/frameBuffer.js:82-195`).  
- Replace direct touches of `state.__resetBaseline` with `sessionStateManager.setResetBaseline(snapshot)` so JSON loads and XYZ loads share the same path (`public/index.js:520-552`).  
- Extend the WS client with `seedSequencing({ nextSeq, lastSeq?, ack, userInteractionCount, simStep })` (idempotent setter + bounds checks) and expose it through `viewerApi` so JSON loads can align counters before `ensureWsInit()` resumes streaming (`public/fairchem_ws_client.js:170-196`).

### Backend Coordination
- Persist `server_seq` and `client_ack` in the JSON. On hydrate, call the new `seedSequencing` and immediately send a `userInteraction` full snapshot with the saved `nextSeq` to reset the server’s state. The backend already resets cached forces when it sees `full_update=true` (`backend_design.md:28-33`), so no schema change is required.  
- When the saved session indicates `wasRunning=true`, queue `startSimulation` after local timeline playback finishes, using the saved simulation parameters (temperature, friction, relax tolerances) captured from `SessionSnapshot.timeline.pendingSimParams`.

## JSON Session Format
Example payload (pretty-printed for readability):

```json
{
  "schemaVersion": 1,
  "savedAt": "2025-02-12T18:05:42.611Z",
  "source": { "kind": "xyz", "label": "molecules/roy.xyz" },
  "viewer": {
    "elements": ["C", "H", "O"],
    "positions": [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
    "velocities": [[0, 0, 0], [0, 0, 0]],
    "cell": [[8.6, 0, 0], [0, 8.6, 0], [0, 0, 8.6]],
    "showCell": true,
    "showGhostCells": true
  },
  "energyPlot": {
    "series": [
      { "energy": -152.34, "kind": "idle" },
      { "energy": -153.01, "kind": "md" }
    ],
    "markerIndex": 1
  },
  "timeline": {
    "capacity": 500,
    "frames": [
      {
        "kind": "md",
        "seq": 412,
        "simStep": 98,
        "userInteractionCount": 37,
        "energy": -153.01,
        "temperature": 305.4,
        "timestamp": 1739383542611,
        "energyIndex": 1,
        "positions": [[0.01, -0.02, 0.0], [1.23, -0.04, 0.01]],
        "velocities": [[0.02, 0.01, 0], [0.03, 0, 0.01]],
        "forces": [[0.1, -0.1, 0], [-0.1, 0.1, 0]]
      }
    ],
    "lastLiveMode": "md",
    "wasRunning": true,
    "pendingSimParams": { "temperature": 305.4, "friction": 0.02, "timestep_fs": 1 }
  },
  "websocket": {
    "nextSeq": 413,
    "lastAck": 409,
    "userInteractionCount": 37,
    "simStep": 98
  }
}
```

## Implementation Plan
1. **State manager scaffolding** – introduce `SessionStateManager`, factor energy plot/frame buffer adapters, and expose serialization hooks through `viewerApi`.  
2. **Loader integration** – update `applyParsedToViewer` and XYZ/SMILES loaders to call `SessionStateManager.loadXYZ` (or `.loadFromSource('smiles')`), ensuring reset baselines and last-source metadata stay in sync (`public/util/moleculeLoader.js:90-200`).  
3. **JSON persistence UI** – add `saveTimeline()` / `loadTimeline()` entry points (e.g., desktop panel buttons, CLI hook) that call `SessionStateManager.saveToFile()` and `.loadFromFile()`.  
4. **WebSocket seeding** – extend `fairchem_ws_client` with `seedSequencing` (accepting `nextSeq`, optional `lastSeq`, ACK + counters) and ensure `ensureWsInit` reuses the seeded counters before sending `full_update`.  
5. **Timeline replay alignment** – when hydrating, push the most recent frame to the viewer via `applyFramePayload(..., { forceAll: true, allowEnergyPlot: true })`, set the energy marker, and conditionally queue `startContinuous` after playback if `wasRunning` was recorded.  
6. **Reset button hookup** – swap the current reset implementation to `SessionStateManager.resetToLastSource()` so it works uniformly for XYZ/SMILES/JSON loads (`public/index.js:2265-2284`, `public/ui/desktopPanel.js:1686-1714`).  
7. **Server parity checks** – verify that rehydrated sessions keep `server_seq - client_ack < max_unacked` by watching Playwright logs and adjust if a preflight ACK is needed.

## Regression Testing Plan
- **Jest unit tests**  
  - `tests/sessionStateManager.spec.js`: round-trip a snapshot (capture → hydrate) and assert geometry, energy series, and timeline frames survive.  
  - `tests/frameBuffer.serialize.spec.js`: verify `exportFrames()`/`importFrames()` preserve signatures and metadata (`public/core/frameBuffer.js`).  
  - `tests/fairchem_ws_seed.spec.js`: seed sequencing, send a `userInteraction`, and confirm the emitted `seq` matches the seed.
- **Playwright**  
  - `ws-timeline-json-load.spec.js`: save a session, reload via JSON, scrub timeline, and resume live streaming without duplicate frames.  
  - `ws-reset-after-json.spec.js`: load JSON, hit reset, and ensure state reverts to the saved frame + energy marker.  
  - Extend an existing MD loop test to assert that `startSimulation` resumes automatically when `wasRunning=true`.
- **Backend (pytest)**  
  - Add a parity test that sends a `full_update` with pre-seeded `seq`/`ack`, ensuring the server accepts the snapshot and continues streaming (`fairchem_local_server2/tests_py`).  
  - Validate that `server_seq` resets when a new client connects with higher `seq` numbers to guard against stale JSON imports.

## Options & Trade-offs
1. **Unified manager (recommended)**  
   - *Pros:* Single entry point for all state mutations; consistent reset behaviour; JSON, XYZ, and SMILES loads share one code path; straightforward to extend with future state (e.g., annotations).  
   - *Cons:* Requires refactoring `index.js` to depend on the manager, touching multiple subsystems at once.
2. **Minimal persistence helpers**  
   - *Pros:* Quicker to implement (serialize `frameBuffer` + energy plot ad-hoc, patch reset button).  
   - *Cons:* Persists fragmentation, harder to test, and JSON loads would still poke private fields in `index.js`, increasing regression risk when refactoring timeline or reset logic later.

The unified manager satisfies the release goals (save/load timelines, align backend sequencing, reset integration) while setting a foundation for VR/AR session portability. Pending approval, we can proceed with the refactor and corresponding tests in Step (6).
