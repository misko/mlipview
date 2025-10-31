# Test Hook Reference

To support reliable Playwright coverage while keeping the production viewer encapsulated, add reusable debug hooks behind `viewerApi` (and document them here when you introduce or update one).

The hooks listed below are **only** for test code. Guard them behind environment checks (`window.__MLIPVIEW_TEST_MODE`, etc.) so they are harmless in production builds.

---

## `viewerApi.debugGhostSnapshot()`

**Purpose:** Surface ghost/periodic bonding details without reaching into `view._internals`.

**Suggested shape:**

```js
viewerApi.debugGhostSnapshot = () => ({
  ghostBondCount: view._internals?.ghostBondGroups
    ? [...view._internals.ghostBondGroups.values()]
        .reduce((sum, grp) => sum + grp.mats.length, 0)
    : 0,
  ghostGroups: view._internals
    ? [...view._internals.ghostBondGroups.entries()].map(([key, grp]) => ({
        key,
        count: grp.mats.length,
      }))
    : [],
});
```

**Used by:** Benzene rotation E2E test—asserts ghost bonds render after rotations.

---

## `viewerApi.debugHighlightState()`

**Purpose:** Query highlight visibility without poking private flags.

**Suggested shape:**

```js
viewerApi.debugHighlightState = () => ({
  atomVisible: !!view._internals?.highlight?.atom?.isVisible,
  bondVisible: !!view._internals?.highlight?.bond?.isVisible,
});
```

**Used by:** Bond-highlight Playwright test (atom vs. bond selection UX).

---

## `viewerApi.debugStreamListenerStats({ reset = false } = {})`

**Purpose:** Surface the attach/detach counts for the MD/relax streaming websocket listeners so tests can assert that repeated start/stop cycles do not leak extra handlers.

**Suggested shape:**

```js
const streamListenerStats = {
  md: { attach: 0, detach: 0, active: 0, maxActive: 0 },
  relax: { attach: 0, detach: 0, active: 0, maxActive: 0 },
};

viewerApi.debugStreamListenerStats = ({ reset = false } = {}) => {
  if (reset) {
    for (const key of Object.keys(streamListenerStats)) {
      Object.assign(streamListenerStats[key], { attach: 0, detach: 0, active: 0, maxActive: 0 });
    }
  }
  return {
    md: { ...streamListenerStats.md },
    relax: { ...streamListenerStats.relax },
  };
};
```

**Used by:** `tests-e2e/x-md-multi-start-stop.spec.js` to confirm triple MD start attempts still leave exactly one active listener.

---

## `viewerApi.debugCameraControls()`

**Purpose:** Track camera detach/attach invocations during drag interactions.

**Suggested shape:**

```js
const cameraStats = { detachCalls: 0, attachCalls: 0 };
const origDetach = camera.detachControl.bind(camera);
const origAttach = camera.attachControl.bind(camera);

camera.detachControl = (...args) => {
  cameraStats.detachCalls += 1;
  return origDetach(...args);
};
camera.attachControl = (...args) => {
  cameraStats.attachCalls += 1;
  return origAttach(...args);
};

viewerApi.debugCameraControls = ({ reset = false } = {}) => {
  if (reset) {
    cameraStats.detachCalls = 0;
    cameraStats.attachCalls = 0;
  }
  return { ...cameraStats };
};
```

**Used by:** Camera-suppression Playwright test to confirm detaching during drag and reattaching afterwards.

---

## `viewerApi.debugBondMetrics()`

**Purpose:** Fetch bond lengths and classifications via public API instead of `state.bonds` internals.

**Suggested shape:**

```js
viewerApi.debugBondMetrics = () => {
  const { elements, positions, bonds } = viewerApi.state;
  const distance = (a, b) => {
    const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return Math.hypot(dx, dy, dz);
  };
  return bonds.map(({ i, j, opacity }) => ({
    i,
    j,
    elements: [elements[i], elements[j]],
    length: distance(positions[i], positions[j]),
    opacity,
  }));
};
```

**Used by:** Benzene rotation test to validate C–C / C–H ranges after manipulating geometry.

---

## `viewerApi.setTestAutoSelectFallback(on = true)`

**Purpose:** Enable or disable the automated “select atom 0 when pick misses” fallback without touching `_internals` from tests.

**Suggested shape:**

```js
viewerApi.setTestAutoSelectFallback = (on = true) => {
  if (view._internals) {
    view._internals._debugAutoSelectFirstOnEmpty = !!on;
    return true;
  }
  return false;
};
```

**Used by:** Camera-suppression Playwright test to trigger drag behavior even when the pointer doesn’t land on an atom.

---

## `viewerApi.debugSelectAtom(index)` / `viewerApi.debugSelectBond({ … })` / `viewerApi.debugGetSelection()`

**Purpose:** Drive the selection service without reaching into private structures. `debugSelectAtom` selects an atom by index, `debugSelectBond` accepts `{ i, j, index?, key?, orientation?, side? }`, and `debugGetSelection` returns the current selection snapshot.

**Suggested shape:**

```js
viewerApi.debugSelectAtom = (index) => {
  selectionService.clickAtom(Number(index));
  return selectionService.get();
};

viewerApi.debugSelectBond = ({ i, j, index = 0, key = `${i}-${j}`, orientation, side } = {}) => {
  selectionService.clickBond({ i, j, index, key });
  if (orientation != null && viewerApi.state.selection?.kind === 'bond') {
    viewerApi.state.selection.data.orientation = orientation;
  } else if (side) {
    viewerApi.state.selection.data.orientation = side === 'i' ? 1 : 0;
  }
  return selectionService.get();
};

viewerApi.debugGetSelection = () => selectionService.get();
```

**Used by:** Bond highlight and benzene rotation Playwright tests to set up selections via public hooks.

---

## `viewerApi.timeline`

**Purpose:** Shallow timeline control surface for Playwright tests and debugging tooling.

**Shape:**

```js
viewerApi.timeline = {
  select: (offset) => handleTimelineOffsetRequest(offset),
  play: (offset) => handleTimelinePlayRequest(offset ?? timelineState.offset),
  pause: () => handleTimelinePauseRequest(),
  live: () => handleTimelineLiveRequest(),
  getState: () => ({
    ...(timelineUi?.getState?.() || {}),
    active: !!timelineState.active,
    playing: !!timelineState.playing,
    offset: timelineState.offset,
  }),
  getSignature: (offset) => frameBuffer.getSignature(offset ?? timelineState.offset),
  getOffsets: () => frameBuffer.listOffsets(),
  getFrameMeta: (offset) => {
    const off = Number.isFinite(offset) ? offset : timelineState.offset;
    if (!Number.isFinite(off)) return null;
    const entry = frameBuffer.getByOffset(off);
    if (!entry) return null;
    return {
      frameId: entry.id,
      offset: off,
      frameIndex: offsetToFrameIndex(off),
    };
  },
};
```

**Used by:** Timeline Playwright suites to assert mode transitions, slider selection, energy marker behaviour, and buffer contents without reaching through DOM internals.

---

## `viewerApi.timelineEditor` (edit mode only)

**Purpose:** Drive the timeline control-message editor without poking DOM internals.

**Shape:**

```js
viewerApi.timelineEditor = {
  refresh: () => timelineEditorPanel.refresh(),
  getState: () => timelineEditorPanel.getState(),
  select: (id) => timelineEditorPanel.select(id),
  getDraft: () => timelineEditorPanel.getCurrentDraft(),
  status: () => ({
    frame: editorStatusFrame,
    selection: editorStatusSelection,
  }),
};
```

**Used by:** Authoring Playwright suites (`ws-timeline-editor-*.spec.js`) to add/remove messages, assert playback configuration, and verify the bottom status HUD updates as frames/selection change.

---

## `viewerApi.debugWsState()`

**Purpose:** Surface the recent websocket state transitions (connecting, reconnect scheduling, open/close, etc.) plus the reconnection banner state so tests can assert precise reconnect behaviour without scraping the DOM.

**Suggested shape:**

```js
viewerApi.debugWsState = () => ({
  reconnect: {
    attempts: reconnectState.attempts,
    nextAttemptAt: reconnectState.nextAttemptAt,
    bannerVisible: reconnectState.bannerVisible,
  },
  pendingResume: pendingResume.kind ? { kind: pendingResume.kind, opts: { ...pendingResume.opts } } : null,
  log: wsStateLog.map((evt) => ({ ...evt })),
});
```

**Used by:** `tests-e2e/x-ws-reconnect.spec.js` to wait for banner visibility changes and verify reconnect/open events.

---

## `viewerApi.forceWsReconnect()`

**Purpose:** Immediately trigger a websocket reconnect attempt (cancelling any backoff timer) so tests can simulate the user pressing the “Reconnect now” CTA.

**Suggested shape:**

```js
viewerApi.forceWsReconnect = () => {
  getWS().reconnectNow?.();
  return true;
};
```

**Used by:** Reconnect Playwright tests to assert manual reconnects clear the banner promptly.

---

## `viewerApi.ws.pauseIncoming(ms = 50)`

**Purpose:** Temporarily suppress inbound WebSocket frames so tests can inject synthetic results without the next live frame immediately overwriting them.

**Suggested shape:**

```js
viewerApi.ws.pauseIncoming = (ms = 50) => {
  const ws = getWS();
  if (ws?.pauseIncoming) {
    ws.pauseIncoming(ms);
    return true;
  }
  return false;
};
```

**Used by:** `tests-e2e/x-bond-rotation.spec.js` right before `injectTestResult`, giving the Playwright assertion ~100 ms of quiet time so the injected `+1 Å` deltas remain visible.

---

## `viewerApi.simulateWsDrop({ failAttempts } = {})`

**Purpose:** Forcefully close the websocket while keeping the client in auto-reconnect mode; optionally pre-program the next N reconnect attempts to fail. This lets tests emulate server restarts mid-run.

**Suggested shape:**

```js
viewerApi.simulateWsDrop = ({ failAttempts = 0 } = {}) => {
  if (Number.isFinite(failAttempts)) window.__MLIPVIEW_WS_FAIL_ATTEMPTS = Math.max(0, failAttempts | 0);
  getWS().forceDisconnect?.('test-drop');
  return true;
};
```

**Used by:** `tests-e2e/x-ws-reconnect.spec.js` to drop the connection during an MD run and verify automatic recovery.

---

## Global Debug Toggles

Not every diagnostic warrants its own hook—some are easier to expose via opt-in globals so manual QA can toggle them quickly:

- `window.__MLIP_DEBUG_STRETCH = true` (or `?bondStretchDebug=1` on the URL) enables detailed logging from the bond service and molecule view. Each recompute prints opacity ranges together with the number of bonds routed to the translucent mesh. Reset the flag (or reload without the query) once you have the output you need.

Keep these toggles behind explicit opt-in so automated suites stay quiet by default.

---

### Authoring Notes

- Keep hook implementations colocated with viewer bootstrap (`public/index.js`) so they can access the live `view` and `state`.
- Wrap each hook registration in a guard:
  ```js
  if (typeof window !== 'undefined' && window.__MLIPVIEW_TEST_MODE) {
    Object.defineProperty(window.viewerApi, 'debugHighlightState', { value: ..., configurable: true });
  }
  ```
- Avoid mutating private structures outside of the hook registration; tests should treat the results as readonly snapshots.
- Update this document whenever you add/remove a hook so teams can reuse existing helpers rather than invent new ones.

---

## XR HUD Babylon Harness (`tests/utils/xrHudTestHarness.js`)

**Purpose:** Shared Babylon GUI + scene stubs for XR HUD unit tests so we avoid duplicating mock setup in every spec.

**Exports:**

- `setupXRHudTestEnv()` – installs global `window`, `document`, and `BABYLON` GUI stubs if absent.
- `resetXRHudSingleton()` – clears cached `window.__XR_HUD_FALLBACK` objects between specs.
- `makeHudScene()` – returns the instrumented mock scene with a camera and render observers.
- `runHudFrames(scene, count)` – plays queued frame callbacks to simulate render ticks.
- `findGuiNode(root, predicate)` – BFS helper to locate GUI nodes by id/text.

**Used by:** `tests/x-xr-energy-hud.dom.spec.js`, `tests/x-xr-hud-buttons.spec.js`, `tests/x-xr-hud-energy-depth.spec.js`, `tests/x-xr-hud-energy-plot-picking.spec.js`, `tests/x-xr-hud-energy-scale.spec.js`, `tests/x-xr-hud-position.spec.js`, `tests/x-xr-hud-texture-clamp.spec.js`.
