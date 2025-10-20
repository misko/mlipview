## Backend geometry debug flags

When diagnosing UMA prediction discrepancies, enable verbose backend geometry logging to capture positions, cell, atomic numbers, dataset, energies, and forces around each UMA predict call.

### Quick reference

- `UMA_GEOM_DEBUG` (backend server env)
  - Scope: Enables geometry debug logs inside the UMA predictor replica (Ray Serve deployment).
  - Effect: Prints one or more lines per predict call with inputs and outputs.
  - Usage: Set in the environment of the backend process: `UMA_GEOM_DEBUG=1`.

- `BACKEND_DEBUG_GEOM` (Playwright harness env)
  - Scope: Convenience switch for e2e runs; the Playwright global setup forwards this to the backend so you don’t have to export `UMA_GEOM_DEBUG` manually.
  - Effect: Also ensures backend stdout/stderr are captured and tailed at the end of the test.
  - Usage: Prefix your Playwright run: `BACKEND_DEBUG_GEOM=1 npx playwright test ...`

### What gets logged

When enabled, the UMA replica emits lines like the following:

- Input snapshot per item:
  - `[UMA][geom] item=0 natoms=3 pos=[[...],[...],[...]] cell=[[...],[...],[...]] Z=[1, 1, 8] dataset=['omol']`

- Output snapshot per item:
  - `[UMA][geom][out] item=0 E=[-2079.8650...] forces=[[...], [...], [...]]`

You may also see timing and device lines useful for performance checks:
- `[UMA] predict called on device cuda size=1`
- `[UMA] predict finished size=1 wall=0.0267s ms/item=26.72`

### Where the logs go

- Playwright e2e: backend logs are piped to `test-ws-e2e.log`. The `tests-e2e/waterRelaxation.spec.js` test prints the tail of this file upon completion to aid debugging in CI.
- Local backend runs: logs print to the backend process stdout/stderr (your terminal or process manager log).

### How to enable

- Playwright (recommended for e2e):
  - Bash:
    - `BACKEND_DEBUG_GEOM=1 npx playwright test tests-e2e/waterRelaxation.spec.js`

- Direct backend launch (outside Playwright):
  - `UMA_GEOM_DEBUG=1 ./mlipview_venv/bin/python -m fairchem_local_server2.serve_ws_app --ngpus 1 --ncpus 2 --nhttp 1 --http-port 8000`

### Troubleshooting tips

- Atomic numbers are zeros (e.g., `Z=[0, 0, 0]`):
  - Cause: Client did not send `atomic_numbers` during initialization.
  - Fix: Ensure the UI path sends `atomic_numbers` (the viewer derives them via `elementToZ` in `public/index.js`), or include them explicitly in `userInteraction({ atomic_numbers, positions, ... })`.

- Backend rejects a request with `[INVALID_ATOM_Z] atomic_numbers must be positive integers; received ...`:
  - Cause: One or more atomic numbers are non-positive.
  - Fix: Correct the `atomic_numbers` in your init; the backend now enforces this in both the UMA predictor path and the ASE worker pool.

- Logs are too noisy / large:
  - Enable only for short, targeted tests; disable for load or long streaming runs.

Notes:
- Geometry debug printing is intended for development and CI diagnostics. It adds I/O but does not change numerical results.
- If you need to grep for a specific field in the Playwright log, e.g., to verify species wiring:
  - `grep -n "Z=\[1, 1, 8\]" test-ws-e2e.log || true`

# Testing guide for the WebSocket + protobuf migration

This document captures the updated testing patterns after migrating to a protobuf-only WebSocket API and cleaning up the frontend test hooks.

## What changed (test-facing)

- REST flows and SIMPLE_CALCULATE/INIT_SYSTEM are removed. Use USER_INTERACTION for initialization and idle computes; use START_SIMULATION/STOP_SIMULATION or requestSingleStep for MD/relax.
- WebSocket client test hooks are simplified:
  - Prefer per-instance test hook: `ws.setTestHook(fn)` (new)
  - Inject decoded frames without protobuf: `ws.injectTestResult(obj)` (new)
  - Single global fallback for legacy tests: `globalThis.__WS_TEST_HOOK__` and `globalThis.__ON_WS_RESULT__`
- Frontend now uses protobuf ESM stubs under `public/proto/fairchem_local_server2/session_pb.js`.

## Frontend: WS client API recap

- Initialization: `ws.userInteraction({ atomic_numbers, positions, velocities?, cell? })`
- Idle compute (forces/energy): send USER_INTERACTION with positions, then `await ws.waitForEnergy()`
- Single simulation step: `await ws.requestSingleStep({ type: 'md' | 'relax', params })`
- Streaming runs: `ws.startSimulation({ type: 'md'|'relax', params }); ws.stopSimulation()`
- Counters for correlation/backpressure: `ws.setCounters({ userInteractionCount, simStep })` and `ws.ack(seq)`

## Test hooks and utilities

- Instance-level (preferred):
  - `ws.setTestHook(fn)` — receives a minimal JSON snapshot of each outgoing ClientAction (type, counters, sim params, counts)
  - `ws.injectTestResult(obj)` — directly fan out a decoded ServerResult-like object to all registered listeners
- Global fallback (legacy support only):
  - `globalThis.__WS_TEST_HOOK__ = fn` — used if no instance hook is set
  - `globalThis.__ON_WS_RESULT__ = fn` — inject decoded frames globally; created automatically when the WS client boots
- Debug toggles:
  - WS logging: `window.__WS_DEBUG_ENABLE__(true)` or `?wsDebug=1`
  - API logging in the app: `window.__MLIPVIEW_DEBUG_API = true` or `?debug=1`

## Jest examples (jsdom/node)

### 1) Assert an outgoing USER_INTERACTION on init

```js
import { getWS } from '../public/fairchem_ws_client.js';

it('sends USER_INTERACTION init with atoms + positions', async () => {
  const ws = getWS();
  const sent = [];
  ws.setTestHook(msg => sent.push(msg));

  // Simulate app init path
  const atomic_numbers = [8, 1, 1];
  const positions = [[0,0,0],[0.96,0,0],[-0.24,0.93,0]];
  ws.userInteraction({ atomic_numbers, positions });

  // Minimal assertion on type and counts
  expect(sent.length).toBeGreaterThan(0);
  const last = sent[sent.length - 1];
  expect(last.type).toBeDefined();
  // If you import the enum, assert equality to ClientAction_Type.USER_INTERACTION
  expect(last.positionsCount).toBe(positions.length);
});
```

### 2) Drive idle compute via waitForEnergy

```js
import { getWS } from '../public/fairchem_ws_client.js';

it('resolves idle compute energy via waitForEnergy', async () => {
  const ws = getWS();
  // Inject an immediate frame to avoid standing up a backend in unit tests
  const p = ws.waitForEnergy({ timeoutMs: 100 });
  ws.injectTestResult({ energy: -1.23, forces: [[0,0,0]], userInteractionCount: 1, simStep: 0 });
  const { energy, forces } = await p;
  expect(energy).toBeCloseTo(-1.23, 6);
  expect(Array.isArray(forces)).toBe(true);
});
```

### 3) Single-step relax

```js
import { getWS } from '../public/fairchem_ws_client.js';

it('requestSingleStep relax resolves with frame', async () => {
  const ws = getWS();
  // Simulate a backend response injected by the test harness
  const resP = new Promise(resolve => {
    ws.onResult(frame => resolve(frame));
  });
  setTimeout(() => {
    ws.injectTestResult({
      positions: [[0,0,0]],
      forces: [[0,0,0]],
      energy: -2.0,
      simStep: 1,
    });
  }, 10);

  const r = await ws.requestSingleStep({ type: 'relax', params: { calculator: 'lj', fmax: 0.1, max_step: 0.2, optimizer: 'bfgs' } });
  expect(typeof r.energy).toBe('number');
  await expect(resP).resolves.toBeDefined();
});
```

Notes:
- In pure unit tests, prefer `injectTestResult()` over standing up the server.
- Use `setTestHook()` to assert the shape and routing of outgoing messages, including counters.

## Playwright examples (browser/e2e)

### 1) Observe outgoing actions

```ts
// In page context
await page.exposeFunction('__REC', (msg) => {
  // record messages in the test runner via page.exposeFunction
});
await page.evaluate(() => {
  const ws = window.__fairchem_ws__;
  ws.setTestHook(window.__REC);
  // Kick off a user interaction init
  ws.userInteraction({ atomic_numbers: [8,1,1], positions: [[0,0,0],[0.96,0,0],[-0.24,0.93,0]] });
});
```

### 2) Inject a frame and wait for energy

```ts
await page.evaluate(async () => {
  const ws = window.__fairchem_ws__;
  const p = ws.waitForEnergy({ timeoutMs: 1000 });
  ws.injectTestResult({ energy: -1.0, forces: [[0,0,0]], simStep: 0 });
  const r = await p;
  // r.energy is -1.0
});
```

### 3) Streamed runs

```ts
await page.evaluate(() => {
  const ws = window.__fairchem_ws__;
  ws.startSimulation({ type: 'md', params: { calculator: 'uma', temperature: 300, timestep_fs: 1.0, friction: 0.02 } });
});
// Optionally listen for frames via global fallback
await page.evaluate(() => {
  window.__ON_WS_RESULT__ = (obj) => { /* collect frames */ };
});
// Stop later:
await page.evaluate(() => window.__fairchem_ws__.stopSimulation());
```

## Backend tests (Python) — update notes

- The protobuf schema removed legacy types. Replace any INIT_SYSTEM/SIMPLE_CALCULATE usage:
  - Initialize by sending USER_INTERACTION with atomic_numbers + positions (and optional cell/velocities).
  - To trigger idle compute, send a USER_INTERACTION with positions; the server responds with energy/forces and omits positions.
  - For single steps, use the MD/RELAX pathways (the worker pool provides `run_md`, `run_relax`, and `run_simple` internally).
- WebSocket is protobuf-only. Text frames can be ignored; ACKs are sent as normal ClientAction messages with type=PING and ack set.

## Common pitfalls and tips

- Counters: set `ws.setCounters({ userInteractionCount, simStep })` before `startSimulation` to keep server-side echo aligned for gating.
- Stale-frame gating: if your UI discards frames based on `userInteractionCount`, make sure tests set the same count before injecting or asserting.
- Debugging WS:
  - Enable noisy logs with `window.__WS_DEBUG_ENABLE__(true)` or `?wsDebug=1` in the test URL.
  - For API-level logs in the viewer code, use `window.__MLIPVIEW_DEBUG_API = true` or `?debug=1`.

## Migration checklist for existing tests

- [ ] Remove any references to INIT_SYSTEM or SIMPLE_CALCULATE
- [ ] Replace init flows with `ws.userInteraction({ atomic_numbers, positions, velocities?, cell? })`
- [ ] For idle force/energy, use `ws.waitForEnergy()` after sending `userInteraction({ positions })`
- [ ] Prefer `ws.setTestHook` for asserting outgoing actions
- [ ] Prefer `ws.injectTestResult` or `globalThis.__ON_WS_RESULT__` for simulating incoming frames
- [ ] If evaluating gating, set `ws.setCounters({ userInteractionCount, simStep })` on the client prior to sending `startSimulation`

---
If you find any gaps or edge cases not covered here, add an example to this file so we can keep tests consistent with the new protocol.
