/** @jest-environment jsdom */
/**
 * Mock browser test: ?autoMD=0 and a user interaction results in multiple plotted energies
 */
import { jest } from '@jest/globals';
import { getWS } from '../public/fairchem_ws_client.js';

describe('autoMD=0 interaction energies', () => {
  beforeEach(() => {
    if (!global.window) global.window = global;
  });

  test('after baseline, a drag interaction yields multiple energy values', async () => {
    window.__MLIPVIEW_TEST_MODE = true;
    window.location = { protocol: 'http:', host: '127.0.0.1:4000', search: '?autoMD=0' };
    window.__MLIPVIEW_SERVER = 'ws://127.0.0.1:8000';
    // Stub WS and respond to each USER_INTERACTION with a distinct energy
    const origWS = global.WebSocket;
    class FakeWS {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen && this.onopen();
        }, 0);
      }
      send() {}
      close() {}
      onopen() {}
      onmessage() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;
    const ws = getWS();
    let __uiCount = 0;
    ws.setTestHook((msg) => {
      try {
        if (msg && msg.type != null && msg.type !== 'INIT_SYSTEM') {
          __uiCount++;
          const seq = 100 + __uiCount;
          const e = -0.97 + __uiCount * 0.01;
          setTimeout(
            () =>
              ws.injectTestResult({
                seq,
                energy: e,
                forces: [
                  [0, 0, 0],
                  [0, 0, 0],
                ],
              }),
            0
          );
        }
      } catch {}
    });

    document.body.innerHTML =
      '<canvas id="energyCanvas" width="200" height="50"></canvas><div id="energyLabel"></div>';
    const viewer = document.createElement('canvas');
    viewer.id = 'viewer';
    viewer.addEventListener = () => {};
    document.body.appendChild(viewer);

    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(viewer, {
      elements: ['H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [],
    });

    // Baseline energy (auto-emitted via hook). Wait until plotted.
    const waitUntil = async (fn, timeoutMs = 200) => {
      const t0 = Date.now();
      while (Date.now() - t0 < timeoutMs) {
        if (fn()) return true;
        await new Promise((r) => setTimeout(r, 5));
      }
      return false;
    };
    const ok = await waitUntil(() => api.debugEnergySeriesLength() >= 1, 300);
    expect(ok).toBe(true);
    const baseLen = api.debugEnergySeriesLength();

    // Simulate a UI drag sequence that triggers computeForces calls (and thus USER_INTERACTION+energy)
    expect(typeof api.manipulation.beginDrag).toBe('function');
    expect(typeof api.manipulation.updateDrag).toBe('function');
    expect(typeof api.manipulation.endDrag).toBe('function');
    // Begin drag (provide a simple intersector)
    api.manipulation.beginDrag(() => ({ x: 0, y: 0, z: 0 }));
    // Update drag to new position (will ff.computeForces when idle)
    api.manipulation.updateDrag(() => ({ x: 1.2, y: 0, z: 0 }));
    // Allow first request/response cycle to complete
    await new Promise((r) => setTimeout(r, 60));
    // End drag (will ff.computeForces again when idle)
    api.manipulation.endDrag();
    // Wait for at least one more point to appear
    const ok2 = await (async () => {
      const t0 = Date.now();
      while (Date.now() - t0 < 400) {
        if (api.debugEnergySeriesLength() >= baseLen + 1) return true;
        await new Promise((r) => setTimeout(r, 20));
      }
      return false;
    })();
    expect(ok2).toBe(true);
    global.WebSocket = origWS;
  });
});
