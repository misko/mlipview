/** @jest-environment jsdom */
/**
 * Mock browser test: ?autoMD=0 disables auto-start MD and baseline energy is plotted from WS
 */
import { jest } from '@jest/globals';
import { stubWebSocketAndHook } from './utils/wsTestStub.js';

describe('autoMD=0 baseline energy', () => {
  beforeEach(() => { if (!global.window) global.window = global; });

  test('loads with autoMD=0 and plots baseline energy from WS', async () => {
    // Seed jsdom env and WS stub
    window.__MLIPVIEW_TEST_MODE = true; // disable real render loop
    window.location = { protocol: 'http:', host: '127.0.0.1:4000', search: '?autoMD=0' };
    window.__MLIPVIEW_SERVER = 'ws://127.0.0.1:8000';
    const { emit } = stubWebSocketAndHook();
    // Auto-respond with an energy frame when baseline USER_INTERACTION is sent
    const __origHook = window.__WS_TEST_HOOK__;
    window.__WS_TEST_HOOK__ = (msg)=>{
      try { __origHook && __origHook(msg); } catch {}
      try {
        if (msg && (msg.type === 'USER_INTERACTION' || msg.type === 1)) {
          // simulate immediate simple_calculate energy
          setTimeout(()=> emit({ seq: 1, energy: -1.234, forces: [[0,0,0],[0,0,0]] }), 0);
        }
      } catch {}
    };

  // Minimal DOM for energy canvas + viewer canvas with addEventListener
  document.body.innerHTML = '<canvas id="energyCanvas" width="200" height="50"></canvas><div id="energyLabel"></div>';
  const viewer = document.createElement('canvas'); viewer.id='viewer'; viewer.addEventListener=()=>{}; document.body.appendChild(viewer);

  const { initNewViewer } = await import('../public/index.js');
  const api = await initNewViewer(viewer, { elements:['H','H'], positions:[{x:0,y:0,z:0},{x:1, y:0, z:0}], bonds:[] });

    // Give the event loop a tick for handlers
    await new Promise(r=> setTimeout(r, 5));

    // Expect energy series to have at least one point
    expect(typeof api.debugEnergySeriesLength).toBe('function');
    expect(api.debugEnergySeriesLength()).toBeGreaterThanOrEqual(1);
  });
});
