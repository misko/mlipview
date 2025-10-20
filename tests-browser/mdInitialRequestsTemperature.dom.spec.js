/** @jest-environment jsdom */

// Ensures the first 10 MD API temperatures sent after page load reflect the XYZ temperature if present.

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';
import { initTemperatureSlider } from '../public/ui/temperatureSlider.js';

function makeViewer() {
  const listeners = {};
  const bus = {
    on: (ev, fn) => {
      (listeners[ev] || (listeners[ev] = [])).push(fn);
    },
    emit: (ev) => {
      (listeners[ev] || []).forEach((fn) => fn());
    },
  };
  const state = {
    bus,
    elements: ['H', 'H', 'O'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: 0, y: 1, z: 0 },
    ],
    bonds: [],
    cell: {
      a: { x: 1, y: 0, z: 0 },
      b: { x: 0, y: 1, z: 0 },
      c: { x: 0, y: 0, z: 1 },
      enabled: false,
      originOffset: { x: 0, y: 0, z: 0 },
    },
    showCell: false,
    markCellChanged() {
      bus.emit('cellChanged');
    },
    markPositionsChanged() {
      bus.emit('positionsChanged');
    },
    markBondsChanged() {
      bus.emit('bondsChanged');
    },
    dynamics: {},
  };
  return { state, recomputeBonds: () => {} };
}

// Fixture with temp header
function xyzWithTemp(T = 500) {
  return `3\n temp=${T}K\nH 0 0 0\nH 1 0 0\nO 0 1 0\n`;
}

describe('MD initial requests use XYZ temperature when present', () => {
  beforeEach(() => {
    document.body.innerHTML = '<div id="hud"></div>';
    delete window.__MLIP_TARGET_TEMPERATURE;
    window.__MLIP_DEBUG_MD_TEMP = true;
  });

  test('first 10 MD calls use 500K', async () => {
    // Ensure viewer runs in test mode (bypass focus gating) and auto-MD is disabled
    window.__MLIPVIEW_TEST_MODE = true;
    window.__MLIPVIEW_NO_AUTO_MD = true;
    const xyz = xyzWithTemp(500);
    const parsed = parseXYZ(xyz);
    const viewer = makeViewer();
    // Build slider, then apply parsed so global is set and a sync event fires
    const hud = document.getElementById('hud');
    initTemperatureSlider({ hudEl: hud, getViewer: () => viewer });
    applyParsedToViewer(viewer, parsed);

    // Spy on fetch to intercept MD requests and record their temperatures
    const origFetch = global.fetch || window.fetch;
    const temps = [];
    global.fetch = async (url, opts) => {
      try {
        if (typeof opts?.body === 'string') {
          const body = JSON.parse(opts.body);
          if (url && /\/md\b/.test(String(url)) && typeof body.temperature !== 'undefined') {
            temps.push(body.temperature);
          }
        }
      } catch {}
      // Return a tiny valid MD response structure
      const response = {
        positions: viewer.state.positions.map((p) => [p.x, p.y, p.z]),
        forces: viewer.state.positions.map(() => [0, 0, 0]),
        final_energy: 0,
        velocities: viewer.state.positions.map(() => [0, 0, 0]),
        temperature: 500,
      };
      return {
        ok: true,
        status: 200,
        json: async () => response,
        text: async () => JSON.stringify(response),
      };
    };

    try {
      // Import index.js and run a short MD sequence to generate requests
      const mod = await import('../public/index.js');
      // Canvas stub with minimal DOM-like APIs required by scene.js
      const canvas = {
        width: 800,
        height: 600,
        getContext: () => ({}),
        addEventListener: () => {},
        removeEventListener: () => {},
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
      };
      // Initialize a viewer API instance
      const api = await mod.initNewViewer(canvas, {
        elements: viewer.state.elements,
        positions: viewer.state.positions,
        bonds: [],
      });
      // Inject our state (for this test we can adopt the viewer state's positions/elements)
      api.state.elements = viewer.state.elements;
      api.state.positions = viewer.state.positions;
      // Start MD with 10 steps; MD loop will read window.__MLIP_TARGET_TEMPERATURE each step
      await api.startMDContinuous({ steps: 10 });

      // Assert first 10 captured temperatures equal 500
      expect(temps.length).toBeGreaterThanOrEqual(10);
      for (let i = 0; i < 10; i++) {
        expect(temps[i]).toBe(500);
      }
    } finally {
      if (origFetch) global.fetch = origFetch;
    }
  });
});
