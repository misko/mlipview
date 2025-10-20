/** @jest-environment jsdom */
// Test: Temperature slider and MD parameter stay in sync on initial load and after updates (WS-only).

beforeAll(() => {
  // Minimal BABYLON stubs
  global.BABYLON = global.BABYLON || {
    Engine: function () {
      this.runRenderLoop = () => {};
    },
    Scene: function () {
      this.onPointerObservable = {};
      this.render = () => {};
    },
    Color3: function () {},
    MeshBuilder: { CreateSphere: () => ({}) },
    StandardMaterial: function () {},
    ArcRotateCamera: function () {
      this.attachControl = () => {};
    },
  };
});

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => {} },
    scene: { render: () => {}, onPointerObservable: { add: () => {} } },
    camera: {},
  }),
}));
import { getWS } from '../public/fairchem_ws_client.js';

async function setup() {
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"></div>`;
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
  const sent = [];
  ws.setTestHook((m) => sent.push(m));
  const { initNewViewer } = await import('../public/index.js');
  const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
  const canvas = document.getElementById('viewer');
  const hud = document.querySelector('.hud');
  const api = await initNewViewer(canvas, {
    elements: ['O', 'H', 'H'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: 0, y: 1, z: 0 },
    ],
    bonds: [],
  });
  window.viewerApi = api;
  initTemperatureSlider({ hudEl: hud, getViewer: () => api });
  return {
    api,
    ws,
    sent,
    restore: () => {
      global.WebSocket = origWS;
    },
  };
}

describe('Temperature slider sync', () => {
  test('initial load uses 1500K and mdStep payload matches', async () => {
    const { api, ws, sent, restore } = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    expect(slider).toBeTruthy();
    expect(label.textContent).toContain('1500');
    // Trigger an MD step; verify outgoing START_SIMULATION temperature
    await Promise.resolve().then(() => {}); // tick for ws connect
    const p = api.mdStep({ temperature: window.__MLIP_TARGET_TEMPERATURE });
    // Simulate a response frame to resolve mdStep
    ws.injectTestResult({
      positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
      forces: [],
      energy: -1.0,
      temperature: window.__MLIP_TARGET_TEMPERATURE,
    });
    await p;
    const lastStart = sent.reverse().find((m) => m && m.type != null);
    expect(sent.length).toBeGreaterThan(0);
    // Find a START_SIMULATION with simulationParams.temperature
    const startMsg = sent.find(
      (m) => m.simulationParams && typeof m.simulationParams.temperature === 'number'
    );
    expect(startMsg).toBeTruthy();
    expect(startMsg.simulationParams.temperature).toBe(window.__MLIP_TARGET_TEMPERATURE);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(1500);
    restore();
  });

  test('moving slider updates __MLIP_TARGET_TEMPERATURE and mdStep payload follows', async () => {
    const { api, ws, sent, restore } = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    // Move slider to max (near 3000)
    slider.value = slider.max;
    slider.dispatchEvent(new Event('input'));
    const target = window.__MLIP_TARGET_TEMPERATURE;
    expect(label.textContent).toContain(String(target));
    const p = api.mdStep({ temperature: target });
    ws.injectTestResult({
      positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
      forces: [],
      energy: -1.0,
      temperature: target,
    });
    await p;
    const startMsg = sent.find(
      (m) => m.simulationParams && typeof m.simulationParams.temperature === 'number'
    );
    expect(startMsg).toBeTruthy();
    expect(startMsg.simulationParams.temperature).toBe(target);
    restore();
  });
});
