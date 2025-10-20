/** @jest-environment jsdom */
// Verifies that the instantaneous temperature from WS updates the HUD element #instTemp.

import { getWS } from '../public/fairchem_ws_client.js';

beforeAll(() => {
  if (!global.BABYLON) {
    global.BABYLON = {
      TransformNode: class {},
      MeshBuilder: {
        CreateCylinder: () => ({
          dispose() {},
          setEnabled() {},
          position: { set() {} },
          rotationQuaternion: null,
          scaling: {},
          isPickable: false,
          visibility: 1,
        }),
      },
      StandardMaterial: class {},
      Color3: class {},
      Vector3: class {
        constructor(x = 0, y = 0, z = 0) {
          this.x = x;
          this.y = y;
          this.z = z;
        }
      },
      Quaternion: class {},
      Scene: class {},
    };
  }
});

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: {
        _l: [],
        add(fn) {
          this._l.push(fn);
        },
      },
    },
    camera: { attachControl: () => {} },
  }),
}));

async function setup() {
  // Minimal DOM similar to app
  document.body.innerHTML = '';
  const canvas = document.createElement('canvas');
  canvas.id = 'viewer';
  canvas.addEventListener = () => {};
  document.body.appendChild(canvas);
  const hud = document.createElement('div');
  hud.className = 'hud';
  document.body.appendChild(hud);
  const inst = document.createElement('span');
  inst.id = 'instTemp';
  inst.textContent = 'T: --.- K';
  hud.appendChild(inst);
  const energyWrapper = document.createElement('div');
  energyWrapper.id = 'energyPlot';
  document.body.appendChild(energyWrapper);
  const energyCanvas = document.createElement('canvas');
  energyCanvas.id = 'energyCanvas';
  energyCanvas.width = 260;
  energyCanvas.height = 80;
  energyCanvas.getContext = () => ({
    clearRect() {},
    beginPath() {},
    moveTo() {},
    lineTo() {},
    stroke() {},
    arc() {},
    fill() {},
    fillRect() {},
    strokeStyle: null,
    lineWidth: 1,
    fillStyle: null,
  });
  energyWrapper.appendChild(energyCanvas);
  const energyLabel = document.createElement('div');
  energyLabel.id = 'energyLabel';
  energyWrapper.appendChild(energyLabel);
  const mod = await import('../public/index.js');
  const viewer = await mod.initNewViewer(canvas, {
    elements: [{ Z: 8 }],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  return viewer;
}

describe('Instantaneous MD temperature HUD', () => {
  test('updates #instTemp after mdStep', async () => {
    window.__MLIPVIEW_TEST_MODE = true;
    // Stub WS and capture outgoing messages
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
    const viewer = await setup();
    const el = document.getElementById('instTemp');
    expect(el).toBeTruthy();
    expect(el.textContent).toMatch(/--/);
    // Kick off mdStep; it will subscribe for a single frame via WS
    const p = viewer.mdStep({ temperature: 500 });
    // Allow listener to attach
    await Promise.resolve();
    await new Promise((r) => setTimeout(r, 0));
    // Emit a frame with instantaneous temperature (500 + 12.34)
    const instT = 512.34;
    ws.injectTestResult({
      positions: viewer.state.positions.map((p) => [p.x, p.y, p.z]),
      forces: [[0, 0, 0]],
      velocities: [[0, 0, 0]],
      energy: -1.0,
      temperature: instT,
    });
    await p; // resolve mdStep
    // Outgoing START_SIMULATION captured with temperature (enum is numeric)
    const startMsg = sent.find(
      (m) => m && m.simulationParams && typeof m.simulationParams.temperature === 'number'
    );
    expect(startMsg && startMsg.simulationParams && startMsg.simulationParams.temperature).toBe(
      500
    );
    // HUD should reflect instantaneous value (rounded to 1 decimal)
    expect(el.textContent).toMatch(/T: 512\.3/);
    global.WebSocket = origWS;
  });
});
