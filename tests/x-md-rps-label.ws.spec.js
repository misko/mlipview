/** @jest-environment jsdom */

import { getWS } from '../public/fairchem_ws_client.js';

// Minimal Babylon scaffolding for viewer creation in jsdom.
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
    StandardMaterial: class {
      constructor() {
        this.diffuseColor = {};
        this.emissiveColor = {};
        this.specularColor = {};
      }
    },
    Color3: class {},
    Vector3: class {
      constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
      }
      static Up() {
        return new global.BABYLON.Vector3(0, 1, 0);
      }
    },
    Quaternion: class {
      static Identity() {
        return {};
      }
      static RotationAxis() {
        return {};
      }
    },
    Scene: class {},
    Matrix: class {
      static Compose() {
        return {};
      }
    },
  };
}

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: { add() {} },
    },
    camera: { attachControl: () => {} },
  }),
}));

describe('x-md-rps-label', () => {
  test('streamed MD frames update the RPS HUD label', async () => {
    window.__MLIPVIEW_TEST_MODE = true;
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    window.__MLIP_FEATURES = {
      RELAX_LOOP: false,
      MD_LOOP: true,
      ENERGY_TRACE: false,
      FORCE_VECTORS: false,
    };

    document.body.innerHTML = `
      <canvas id="viewer"></canvas>
      <div class="hud"><span id="rpsLabel">RPS: --</span></div>
      <div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>`;
    const energyCanvas = document.getElementById('energyCanvas');
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

    const OriginalWebSocket = global.WebSocket;
    class FakeSocket {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen?.();
        }, 0);
      }
      send() {}
      close() {}
      onopen() {}
      onmessage() {}
      onerror() {}
    }
    global.WebSocket = FakeSocket;

    const ws = getWS();
    ws.setTestHook(() => {});

    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.96, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });

    const origNow = performance.now;
    let now = 1_000;
    performance.now = () => now;

    try {
      api.startMDContinuous({ steps: 100, temperature: 1500 });
      await new Promise((r) => setTimeout(r, 0));

      for (let i = 0; i < 5; i++) {
        now += 100;
        ws.injectTestResult({
          positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
          forces: [],
          energy: -1.0,
          temperature: 1500,
        });
        // eslint-disable-next-line no-await-in-loop
        await new Promise((r) => setTimeout(r, 0));
      }
      const label = document.getElementById('rpsLabel');
      expect(label.textContent).toMatch(/^RPS:\s+\d+\.\d/);
      api.stopSimulation();
    } finally {
      performance.now = origNow;
      global.WebSocket = OriginalWebSocket;
    }
  });
});
