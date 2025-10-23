/** @jest-environment jsdom */
// Test: friction slider drives friction field in /serve/md POST bodies with default 0.5.

// Minimal BABYLON stubs (align with other jsdom tests)
beforeAll(() => {
  if (!global.BABYLON) {
    global.BABYLON = {
      TransformNode: class { },
      MeshBuilder: {
        CreateCylinder: () => ({
          dispose() { },
          setEnabled() { },
          position: { set() { } },
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
      Color3: class { },
      Vector3: class {
        constructor(x = 0, y = 0, z = 0) {
          this.x = x;
          this.y = y;
          this.z = y;
        }
        static Up() {
          return new global.BABYLON.Vector3(0, 1, 0);
        }
        static Dot(a, b) {
          return a.x * b.x + a.y * b.y + a.z * b.z;
        }
        static Cross(a, b) {
          return new global.BABYLON.Vector3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
          );
        }
        normalize() {
          const m = Math.hypot(this.x, this.y, this.z) || 1;
          this.x /= m;
          this.y /= m;
          this.z /= m;
          return this;
        }
        rotateByQuaternionToRef(_, out) {
          out.x = this.x;
          out.y = this.y;
          out.z = this.z;
          return out;
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
      Scene: class { },
    };
  }
});

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => { } },
    scene: {
      meshes: [],
      render: () => { },
      onPointerObservable: {
        _l: [],
        add(fn) {
          this._l.push(fn);
        },
      },
    },
    camera: { attachControl: () => { } },
  }),
}));

// Capture POST bodies for /serve/md
const mdBodies = [];

global.fetch = async function (url, opts = {}) {
  if (typeof url === 'string' && /\/serve\/md$/.test(url)) {
    try {
      mdBodies.push(JSON.parse(opts.body || '{}'));
    } catch {
      mdBodies.push({ parseError: true, raw: opts.body });
    }
    const last = mdBodies[mdBodies.length - 1] || {};
    // Return a plausible MD response echoing temperature
    const body = JSON.stringify({
      positions: [
        [0, 0, 0],
        [0.96, 0, 0],
        [-0.24, 0.93, 0],
      ],
      velocities: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ],
      forces: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ],
      final_energy: -1.0,
      steps_completed: 1,
      temperature: last.temperature || 0,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (typeof url === 'string' && /\/serve\/simple$/.test(url)) {
    const body = JSON.stringify({ energy: -1.2, forces: [[0, 0, 0]], positions: [[0, 0, 0]] });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (typeof url === 'string' && /\/serve\/relax$/.test(url)) {
    const body = JSON.stringify({
      positions: [
        [0, 0, 0],
        [0.96, 0, 0],
        [-0.24, 0.93, 0],
      ],
      forces: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ],
      final_energy: -1.1,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (typeof url === 'string' && /\/serve\/health$/.test(url)) {
    const body = JSON.stringify({ ok: true });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  throw new Error('Unexpected fetch ' + url);
};

async function setupViewer() {
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  window.__MLIP_FEATURES = {
    RELAX_LOOP: false,
    MD_LOOP: true,
    ENERGY_TRACE: false,
    FORCE_VECTORS: false,
  };
  const canvas = document.createElement('canvas');
  canvas.id = 'viewer';
  document.body.appendChild(canvas);
  const hud = document.createElement('div');
  hud.className = 'hud';
  document.body.appendChild(hud);
  hud.innerHTML = `<span></span>
    <button id="btnRelax"></button>
    <button id="btnRelaxRun"></button>
    <span></span>
    <button id="btnMD"></button>
    <button id="btnMDRun"></button>
    <select id="forceProviderSel"></select>
    <button id="btnCell"></button>
    <button id="btnGhosts"></button>
    <button id="btnToggleForces"></button>
    <span id="status">Ready</span>
    <select id="xrModeSelect"></select>`;
  const energyWrapper = document.createElement('div');
  energyWrapper.id = 'energyPlot';
  document.body.appendChild(energyWrapper);
  const energyCanvas = document.createElement('canvas');
  energyCanvas.id = 'energyCanvas';
  energyCanvas.width = 260;
  energyCanvas.height = 80;
  energyCanvas.getContext = () => ({
    clearRect() { },
    beginPath() { },
    moveTo() { },
    lineTo() { },
    stroke() { },
    arc() { },
    fill() { },
    fillRect() { },
    strokeStyle: null,
    lineWidth: 1,
    fillStyle: null,
  });
  energyWrapper.appendChild(energyCanvas);
  const energyLabel = document.createElement('div');
  energyLabel.id = 'energyLabel';
  energyWrapper.appendChild(energyLabel);
  const mod = await import('../public/index.js');
  const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
  const { initFrictionSlider } = await import('../public/ui/frictionSlider.js');
  window.__MLIP_DEV_MODE = true;
  const viewer = await mod.initNewViewer(canvas, {
    elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 0.96, y: 0, z: 0 },
      { x: -0.24, y: 0.93, z: 0 },
    ],
    bonds: [],
  });
  window.viewerApi = viewer;
  initTemperatureSlider({ hudEl: hud, getViewer: () => viewer });
  initFrictionSlider({ hudEl: hud });
  return viewer;
}

describe('MD friction slider', () => {
  test('default 0.5 friction is sent; slider updates payload', async () => {
    const viewer = await setupViewer();
    const slider = document.getElementById('mdFrictionSlider');
    expect(slider).toBeTruthy();
    mdBodies.length = 0;
    await viewer.mdStep({ temperature: 298 });
    expect(mdBodies.length).toBe(1);
    expect(mdBodies[0].friction).toBeCloseTo(0.5, 5);
    // Change friction to 1.25 and verify
    slider.value = '1.25';
    slider.dispatchEvent(new Event('input'));
    await viewer.mdStep({ temperature: 298 });
    expect(mdBodies.length).toBe(2);
    expect(mdBodies[1].friction).toBeCloseTo(1.25, 5);
  });
});
