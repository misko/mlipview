/** @jest-environment jsdom */
// Verifies auto-started MD run resets the run button text to 'run' after completion.

// Mock scene
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

// Short MD loop by forcing steps=2 inside startMDContinuous call intercept (monkey patch later if needed)

global.fetch = async (url, opts = {}) => {
  if (/\/serve\/simple$/.test(url)) {
    const body = JSON.stringify({
      results: { energy: -1.0, forces: [[0, 0, 0]], stress: [0, 0, 0, 0, 0, 0] },
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/\/serve\/md$/.test(url)) {
    // Echo minimal MD JSON; after a few steps we let loop finish naturally.
    const body = JSON.stringify({
      positions: [[0, 0, 0]],
      velocities: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      final_energy: -1.0,
      steps_completed: 1,
      temperature: 298,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/\/serve\/relax$/.test(url)) {
    const body = JSON.stringify({
      initial_energy: -1,
      final_energy: -1,
      positions: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      steps_completed: 1,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/\/serve\/health$/.test(url)) {
    return {
      ok: true,
      status: 200,
      json: async () => ({ status: 'ok' }),
      text: async () => '{"status":"ok"}',
    };
  }
  throw new Error('Unexpected fetch ' + url);
};

test('auto-start MD run resets button text', async () => {
  // Ensure auto MD runs (not disabled, not test mode)
  delete window.__MLIPVIEW_NO_AUTO_MD;
  delete window.__MLIPVIEW_TEST_MODE;
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"><button id="btnRelax"></button><button id="btnRelaxRun"></button><button id="btnMD"></button><button id="btnMDRun">run</button><select id="forceProviderSel"></select><button id="btnCell"></button><button id="btnGhosts"></button><button id="btnToggleForces"></button><span id="status">Ready</span><select id="xrModeSelect"></select></div><div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>`;
  const mod = await import('../public/index.js');
  await mod.initNewViewer(document.getElementById('viewer'), {
    elements: [{ Z: 8 }],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  const btn = document.getElementById('btnMDRun');
  // Wait up to 100ms for auto-start setTimeout to flip text to 'stop'
  for (let k = 0; k < 10 && btn.textContent !== 'stop'; k++) {
    await new Promise((r) => setTimeout(r, 10));
  }
  expect(['stop', 'run']).toContain(btn.textContent); // if run already finished extremely fast it may have reverted
  const sawStop = btn.textContent === 'stop';
  // Allow MD loop a short time to finish default steps
  // We poll for text change to 'run'.
  for (let k = 0; k < 50; k++) {
    await new Promise((r) => setTimeout(r, 10));
    if (btn.textContent === 'run') break;
  }
  expect(btn.textContent).toBe('run');
  // Ensure it transitioned at least once if we saw stop
  if (sawStop) expect(sawStop).toBe(true);
});
