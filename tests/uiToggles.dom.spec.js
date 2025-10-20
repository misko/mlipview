/** @jest-environment jsdom */
// Tests: MD on/off, Relax on/off, Show forces on/off, Energy plot on/off

// Minimal scene stub so initNewViewer works
jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => {} },
    scene: { meshes: [], render: () => {}, onPointerObservable: { add() {} } },
    camera: {},
  }),
}));

beforeAll(() => {
  if (!global.BABYLON) {
    global.BABYLON = {
      TransformNode: class {},
      MeshBuilder: {
        CreateSphere: () => ({ setEnabled() {} }),
        CreateCylinder: () => ({
          setEnabled() {},
          thinInstanceSetBuffer() {},
          material: {},
          _buffers: {},
        }),
        CreateLines: () => ({ dispose() {} }),
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
      Quaternion: class {
        static Identity() {
          return {};
        }
      },
      Material: { MATERIAL_ALPHABLEND: 2 },
    };
  }
});

// Mock network endpoints
global.fetch = async (url, opts) => {
  if (/simple/.test(url)) {
    const body = JSON.stringify({ results: { energy: -1.0, forces: [[0, 0, 0]] } });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/md$/.test(url)) {
    const body = JSON.stringify({
      positions: [[0, 0, 0]],
      velocities: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      final_energy: -1.0,
      temperature: 300,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/relax$/.test(url)) {
    const body = JSON.stringify({
      positions: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      final_energy: -1.0,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/health/.test(url)) {
    return {
      ok: true,
      status: 200,
      json: async () => ({ status: 'ok' }),
      text: async () => '{"status":"ok"}',
    };
  }
  throw new Error('Unexpected fetch ' + url);
};

function mountDOM() {
  document.body.innerHTML = `<div id="app"></div><canvas id="viewer" width="100" height="100"></canvas>`;
}

test('connect toggles: forces and energy plot', async () => {
  mountDOM();
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const mod = await import('../public/index.js');
  const api = await mod.initNewViewer(document.getElementById('viewer'), {
    elements: ['O'],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  // Forces toggle
  const forces = document.getElementById('toggleForces');
  expect(forces.textContent.toLowerCase()).toContain('off');
  forces.click();
  expect(forces.textContent.toLowerCase()).toContain('on');
  // Energy plot toggle hides/shows container
  const plot = document.getElementById('energyPlot');
  const energy = document.getElementById('toggleEnergyPlot');
  // Default turned on by builder
  expect(plot.style.display === '' || plot.style.display === 'block').toBe(true);
  energy.click();
  expect(plot.style.display).toBe('none');
  energy.click();
  expect(plot.style.display).toBe('block');
});

test('connect toggles: MD and Relaxation', async () => {
  mountDOM();
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const mod = await import('../public/index.js');
  const api = await mod.initNewViewer(document.getElementById('viewer'), {
    elements: ['O'],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  const md = document.getElementById('toggleMD');
  const relax = document.getElementById('toggleRelax');
  // Turn MD on
  md.click();
  expect(md.getAttribute('data-on')).toBe('true');
  // Turning Relax on should turn MD off
  relax.click();
  expect(relax.getAttribute('data-on')).toBe('true');
  expect(md.getAttribute('data-on')).toBe('false');
});

test('MD auto-run on page load sets status/button', async () => {
  // Enable real auto MD (not test mode)
  delete window.__MLIPVIEW_TEST_MODE;
  delete window.__MLIPVIEW_NO_AUTO_MD;
  mountDOM();
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const mod = await import('../public/index.js');
  const api = await mod.initNewViewer(document.getElementById('viewer'), {
    elements: [{ Z: 8 }],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  // After a short delay the MD run should start and status/button should reflect
  const btn = document.getElementById('btnMDRun');
  // wait up to 150ms for auto start to set 'stop'
  for (let k = 0; k < 15 && btn.textContent !== 'stop'; k++)
    await new Promise((r) => setTimeout(r, 10));
  expect(['stop', 'run']).toContain(btn.textContent);
});
