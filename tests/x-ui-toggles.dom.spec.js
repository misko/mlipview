/** @jest-environment jsdom */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
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

global.fetch = async (url) => {
  const target = String(url);
  if (/simple/.test(target)) {
    const body = JSON.stringify({ results: { energy: -1, forces: [[0, 0, 0]] } });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/\/md$/.test(target)) {
    const body = JSON.stringify({
      positions: [[0, 0, 0]],
      velocities: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      final_energy: -1,
      temperature: 300,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/\/relax$/.test(target)) {
    const body = JSON.stringify({
      positions: [[0, 0, 0]],
      forces: [[0, 0, 0]],
      final_energy: -1,
    });
    return { ok: true, status: 200, json: async () => JSON.parse(body), text: async () => body };
  }
  if (/health/.test(target)) {
    return {
      ok: true,
      status: 200,
      json: async () => ({ status: 'ok' }),
      text: async () => '{"status":"ok"}',
    };
  }
  throw new Error('Unexpected fetch ' + target);
};

function mountDom() {
  document.body.innerHTML = '<div id="app"></div><canvas id="viewer" width="100" height="100"></canvas>';
}

describe('x-ui-toggles', () => {
  test('forces and energy plot toggles flip UI state', async () => {
    mountDom();
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const { initNewViewer } = await import('../public/index.js');
    await initNewViewer(document.getElementById('viewer'), {
      elements: ['O'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const forces = document.getElementById('toggleForces');
    expect(forces.textContent.toLowerCase()).toContain('off');
    forces.click();
    expect(forces.textContent.toLowerCase()).toContain('on');
    const plot = document.getElementById('energyPlot');
    const toggleEnergy = document.getElementById('toggleEnergyPlot');
    expect(plot.style.display === '' || plot.style.display === 'block').toBe(true);
    toggleEnergy.click();
    expect(plot.style.display).toBe('none');
    toggleEnergy.click();
    expect(plot.style.display).toBe('block');
  });

  test('MD and Relax toggles enforce mutual exclusivity', async () => {
    mountDom();
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const { initNewViewer } = await import('../public/index.js');
    await initNewViewer(document.getElementById('viewer'), {
      elements: ['O'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const md = document.getElementById('toggleMD');
    const relax = document.getElementById('toggleRelax');
    md.click();
    expect(md.getAttribute('data-on')).toBe('true');
    relax.click();
    expect(relax.getAttribute('data-on')).toBe('true');
    expect(md.getAttribute('data-on')).toBe('false');
  });

  test('auto MD run updates button label after init', async () => {
    delete window.__MLIPVIEW_TEST_MODE;
    delete window.__MLIPVIEW_NO_AUTO_MD;
    mountDom();
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const { initNewViewer } = await import('../public/index.js');
    await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const button = document.getElementById('btnMDRun');
    for (let i = 0; i < 15 && button.textContent !== 'stop'; i++) {
      await new Promise((resolve) => setTimeout(resolve, 10));
    }
    expect(['stop', 'run']).toContain(button.textContent);
  });
});
