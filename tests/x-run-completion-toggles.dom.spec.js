/** @jest-environment jsdom */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: {
        _listeners: [],
        add(fn) {
          this._listeners.push(fn);
        },
      },
    },
    camera: { attachControl: () => {} },
  }),
}));

global.fetch = async (url) => {
  if (/\/serve\/relax$/.test(String(url))) {
    return {
      ok: true,
      status: 200,
      json: async () => ({ positions: [[0, 0, 0]], results: {} }),
    };
  }
  if (/\/serve\/md$/.test(String(url))) {
    return {
      ok: true,
      status: 200,
      json: async () => ({
        positions: [[0, 0, 0]],
        final_energy: -1,
        velocities: [[0, 0, 0]],
        temperature: 300,
      }),
    };
  }
  return { ok: true, status: 200, json: async () => ({ results: {} }) };
};

function wait(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function setup() {
  window.__MLIPVIEW_TEST_MODE = true;
  window.__MLIPVIEW_PANEL_SHORT_RUNS = true;
  window.__MLIPVIEW_NO_AUTO_MD = true;
  document.body.innerHTML = '<canvas id="viewer"></canvas><div id="app"></div>';
  const { initNewViewer } = await import('../public/index.js');
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const viewer = await initNewViewer(document.getElementById('viewer'), {
    elements: [{ Z: 8 }],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  window.viewerApi = viewer;
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  return viewer;
}

async function nextUiTick() {
  await wait(350);
}

describe('x-run-completion-toggles', () => {
  test('MD toggle resets to Off and button text to run', async () => {
    await setup();
    const mdToggle = document.getElementById('toggleMD');
    const mdButton = document.getElementById('btnMDRun');
    mdToggle.click();
    await nextUiTick();
    await wait(900);
    await nextUiTick();
    expect(mdToggle.getAttribute('data-on')).toBe('false');
    expect(mdButton.textContent).toBe('run');
  });

  test('Relax toggle resets to Off and button text to run', async () => {
    await setup();
    const relaxToggle = document.getElementById('toggleRelax');
    const relaxButton = document.getElementById('btnRelaxRun');
    relaxToggle.click();
    await nextUiTick();
    await wait(900);
    await nextUiTick();
    expect(relaxToggle.getAttribute('data-on')).toBe('false');
    expect(relaxButton.textContent).toBe('run');
  });
});
