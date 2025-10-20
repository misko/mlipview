/** @jest-environment jsdom */

// Verify MD/Relax toggles and legacy buttons flip back to Off/'run' when runs complete.

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

// Minimal fetch stub: immediately returns empty ok json to drive baseline energy
global.fetch = async function (url, opts) {
  // Return minimal payloads to drive MD/Relax steps
  if (/relax/.test(String(url))) {
    return { ok: true, status: 200, json: async () => ({ positions: [[0, 0, 0]], results: {} }) };
  }
  if (/md/.test(String(url))) {
    return {
      ok: true,
      status: 200,
      json: async () => ({
        positions: [[0, 0, 0]],
        final_energy: -1.0,
        velocities: [[0, 0, 0]],
        temperature: 300,
      }),
    };
  }
  return { ok: true, status: 200, json: async () => ({ results: {} }) };
};

function wait(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

async function setup() {
  // Disable auto-MD for deterministic tests
  window.__MLIPVIEW_TEST_MODE = true;
  window.__MLIPVIEW_PANEL_SHORT_RUNS = true;
  window.__MLIPVIEW_NO_AUTO_MD = true;
  document.body.innerHTML = `
    <canvas id="viewer"></canvas>
    <div id="app"></div>
  `;
  const { initNewViewer } = await import('../public/index.js');
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const canvas = document.getElementById('viewer');
  const viewer = await initNewViewer(canvas, {
    elements: [{ Z: 8 }],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  window.viewerApi = viewer;
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  return viewer;
}

// Helper: wait until UI sync interval has a chance to update
async function nextUiTick() {
  await wait(350);
}

describe('Run completion flips toggles Off', () => {
  test('MD toggle returns to Off and legacy btn to run', async () => {
    const v = await setup();
    // Force very short MD run: 2 steps
    setTimeout(() => {
      window.__MLIPVIEW_NO_AUTO_MD = true;
    }, 0);
    const md = document.getElementById('toggleMD');
    expect(md).toBeTruthy();
    // Click to start MD (panel wiring awaits completion and then flips Off)
    md.click();
    await nextUiTick(); // start
    await wait(900); // allow mocked run to finish and UI sync to reflect
    await nextUiTick(); // sync back
    expect(md.getAttribute('data-on')).toBe('false');
    const legacy = document.getElementById('btnMDRun');
    expect(legacy && legacy.textContent).toBe('run');
  });

  test('Relax toggle returns to Off and legacy btn to run', async () => {
    const v = await setup();
    const relax = document.getElementById('toggleRelax');
    expect(relax).toBeTruthy();
    relax.click();
    await nextUiTick();
    await wait(900);
    await nextUiTick();
    expect(relax.getAttribute('data-on')).toBe('false');
    const legacy = document.getElementById('btnRelaxRun');
    expect(legacy && legacy.textContent).toBe('run');
  });
});
