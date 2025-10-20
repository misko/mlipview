/** @jest-environment jsdom */

// Verify that clicking MD/Relax starts full runs (not limited to ~5 steps)
// and that programmatic completion flips UI via polling.

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

// Mock server endpoints: respond with progressing data for multiple steps
let mdStepCount = 0;
let relaxStepCount = 0;
global.fetch = async function (url, opts) {
  const u = String(url);
  if (/md/.test(u)) {
    mdStepCount++;
    return {
      ok: true,
      status: 200,
      json: async () => ({
        positions: [[0, 0, 0]],
        final_energy: -1.0 - mdStepCount,
        velocities: [[0, 0, 0]],
        temperature: 300,
      }),
    };
  }
  if (/relax/.test(u)) {
    relaxStepCount++;
    return {
      ok: true,
      status: 200,
      json: async () => ({
        positions: [[0, 0, 0]],
        results: { energy: -2.0 - relaxStepCount, forces: [[0, 0, 0]] },
      }),
    };
  }
  if (/simple/.test(u)) {
    return {
      ok: true,
      status: 200,
      json: async () => ({ results: { energy: -3.0, forces: [[0, 0, 0]] } }),
    };
  }
  return { ok: true, status: 200, json: async () => ({}) };
};

function wait(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

async function setup() {
  window.__MLIPVIEW_TEST_MODE = true;
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

async function nextUiTick() {
  await wait(300);
}

describe('Panel starts full runs', () => {
  test('MD starts continuous run (more than 5 steps)', async () => {
    mdStepCount = 0;
    await setup();
    const md = document.getElementById('toggleMD');
    md.click(); // start
    await wait(900); // let several steps accumulate
    expect(mdStepCount).toBeGreaterThan(5);
  });

  test('Relax starts continuous run (more than 5 steps)', async () => {
    relaxStepCount = 0;
    await setup();
    // ensure MD not running
    try {
      window.viewerApi.stopSimulation();
    } catch {}
    const rx = document.getElementById('toggleRelax');
    rx.click();
    await wait(900);
    expect(relaxStepCount).toBeGreaterThan(5);
  });
});
