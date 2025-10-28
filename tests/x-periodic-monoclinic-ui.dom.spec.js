/** @jest-environment jsdom */

function wait(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

describe('x-periodic-monoclinic-ui', () => {
  async function setup() {
    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = '<div id="app"></div>';
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const listeners = {};
    const bus = {
      on: (ev, fn) => {
        (listeners[ev] || (listeners[ev] = [])).push(fn);
      },
      emit: (ev) => {
        (listeners[ev] || []).forEach((fn) => fn());
      },
    };
    const state = {
      bus,
      showCell: false,
      showGhostCells: false,
      cell: {
        a: { x: 10, y: 0, z: 0 },
        b: { x: 0, y: 12, z: 0 },
        c: { x: 0, y: 0, z: 15 },
        originOffset: { x: 0, y: 0, z: 0 },
        enabled: true,
      },
      toggleCellVisibilityEnhanced() {
        this.showCell = !this.showCell;
        this.bus.emit('cellChanged');
      },
      toggleGhostCells() {
        this.showGhostCells = !this.showGhostCells;
        this.bus.emit('cellChanged');
      },
      markCellChanged() {
        this.bus.emit('cellChanged');
      },
    };
    window.viewerApi = { state, setForceVectorsEnabled: () => {} };
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    await wait(0);
    return state;
  }

  test('Periodic tab renders controls', async () => {
    await setup();
    const section = document.getElementById('section-periodic');
    expect(section).toBeTruthy();
    const toggle = section.querySelector('#togglePBCHeader');
    const fieldA = section.querySelector('#cellA');
    const fieldB = section.querySelector('#cellB');
    const fieldC = section.querySelector('#cellC');
    const fieldBeta = section.querySelector('#cellBeta');
    expect(toggle && fieldA && fieldB && fieldC && fieldBeta).toBeTruthy();
  });

  test('enabling PBC unlocks nudgers and updates cell', async () => {
    const state = await setup();
    const section = document.getElementById('section-periodic');
    const toggle = section.querySelector('#togglePBCHeader');
    const buttons = ['#cellAPlus', '#cellBPlus', '#cellCPlus', '#cellBetaPlus'].map((sel) =>
      section.querySelector(sel)
    );
    expect(buttons.every((btn) => btn.disabled)).toBe(true);
    toggle.click();
    expect(state.showCell).toBe(true);
    expect(buttons.some((btn) => !btn.disabled)).toBe(true);
    buttons.forEach((btn) => btn.click());
    const cVec = state.cell.c;
    const cLen = Math.hypot(cVec.x, cVec.y, cVec.z);
    expect(cLen).toBeGreaterThan(0);
  });

  test('beta hold clamps at < 180°', async () => {
    const state = await setup();
    jest.useFakeTimers();
    const section = document.getElementById('section-periodic');
    const toggle = section.querySelector('#togglePBCHeader');
    toggle.click();
    const beta = section.querySelector('#cellBeta');
    const betaPlus = section.querySelector('#cellBetaPlus');
    beta.textContent = '178';
    betaPlus.dispatchEvent(new Event('pointerdown', { bubbles: true }));
    jest.advanceTimersByTime(400);
    window.dispatchEvent(new Event('pointerup', { bubbles: true }));
    jest.runOnlyPendingTimers();
    const value = parseFloat(beta.textContent);
    expect(value).toBeLessThanOrEqual(179);
    expect(state.cell.enabled).toBe(true);
    jest.useRealTimers();
  });
});
