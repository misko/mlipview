/** @jest-environment jsdom */

// Desktop UI bond selection regression test: ensure picking a bond selects it,
// shows highlight, updates Selection panel details, and rotate +/- buttons work.

function wait(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

describe('desktop bond selection via picking', () => {
  async function setup() {
    // Minimal DOM
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    // Enable test mode to relax rendering and forces visualization defaults
    window.__MLIPVIEW_TEST_MODE = true;
    // Disable touch controls to avoid synthetic pointer events interfering with desktop pick
    window.__MLIPVIEW_NO_TOUCH = true;
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    // Expose for desktop panel helpers
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    // Allow initial listeners to attach and selection UI to render
    await wait(0);
    return viewer;
  }

  test('picking a bond highlights and enables rotate in Selection panel', async () => {
    const v = await setup();
    const canvas = document.getElementById('viewer');
    // Prepare a bond pick result: return the bond master mesh and index 0
    const bondGroups = v.view._internals.bondGroups;
    const firstGroup = Array.from(bondGroups.values())[0];
    expect(firstGroup).toBeTruthy();
    // Ensure no stray multiPick path interferes
    v.scene.multiPick = undefined;
    v.scene.pick = () => ({ hit: true, pickedMesh: firstGroup.master, thinInstanceIndex: 0 });
    // Dispatch a pointerdown to trigger picking service
    await wait(0);
    const evt =
      typeof PointerEvent !== 'undefined'
        ? new PointerEvent('pointerdown', { clientX: 10, clientY: 10 })
        : new Event('pointerdown');
    canvas.dispatchEvent(evt);
    // Selection state should be bond
    expect(v.selection.get().kind).toBe('bond');
    // Highlight bond mesh should be visible
    const hl = v.view._internals.highlight;
    expect(hl && hl.bond && hl.bond.isVisible).toBe(true);
    // Selection panel elements
    const name = document.getElementById('selElementName');
    const pos = document.getElementById('selPosition');
    const bondLen = document.getElementById('bondLength');
    expect((name.textContent || '').toLowerCase()).toContain('carbon');
    expect((name.textContent || '').toLowerCase()).toContain('hydrogen');
    expect(pos.textContent).toContain(') – (');
    expect(bondLen.textContent).toContain('Å');
    // Rotate buttons should be visible; clicking should increment positions version
    const rotateBtns = document.getElementById('rotateBtns');
    const minus = document.getElementById('bondRotMinus');
    const plus = document.getElementById('bondRotPlus');
    expect(rotateBtns && rotateBtns.style.display).toBe('inline-flex');
    const before = v.state.versions.positions;
    plus.click();
    // Allow any async microtasks
    await wait(0);
    expect(v.state.versions.positions).toBeGreaterThan(before);
    // Clicking minus should also rotate
    const before2 = v.state.versions.positions;
    minus.click();
    await wait(0);
    expect(v.state.versions.positions).toBeGreaterThan(before2);
  });
});
