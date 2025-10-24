import { test, expect } from './fixtures.js';

const CC_MIN = 0.9;
const CC_MAX = 2.0;
const CH_MIN = 0.7;
const CH_MAX = 1.5;

test('x-benzene bond rotation keeps lengths sane', async ({ page, loadViewerPage }) => {
  await loadViewerPage({
    query: { mol: 'molecules/benzene.xyz', autoMD: 0 },
    disableAutoMd: true,
  });

  const result = await page.evaluate(async () => {
    const api = window.viewerApi;
    const { state, manipulation, view } = api;

    if (typeof state.toggleCellVisibilityEnhanced === 'function' && !state.showCell) {
      state.toggleCellVisibilityEnhanced();
    }
    if (typeof state.toggleGhostCells === 'function' && !state.showGhostCells) {
      state.toggleGhostCells();
    }
    state.markCellChanged?.();

    const bonds = Array.isArray(state.bonds) ? state.bonds : [];
    const ccBond = bonds.find((b) => state.elements[b.i] === 'C' && state.elements[b.j] === 'C');
    if (!ccBond) return { cc: [], ch: [], ghostBondCount: 0 };

    const bondIdx = bonds.indexOf(ccBond);
    api.debugSelectBond({
      i: ccBond.i,
      j: ccBond.j,
      index: bondIdx >= 0 ? bondIdx : 0,
      key: ccBond.key ?? `${ccBond.i}-${ccBond.j}`,
      orientation: 0,
    });

    const angles = [Math.PI / 18, Math.PI / 18, -Math.PI / 12];
    for (const angle of angles) {
      manipulation.rotateBond(angle);
    }
    view.rebuildGhosts?.();

    await new Promise((resolve) => setTimeout(resolve, 20));

    const snapshot = api.debugGhostSnapshot();
    const metrics = api.debugBondMetrics();

    const cc = metrics
      .filter((m) => m.elements[0] === 'C' && m.elements[1] === 'C')
      .map((m) => m.length);
    const ch = metrics
      .filter((m) => m.elements.includes('C') && m.elements.includes('H'))
      .map((m) => m.length);

    return { cc, ch, ghostBondCount: snapshot.ghostBondCount };
  });

  expect(result.cc.length).toBeGreaterThan(0);
  expect(result.ch.length).toBeGreaterThan(0);
  expect(result.ghostBondCount).toBeGreaterThan(0);

  for (const len of result.cc) {
    expect(len).toBeGreaterThanOrEqual(CC_MIN);
    expect(len).toBeLessThanOrEqual(CC_MAX);
  }
  for (const len of result.ch) {
    expect(len).toBeGreaterThanOrEqual(CH_MIN);
    expect(len).toBeLessThanOrEqual(CH_MAX);
  }
});
