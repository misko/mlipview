import { test, expect } from './fixtures.js';

test('x-bond highlight toggles between atom and bond selection', async ({
  page,
  loadViewerPage,
}) => {
  await loadViewerPage({
    query: { mol: 'molecules/water.xyz', autoMD: 0 },
    disableAutoMd: true,
  });

  const states = await page.evaluate(() => {
    const api = window.viewerApi;
    const bonds = Array.isArray(api.state.bonds) ? api.state.bonds : [];
    if (bonds.length === 0) return null;

    api.debugSelectAtom(0);
    const afterAtom = api.debugHighlightState();

    const first = bonds[0];
    api.debugSelectBond({
      i: first.i,
      j: first.j,
      key: first.key ?? `${first.i}-${first.j}`,
      index: first.index ?? 0,
    });
    const afterBond = api.debugHighlightState();
    return { afterAtom, afterBond };
  });

  expect(states).not.toBeNull();
  expect(states.afterAtom.atomVisible).toBe(true);
  expect(states.afterAtom.bondVisible).toBe(false);
  expect(states.afterBond.atomVisible).toBe(false);
  expect(states.afterBond.bondVisible).toBe(true);
});
