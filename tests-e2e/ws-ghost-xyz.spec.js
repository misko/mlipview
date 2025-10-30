import { Buffer } from 'node:buffer';
import { test, expect } from './fixtures.js';

function encodeXYZ(text) {
  return Buffer.from(text, 'utf8').toString('base64');
}

async function waitForBondState(page) {
  await page.waitForFunction(() => {
    const api = window.viewerApi;
    if (!api) return false;
    const st = api.state;
    return Array.isArray(st?.bonds) && Array.isArray(st?.ghostBondMeta);
  });

  // Re-run the recompute to make sure logs/ghosts are up to date.
  await page.evaluate(() => {
    try {
      window.viewerApi.recomputeBonds?.();
      window.viewerApi.view?.rebuildGhosts?.();
    } catch (err) {
      console.warn('[ws-ghost-xyz][recompute]', err);
    }
  });

  return await page.evaluate(() => {
    const { bonds, ghostBondMeta, showGhostCells } = window.viewerApi.state;
    return {
      bondCount: Array.isArray(bonds) ? bonds.length : 0,
      ghostCount: Array.isArray(ghostBondMeta) ? ghostBondMeta.length : 0,
      ghostsVisible: !!showGhostCells,
    };
  });
}

test.describe('XYZ periodic ghost bonds', () => {
  test('atoms near opposing faces yield two ghost bonds', async ({ page, loadViewerPage }) => {
    const xyz = [
      '2',
      'Lattice=10,0,0;0,10,0;0,0,10 origin=-5,-5,-5',
      'C   9.500000   0.000000   0.000000',
      'C  -9.500000   0.000000   0.000000',
      '',
    ].join('\n');

    await loadViewerPage({
      query: { autoMD: 0, molxyz: encodeXYZ(xyz) },
      disableAutoMd: true,
    });

    const summary = await waitForBondState(page);
    expect(summary.bondCount).toBe(0);
    expect(summary.ghostCount).toBe(2);
    expect(summary.ghostsVisible).toBe(true);
  });

  test('adjacent atoms yield one primary bond plus six ghosts', async ({ page, loadViewerPage }) => {
    const xyz = [
      '2',
      'Lattice=10,0,0;0,10,0;0,0,10 origin=-5,-5,-5',
      'C   4.500000   0.000000   0.000000',
      'C   5.500000   0.000000   0.000000',
      '',
    ].join('\n');

    await loadViewerPage({
      query: { autoMD: 0, molxyz: encodeXYZ(xyz) },
      disableAutoMd: true,
    });

    const summary = await waitForBondState(page);
    expect(summary.bondCount).toBe(1);
    expect(summary.ghostCount).toBe(6);
    expect(summary.ghostsVisible).toBe(true);
  });
});
