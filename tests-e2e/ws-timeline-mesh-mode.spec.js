import { test, expect } from './fixtures.js';

async function waitForTimelineBuffer(page, targetSize) {
  await page.waitForFunction(
    (target) => {
      const stats = window.viewerApi.timeline.bufferStats();
      return stats && stats.size >= target;
    },
    targetSize,
    { timeout: 45_000 }
  );
}

async function pickHistoricalOffset(page) {
  return await page.evaluate(() => {
    const offsets = window.viewerApi.timeline.getOffsets();
    if (!offsets || !offsets.length) return null;
    const candidates = offsets.filter((o) => o < -1);
    if (!candidates.length) return offsets[offsets.length - 1];
    return candidates[Math.floor(candidates.length / 2)];
  });
}

test.describe('timeline mesh mode transitions', () => {
  test('fades migrate instances to soft masters and restore on clear', async ({ page, loadViewerPage }) => {
    test.setTimeout(150_000);

    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

    await page.evaluate(() => {
      window.viewerApi.startMDContinuous({ steps: 240, temperature: 1600 });
    });

    await waitForTimelineBuffer(page, 40);

    const offset = await pickHistoricalOffset(page);
    expect(offset).not.toBeNull();

    await page.evaluate((value) => window.viewerApi.timeline.select(value), offset);
    await page.waitForFunction(() => window.viewerApi.timeline.getState().active, null, { timeout: 12_000 });

    const baselineModes = await page.evaluate(() => ({
      atoms: window.viewerApi.view.getAtomMeshModes(),
      bonds: window.viewerApi.view.getBondMeshModes(),
    }));
    const baselineAtomSoft = baselineModes.atoms.current.filter((mode) => mode === 'soft').length;
    const baselineBondSoft = baselineModes.bonds.current.filter((mode) => mode === 'soft').length;

    await page.evaluate(() => {
      window.viewerApi.view.setOpacityMask({
        focus: { atoms: [] },
        mode: { focus: 'solid', background: 'soft' },
        applyTo: { primaryAtoms: true, primaryBonds: true },
      });
    });

    await page.waitForFunction(
      () => window.viewerApi.view.getBondMeshModes().current.some((mode) => mode === 'soft'),
      null,
      { timeout: 5_000 }
    );
    await page.waitForFunction(
      () => window.viewerApi.view.getAtomMeshModes().current.some((mode) => mode === 'soft'),
      null,
      { timeout: 5_000 }
    );

    const fadedModes = await page.evaluate(() => ({
      atoms: window.viewerApi.view.getAtomMeshModes().current,
      bonds: window.viewerApi.view.getBondMeshModes().current,
    }));
    const fadedAtomSoft = fadedModes.atoms.filter((mode) => mode === 'soft').length;
    const fadedBondSoft = fadedModes.bonds.filter((mode) => mode === 'soft').length;
    expect(fadedAtomSoft).toBeGreaterThanOrEqual(baselineAtomSoft);
    expect(fadedBondSoft).toBeGreaterThanOrEqual(baselineBondSoft);
    expect(fadedAtomSoft).toBeGreaterThan(0);
    expect(fadedBondSoft).toBeGreaterThan(0);

    await page.evaluate(() => window.viewerApi.view.setOpacityMask(null));

    await page.waitForFunction(
      ({ atomSoftCount, bondSoftCount }) => {
        const atomCurrent = window.viewerApi.view.getAtomMeshModes().current;
        const bondCurrent = window.viewerApi.view.getBondMeshModes().current;
        const atomSoft = atomCurrent.filter((mode) => mode === 'soft').length;
        const bondSoft = bondCurrent.filter((mode) => mode === 'soft').length;
        return atomSoft === atomSoftCount && bondSoft === bondSoftCount;
      },
      { atomSoftCount: baselineAtomSoft, bondSoftCount: baselineBondSoft },
      { timeout: 5_000 }
    );

    await page.evaluate(() => window.viewerApi.timeline.live());
    await page.waitForFunction(() => !window.viewerApi.timeline.getState().active, null, { timeout: 12_000 });

    const finalModes = await page.evaluate(() => ({
      atoms: window.viewerApi.view.getAtomMeshModes().current,
      bonds: window.viewerApi.view.getBondMeshModes().current,
    }));
    const finalAtomSoft = finalModes.atoms.filter((mode) => mode === 'soft').length;
    const finalBondSoft = finalModes.bonds.filter((mode) => mode === 'soft').length;
    expect(finalAtomSoft).toBe(baselineAtomSoft);
    expect(finalBondSoft).toBe(baselineBondSoft);

    await page.evaluate(() => window.viewerApi.stopSimulation());
  });
});
