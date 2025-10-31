import { test, expect } from './fixtures.js';
import {
  waitForTimelineBuffer,
  computeTimelineOffset,
  selectTimelineOffset,
  waitForTimelineState,
} from './utils/timeline.js';

test.describe('timeline mesh mode transitions', () => {
  test('fades migrate instances to soft masters and restore on clear', async ({ page, loadViewerPage }) => {
    test.setTimeout(150_000);

    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

    await page.evaluate(() => {
      window.viewerApi.startMDContinuous({ steps: 240, temperature: 1600 });
    });

    await waitForTimelineBuffer(page, 40, { label: 'mesh-buffer' });

    const offset = await computeTimelineOffset(page, { mode: 'middle' });
    expect(offset).not.toBeNull();

    await selectTimelineOffset(page, offset, { label: 'mesh-select' });
    await waitForTimelineState(
      page,
      () => {
        const status =
          window.viewerApi?.timeline?.getStatus?.() ??
          window.viewerApi?.timeline?.getState?.() ??
          null;
        return !!status && !!status.active;
      },
      null,
      { label: 'mesh-active' }
    );

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

    await page.evaluate(async () => {
      await (window.viewerApi.timeline.live?.() || Promise.resolve());
    });
    await waitForTimelineState(
      page,
      () => {
        const status =
          window.viewerApi?.timeline?.getStatus?.() ??
          window.viewerApi?.timeline?.getState?.() ??
          null;
        return !!status && !status.active;
      },
      null,
      { label: 'mesh-return-live' }
    );

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
