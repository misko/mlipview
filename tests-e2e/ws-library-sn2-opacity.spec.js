import { test, expect } from './fixtures.js';

async function waitForOffsets(page, minCount) {
  await page.waitForFunction(
    (need) => {
      const api = window.viewerApi;
      if (!api?.timeline?.getOffsets) return false;
      const offsets = api.timeline.getOffsets();
      return Array.isArray(offsets) && offsets.length >= need;
    },
    minCount,
    { timeout: 20_000 }
  );
}

test.describe('library SN2 opacity', () => {
  test('earliest frames keep atoms opaque', async ({ page, loadViewerPage }) => {
    test.setTimeout(60_000);

    await loadViewerPage({ disableAutoMd: false });

    await page.evaluate(async () => {
      const api = window.viewerApi;
      if (!api?.session?.loadFromLibrary) throw new Error('loadFromLibrary unavailable');
      await api.session.loadFromLibrary('sn2');
      api.timeline.pause();
    });

    await waitForOffsets(page, 4);

    const earliestOffsets = await page.evaluate(() => {
      const api = window.viewerApi;
      const offsets = api.timeline.getOffsets();
      if (!Array.isArray(offsets) || offsets.length < 4) return [];
      return offsets.slice(-4);
    });

    expect(Array.isArray(earliestOffsets)).toBe(true);
    expect(earliestOffsets.length).toBeGreaterThanOrEqual(4);

    for (const offset of earliestOffsets) {
      await page.evaluate((value) => {
        window.viewerApi.timeline.select(value);
      }, offset);

      await page.waitForFunction(
        (expected) => window.viewerApi.timeline.getState().offset === expected,
        offset,
        { timeout: 10_000 }
      );

      const { atomModes, offsetConfirmed } = await page.evaluate(() => {
        const api = window.viewerApi;
        const modes = api.view?.getAtomMeshModes?.();
        return {
          atomModes: Array.isArray(modes?.current) ? modes.current.slice() : [],
          offsetConfirmed: api.timeline.getState().offset,
        };
      });

      expect(offsetConfirmed).toBe(offset);
      expect(atomModes.length).toBeGreaterThan(0);
      expect(atomModes.every((mode) => mode === 'solid')).toBe(true);
    }
  });
});
