import { test, expect } from './fixtures.js';

async function exposeHighlightState(page) {
  await page.addInitScript(() => {
    window.__MLIP_E2E__ = window.__MLIP_E2E__ || {};
    window.__MLIP_E2E__.getHighlightState = () => {
      try {
        return typeof window.viewerApi?.debugHighlightState === 'function'
          ? window.viewerApi.debugHighlightState()
          : null;
      } catch {
        return null;
      }
    };
  });
}

test.describe('x-benzene-selection', () => {
  test('loading benzene leaves highlight hidden', async ({ page }) => {
    await exposeHighlightState(page);
    await page.goto('/?debug=1');

    await page.waitForFunction(() => !!window.viewerApi, { timeout: 10000 });

    // Open molecule selector via UI if present and choose benzene; fallback to API.
    const molBtn = page.locator('#btnMolecules');
    if (await molBtn.count()) {
      await molBtn.click();
      const benzeneOption = page.locator('text=/benzene/i').first();
      await benzeneOption.waitFor({ state: 'visible', timeout: 5000 }).catch(() => {});
      if (await benzeneOption.count()) {
        await benzeneOption.click();
      }
    } else {
      await page.evaluate(async () => {
        try {
          const resp = await fetch('/molecules/benzene.xyz');
          const txt = await resp.text();
          const { parseXYZ } = await import('../public/util/xyzLoader.js');
          const { applyParsedToViewer } = await import('../public/util/moleculeLoader.js');
          const parsed = parseXYZ(txt);
          applyParsedToViewer(window.viewerApi, parsed);
        } catch {}
      });
    }

    await page.waitForTimeout(1000);

    const highlight = await page.evaluate(
      () => window.__MLIP_E2E__?.getHighlightState && window.__MLIP_E2E__.getHighlightState()
    );
    expect(highlight).toBeTruthy();
    expect(highlight.atomVisible).toBe(false);
    expect(highlight.bondVisible).toBe(false);
  });
});

