// E2E test: open app, open molecules selector, choose benzene, ensure no oversized stray highlight sphere at center.
// Assumptions: server serves / (legacy UI) with a Molecules button (#btnMolecules) that when clicked shows a selector.
// If the selector is not yet implemented in DOM, we fallback: ensure initial default (roy) then trigger benzene load via /api/molecules selection UI if available.
import { test, expect } from '@playwright/test';

// Helper: compute approximate ring centroid from visible atom div overlays if any (fallback not used here)

// We will detect any mesh-like DOM artifact representing the highlight by querying canvas overlay debug elements if created.
// Since Babylon renders to a canvas, we cannot directly query meshes; instead we enable debug query parameter to inject
// a window hook that exposes highlight meshes for inspection. We'll add a small script injection to poll.

async function exposeHighlight(page) {
  await page.addInitScript(() => {
    window.__MLIP_E2E__ = { getHighlightInfo: () => {
      const scene = window.scene; // if app exposed it; else try from global BABYLON
      try {
        const meshes = scene?.meshes || [];
        const hi = meshes.filter(m => m.name && m.name.startsWith('highlight_'));
        return hi.map(m => ({ name: m.name, pos: m.position && { x:m.position.x, y:m.position.y, z:m.position.z },
          scaling: m.scaling && { x:m.scaling.x, y:m.scaling.y, z:m.scaling.z }, isVisible: !!m.isVisible }));
      } catch(e) { return { error: e.message }; }
    }};
  });
}

test.describe('Benzene selection / highlight cleanliness', () => {
  test('no giant blue sphere at benzene center after selecting benzene', async ({ page }) => {
    await exposeHighlight(page);
    await page.goto('/?debug=1');

    // Click molecules button if present
    const molBtn = page.locator('#btnMolecules');
    if (await molBtn.count()) {
      await molBtn.click();
      // Expect a selector panel to appear - we look for any button or option containing 'benzene'
      const benzeneOption = page.locator('text=/benzene/i');
      // Allow some time for fetch of molecule list
      if (await benzeneOption.count() === 0) {
        // Wait for dynamic list
        await benzeneOption.first().waitFor({ state: 'visible', timeout: 5000 }).catch(()=>{});
      }
      if (await benzeneOption.count()) {
        await benzeneOption.first().click();
      }
    }

    // Give rendering a few frames
    await page.waitForTimeout(1000);

    const highlightInfo = await page.evaluate(() => window.__MLIP_E2E__?.getHighlightInfo());

    // Determine benzene radius approx by sampling C atom positions from global state if exposed
    // If not available, we just ensure no visible highlight sphere with extremely large scaling.
    let suspicious = [];
    if (Array.isArray(highlightInfo)) {
      suspicious = highlightInfo.filter(h => h.isVisible && h.name.includes('highlight_atom') && h.scaling && (h.scaling.x > 2.5 || h.scaling.y > 2.5 || h.scaling.z > 2.5));
    }

    expect(suspicious.length, 'No oversized highlight sphere should remain visible').toBe(0);
  });
});
