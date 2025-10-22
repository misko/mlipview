import { test, expect } from './fixtures.js';

// Deterministic energy-step increment test.
// Rationale: Pointer-based drag in headless browsers can be flaky due to lack of real
// GPU/layout timing and the ray->plane projection producing zero deltas. Instead we:
//  1. Wait for viewer initialization & default molecule load.
//  2. Programmatically select atom 0 via selection service.
//  3. Mutate its x coordinate by a small amount and call markPositionsChanged().
//     This triggers the debounced positionsChanged listener in index.js which
//     recomputes forces and records an interaction ('posChange'), incrementing steps.
// This keeps the test focused on the contract: geometry change -> energy steps increase.

test.describe('Energy plot updates on atom drag (full page)', () => {
  test('mutating an atom position increments energy steps', async ({ page }) => {
    await page.goto('/?debug=1');
    await page.waitForFunction(() => !!window._viewer, null, { timeout: 10000 });
    await page.waitForSelector('#energyLabel');

    const parseSteps = (txt) => {
      const m = txt.match(/steps=(\d+)/);
      return m ? +m[1] : 0;
    };
    const initialSteps = parseSteps(await page.locator('#energyLabel').innerText());

    // Perform deterministic position mutation inside the page context
    await page.evaluate(() => {
      const v = window._viewer;
      if (!v || !v.state?.positions?.length) return;
      // Select first atom using public selection service API if available
      try {
        v.selection.clickAtom?.(0);
      } catch { }
      // Shift x coordinate slightly
      v.state.positions[0].x += 0.4;
      // Notify system of position change; energy debounced listener will fire
      v.state.markPositionsChanged();
    });

    // Allow debounce (50ms) + force compute + draw
    let afterSteps = initialSteps;
    for (let t = 0; t < 12; t++) {
      await page.waitForTimeout(50);
      afterSteps = parseSteps(await page.locator('#energyLabel').innerText());
      if (afterSteps > initialSteps) break;
    }
    expect(afterSteps).toBeGreaterThan(initialSteps);
  });
});
