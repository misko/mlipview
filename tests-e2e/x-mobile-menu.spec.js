import { test, expect } from './fixtures.js';

const MOBILE_VIEWPORT = { width: 540, height: 900 };
const MOBILE_USER_AGENT =
  'Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.0 Mobile/15E148 Safari/604.1';

test.describe('mobile responsive menu', () => {
  test.use({
    viewport: MOBILE_VIEWPORT,
    hasTouch: true,
    userAgent: MOBILE_USER_AGENT,
  });

  test('x-mobile top bar replaces desktop panel on narrow screens', async ({ page, loadViewerPage }) => {
    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      disableAutoMd: true,
    });

    const topBar = page.locator('#mobileTopBar');
    await expect(topBar).toBeVisible();
    await expect(page.locator('#controlPanel')).toBeHidden();

    const simTab = page.locator('#mobileTab-simulation');
    await simTab.click();
    await expect(topBar).toHaveAttribute('data-open', 'true');
    const sheet = topBar.locator('.mobile-sheet');
    await expect(sheet).toBeVisible();
    await expect(page.locator('#controlPanel')).toHaveAttribute('data-mobile-open', 'true');
    const childCount = await sheet.locator(':scope > *').count();
    expect(childCount).toBeGreaterThan(0);
  });
});

test.describe('desktop layout on wide screens', () => {
  test.use({
    viewport: { width: 1280, height: 900 },
    hasTouch: false,
  });

  test('x-mobile responsive reverts to desktop panel', async ({ page, loadViewerPage }) => {
    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      disableAutoMd: true,
    });
    await expect(page.locator('#mobileTopBar')).toBeHidden();
    await expect(page.locator('#controlPanel')).toBeVisible();
    await expect(page.locator('#controlPanel .panel-section').first()).toBeVisible();
  });
});
