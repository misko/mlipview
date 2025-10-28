import { test, expect } from './fixtures.js';

const MOBILE_VIEWPORT = { width: 540, height: 900 };
const MOBILE_USER_AGENT =
  'Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.0 Mobile/15E148 Safari/604.1';

test.describe('mobile wiring smoke', () => {
  test.use({
    viewport: MOBILE_VIEWPORT,
    hasTouch: true,
    userAgent: MOBILE_USER_AGENT,
  });

  test('x-mobile viewer wiring initializes touch controls and top bar', async ({
    page,
    loadViewerPage,
  }) => {
    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      disableAutoMd: true,
    });

    const flags = await page.evaluate(() => ({
      hasApi: typeof window.viewerApi === 'object',
      touchInstalled: !!window.__MLIPVIEW_TOUCH_INSTALLED,
      topBarExists: !!document.getElementById('mobileTopBar'),
      panelMode: document.getElementById('controlPanel')?.getAttribute('data-mode') || null,
    }));

    expect(flags.hasApi).toBe(true);
    expect(flags.touchInstalled).toBe(true);
    expect(flags.topBarExists).toBe(true);
    expect(flags.panelMode).toBe('mobile');
  });
});
