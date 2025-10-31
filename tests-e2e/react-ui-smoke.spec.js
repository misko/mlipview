import { test, expect } from './fixtures.js';

test.describe('React UI smoke', () => {
  test('loads viewer without console errors', async ({ page, loadViewerPage }) => {
    const errors = [];
    const ignoredErrorPatterns = [/WebXR can only be served over HTTPS/i];

    page.on('console', (msg) => {
      if (msg.type() === 'error') {
        errors.push(msg.text());
      }
    });
    page.on('pageerror', (err) => {
      errors.push(err?.message || String(err));
    });

    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', ui: 'react' },
      waitForViewer: true,
    });

    const status = await page.evaluate(() => ({
      hasViewerApi: typeof window.viewerApi === 'object' && !!window.viewerApi,
      defaultLoaded: window.__MLIP_DEFAULT_LOADED === true,
    }));

    const relevantErrors = errors.filter(
      (msg) => !ignoredErrorPatterns.some((pattern) => pattern.test(msg)),
    );

    expect(status.hasViewerApi).toBe(true);
    expect(status.defaultLoaded).toBe(true);
    expect(relevantErrors).toEqual([]);
  });
});
