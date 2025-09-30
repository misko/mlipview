import { test, expect } from '@playwright/test';

test.describe('XR dropdown presence', () => {
  test('XR mode select exists with expected options', async ({ page }) => {
    await page.goto('/');
    await page.waitForSelector('#xrModeSelect');
    const options = await page.$$eval('#xrModeSelect option', els => els.map(e=>e.value));
    expect(options).toEqual(expect.arrayContaining(['none','vr','ar']));
    const current = await page.$eval('#xrModeSelect', el=>el.value);
    expect(['none','vr','ar']).toContain(current);
  });
});
