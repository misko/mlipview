import { test as base, expect } from '@playwright/test';

// Shared fixture to mirror browser console, page errors, and failed requests
// to the test runner output for all tests that import from this module.
export const test = base.extend({
  context: async ({ context }, use) => {
    const attachListeners = (page) => {
      page.on('console', (msg) => {
        // eslint-disable-next-line no-console
        console.log(`[browser:${msg.type()}] ${msg.text()}`);
      });
      page.on('pageerror', (err) => {
        // eslint-disable-next-line no-console
        console.log(`[pageerror] ${err?.message || String(err)}`);
      });
      page.on('requestfailed', (req) => {
        const failure = (req.failure && req.failure()) || {};
        // eslint-disable-next-line no-console
        console.log(
          `[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`
        );
      });
    };

    // Attach to any pages already created for this context
    try {
      const pages = typeof context.pages === 'function' ? context.pages() : [];
      for (const p of pages) attachListeners(p);
    } catch { }

    // Attach for all subsequently opened pages
    context.on('page', (page) => attachListeners(page));

    await use(context);
  },
});

export { expect };
