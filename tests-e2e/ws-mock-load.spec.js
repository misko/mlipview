/**
 * E2E (browser-only): load ws-test.html with a mocked WebSocket and
 * verify protobuf stubs initialize without errors.
 *
 * What this test covers
 * - Boots the minimal harness (ws-test.html) under bundling.
 * - Mocks window.WebSocket so the page can initialize the WS client
 *   without a real backend connection.
 * - Streams browser console output into the Playwright log and collects
 *   console errors to assert against known protobuf init failures.
 * - Asserts that the bundled client API initializes (window.__WS_API__)
 *   and the page signals readiness (window.__WS_READY__).
 *
 * Expectations
 * - No protobuf initialization errors in console (e.g., missing stubs,
 *   google-protobuf issues, or goog exportSymbol errors).
 * - Page reaches a ready state with the client API exposed.
 *
 * Why this exists
 * - Provides a lightweight smoke test for the frontend protobuf wiring
 *   and module loading, independent of the backend server.
 */
// Minimal Playwright test to load ws-test.html with a mocked WebSocket
// Ensures google-protobuf and session_pb.js initialize without errors

import { test, expect } from '@playwright/test';

test.describe('WS mock page load', () => {
  test('loads ws-test.html without protobuf init errors', async ({ page, baseURL }) => {
    test.setTimeout(20000);

    // Stream full browser console to test output and collect errors for assertions
    const consoleErrors = [];
    page.on('console', (msg) => {
      const text = msg.text();
      if (msg.type() === 'error') consoleErrors.push(text);
      // Always print to Playwright runner output for debugging
      // Prefix with [browser:<type>] so it's easy to scan in CI logs
      // Avoid awaiting msg.args()/jsonValue() to keep this fast
      // eslint-disable-next-line no-console
      console.log(`[browser:${msg.type()}] ${text}`);
    });
    page.on('pageerror', (err) => {
      const text = (err && (err.message || String(err))) || 'unknown pageerror';
      consoleErrors.push(text);
      // eslint-disable-next-line no-console
      console.log(`[pageerror] ${text}`);
    });
    page.on('requestfailed', (req) => {
      const failure = (req.failure && req.failure()) || {};
      // eslint-disable-next-line no-console
      console.log(
        `[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`
      );
    });

    // Mock WebSocket so page can call ensureConnected/init without hitting a real backend
    await page.addInitScript(() => {
      // Ensure secure-ish origin behavior warning is not fatal; we just observe console
      window.__MLIPVIEW_TEST_MODE = true;
      window.__MLIPVIEW_SERVER = 'http://localhost:8000';

      class MockWS {
        constructor(url) {
          this._url = url;
          this.readyState = 0;
          setTimeout(() => {
            this.readyState = 1;
            this.onopen && this.onopen({});
          }, 10);
        }
        set binaryType(_) { }
        send(_) {
          /* ignore */
        }
        close() {
          this.readyState = 3;
          this.onclose && this.onclose({});
        }
        onopen() { }
        onmessage() { }
        onerror() { }
        onclose() { }
      }
      window.WebSocket = MockWS;
    });

    // Navigate to the test harness
    await page.goto(`${baseURL}/ws-test.html?sim=0`);

    // With bundling, we don't expose goog/jspb globals; instead, assert the client API boots
    await page.waitForFunction(() => typeof window.__WS_API__ === 'object', { timeout: 10000 });

    // The page sets __WS_READY__ = true after initializing the ws client
    await page.waitForFunction(() => !!window.__WS_READY__, { timeout: 10000 });

    // Assert no protobuf initialization error messages
    const errBlob = consoleErrors.join('\n');
    expect(errBlob).not.toMatch(
      /protobuf-missing|Failed to load google-protobuf|Failed to load session_pb|goog\.exportSymbol is not a function/
    );
  });
});
