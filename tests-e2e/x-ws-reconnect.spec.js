import { test, expect } from './fixtures.js';

const BACKEND_URL = 'http://127.0.0.1:8000';

async function ensureCudaBackend(page) {
  const health = await page.request.get(`${BACKEND_URL}/serve/health`);
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const info = await health.json();
  if (!info || String(info.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');
}

test.describe('WebSocket reconnect handling', () => {
  test('shows reconnect banner with countdown and manual reconnect', async ({ page, loadViewerPage, ensureWsReady }) => {
    test.setTimeout(60_000);
    await ensureCudaBackend(page);

    await loadViewerPage({
      server: BACKEND_URL,
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      testMode: true,
      disableAutoMd: true,
    });

    expect(await ensureWsReady()).toBeTruthy();

    const bannerLocator = page.locator('[data-banner="ws-reconnect"]');
    const countdownLocator = page.locator('#wsReconnectCountdown');

    await page.evaluate(() => {
      window.viewerApi.simulateWsDrop({ failAttempts: 2 });
    });

    await page.waitForFunction(() => window.viewerApi.debugWsState().reconnect.bannerVisible === true, { timeout: 20000 });
    await expect(bannerLocator).toBeVisible({ timeout: 20000 });
    const initial = await countdownLocator.innerText();
    await page.waitForTimeout(1200);
    const later = await countdownLocator.innerText();
    expect(Number.parseInt(later.replace(/\D+/g, ''), 10)).toBeLessThanOrEqual(Number.parseInt(initial.replace(/\D+/g, ''), 10));

    await page.click('#wsReconnectNow');
    await page.waitForFunction(() => window.viewerApi.debugWsState().reconnect.bannerVisible === false, { timeout: 20000 });

    await page.waitForFunction(() => {
      const info = window.viewerApi.debugWsState();
      return info.log.some((evt) => evt.type === 'open' && evt.isReconnect === true);
    });
  });

  test('recovers an in-flight MD run after a simulated server restart', async ({ page, loadViewerPage, ensureWsReady }) => {
    test.setTimeout(60_000);
    await ensureCudaBackend(page);

    await loadViewerPage({
      server: BACKEND_URL,
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      testMode: true,
      disableAutoMd: true,
    });

    expect(await ensureWsReady()).toBeTruthy();

    await page.evaluate(() => {
      const ws = window.__fairchem_ws__;
      window.__MLIPVIEW_WS_FRAMES = [];
      window.__WS_FRAME_OFF__?.();
      window.__WS_FRAME_OFF__ = ws.onResult((res) => {
        try {
          if (res && typeof res.seq === 'number') {
            window.__MLIPVIEW_WS_FRAMES.push({ seq: res.seq, simStep: res.simStep, timestamp: Date.now() });
          }
        } catch {}
      });
      ws.setTestHook?.((msg) => {
        window.__MLIPVIEW_WS_OUTBOUND ||= [];
        window.__MLIPVIEW_WS_OUTBOUND.push({ ...msg, timestamp: Date.now() });
      });
    });

    await page.evaluate(async () => {
      await window.viewerApi.startMDContinuous({ steps: 400, temperature: 900 });
    });

    await page.waitForFunction(() => (window.__MLIPVIEW_WS_FRAMES || []).filter((f) => Number.isFinite(f.simStep)).length >= 8, { timeout: 20000 });

    const before = await page.evaluate(() => ({
      outbound: (window.__MLIPVIEW_WS_OUTBOUND || []).filter((m) => m.type === 'USER_INTERACTION').length,
      frames: (window.__MLIPVIEW_WS_FRAMES || []).length,
    }));
    const beforeOutbound = before.outbound;
    const beforeFrames = before.frames;

    const startTime = Date.now();
    await page.evaluate(() => {
      window.viewerApi.simulateWsDrop({ failAttempts: 1 });
    });

    await page.waitForFunction((startedAt) => {
      const info = window.viewerApi.debugWsState();
      return info.log.some((evt) => evt.type === 'open' && evt.isReconnect === true && evt.timestamp > startedAt);
    }, startTime, { timeout: 20000 });

    await page.waitForFunction((prevOutbound) => {
      const outbound = (window.__MLIPVIEW_WS_OUTBOUND || []).filter((m) => m.type === 'USER_INTERACTION').length;
      return outbound > prevOutbound;
    }, beforeOutbound, { timeout: 20000 });

    await page.waitForFunction((prevFrames) => {
      const frames = (window.__MLIPVIEW_WS_FRAMES || []).length;
      return frames >= prevFrames + 5;
    }, beforeFrames, { timeout: 20000 });

    await page.waitForFunction(() => window.viewerApi.getMetrics()?.running === 'md', { timeout: 10000 });

    await page.evaluate(() => {
      window.__WS_FRAME_OFF__?.();
      window.__WS_FRAME_OFF__ = null;
    });
  });

  test('limits connect attempts over 8 seconds for unreachable server', async ({ page, loadViewerPage }) => {
    test.setTimeout(60_000);

    await loadViewerPage({
      server: 'http://127.0.0.1:65532',
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      testMode: true,
      disableAutoMd: true,
      waitForViewer: false,
      extraInit: () => {
        window.__WS_STATE_TEST_LOG__ = [];
        const prev = typeof window.__WS_ON_STATE__ === 'function' ? window.__WS_ON_STATE__ : null;
        window.__WS_ON_STATE__ = function patchedWsState(evt) {
          window.__WS_STATE_TEST_LOG__.push({ ...(evt || {}), timestamp: Date.now() });
          if (prev) {
            try { prev(evt); } catch { }
          }
        };
      },
    });

    await page.waitForFunction(() => Array.isArray(window.__WS_STATE_TEST_LOG__) && window.__WS_STATE_TEST_LOG__.some((evt) => evt.type === 'error' || evt.type === 'close'), { timeout: 30000 });

    await page.waitForTimeout(8200);

    const summary = await page.evaluate(() => {
      const log = window.__WS_STATE_TEST_LOG__ || [];
      const connectingEvents = log.filter((evt) => evt.type === 'connecting').map((evt) => ({ attempts: evt.attempts || evt.details?.attempts || evt.attempt || evt.value || 0, timestamp: evt.timestamp }));
      const scheduled = log.filter((evt) => evt.type === 'reconnect-scheduled');
      return {
        connectingCount: connectingEvents.length,
        connectingAttempts: connectingEvents.map((evt) => evt.attempts),
        firstConnectingTs: connectingEvents.length ? connectingEvents[0].timestamp : null,
        lastConnectingTs: connectingEvents.length ? connectingEvents[connectingEvents.length - 1].timestamp : null,
        maxAttempts: connectingEvents.length ? Math.max(...connectingEvents.map((evt) => evt.attempts)) : 0,
        scheduledCount: scheduled.length,
      };
    });

    expect(summary.connectingCount).toBeLessThanOrEqual(3);
    expect(summary.connectingCount).toBeGreaterThanOrEqual(2);
    if (summary.firstConnectingTs != null && summary.lastConnectingTs != null) {
      expect(summary.lastConnectingTs - summary.firstConnectingTs).toBeLessThanOrEqual(8200);
    }
    expect(summary.maxAttempts).toBeLessThanOrEqual(3);
    expect(summary.scheduledCount).toBeGreaterThanOrEqual(1);
  });
});
