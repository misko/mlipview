import { test as base, expect } from '@playwright/test';

const DEFAULT_SERVER_URL = 'http://127.0.0.1:8000';

function attachLogging(context) {
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

  try {
    const pages = typeof context.pages === 'function' ? context.pages() : [];
    for (const p of pages) attachListeners(p);
  } catch { }

  context.on('page', (page) => attachListeners(page));
}

const DEFAULT_RUNTIME_CONFIG = {
  wsConnectTimeoutMs: 15000,
};

async function applyRuntimeInit(page, {
  server = DEFAULT_SERVER_URL,
  testMode = true,
  disableAutoMd = true,
  configOverrides = null,
} = {}) {
  const mergedOverrides = { ...DEFAULT_RUNTIME_CONFIG, ...(configOverrides || {}) };
  await page.addInitScript(
    (
      { server, testMode, disableAutoMd, configOverrides }
    ) => {
      if (server) window.__MLIPVIEW_SERVER = server;
      if (typeof testMode === 'boolean') window.__MLIPVIEW_TEST_MODE = !!testMode;
      if (disableAutoMd) window.__MLIPVIEW_NO_AUTO_MD = true;
      window.__MLIP_CONFIG = window.__MLIP_CONFIG || {};
      if (configOverrides && typeof configOverrides === 'object') {
        Object.assign(window.__MLIP_CONFIG, configOverrides);
      }
    },
    { server, testMode, disableAutoMd, configOverrides: mergedOverrides }
  );
}

function buildQuery(defaults, overrides) {
  const params = new URLSearchParams();
  const merged = { ...defaults, ...(overrides || {}) };
  Object.entries(merged).forEach(([key, value]) => {
    if (value === undefined || value === null) return;
    params.set(key, String(value));
  });
  return params.toString();
}

// Shared fixture to mirror browser console, page errors, and failed requests
// to the test runner output for all tests that import from this module.
export const test = base.extend({
  context: async ({ context }, use) => {
    attachLogging(context);
    await use(context);
  },

  loadViewerPage: async ({ page, baseURL }, use) => {
  const loader = async (
    {
      query = {},
      server = DEFAULT_SERVER_URL,
      testMode = true,
      disableAutoMd = true,
      configOverrides = null,
      waitForViewer = true,
      extraInit = null,
    } = {}
  ) => {
      await applyRuntimeInit(page, {
        server,
        testMode,
        disableAutoMd,
        configOverrides,
      });

      if (typeof extraInit === 'function') {
        await page.addInitScript(extraInit);
      }

      const qs = buildQuery({ debug: 1 }, query);
      const target = `${baseURL || ''}/index.html${qs ? `?${qs}` : ''}`;
      await page.goto(target);

      if (waitForViewer) {
        await page.waitForFunction(
          () => !!window.viewerApi && window.__MLIP_DEFAULT_LOADED === true,
          null,
          { timeout: 45_000 }
        );
        await page.evaluate(async () => {
          const ws = window.__fairchem_ws__ || window.__WS_API__;
          if (!ws || typeof ws.ensureConnected !== 'function') return;
          const tryConnect = ws.ensureConnected({ timeoutMs: 8000 }).catch(() => false);
          await Promise.race([
            tryConnect,
            new Promise(resolve => setTimeout(resolve, 8000, false)),
          ]);
          try { await window.viewerApi?.baselineEnergy?.(); } catch { }
          try { ws.waitForEnergy({ timeoutMs: 10000 }).catch(() => { }); } catch { }
        });
      }
    };

    await use(loader);
  },

  loadWsHarnessPage: async ({ page, baseURL }, use) => {
    const loader = async (
      {
        query = {},
        server = DEFAULT_SERVER_URL,
        testMode = true,
        disableAutoMd = true,
        configOverrides = null,
        waitForReady = true,
        extraInit = null,
      } = {}
    ) => {
      await applyRuntimeInit(page, {
        server,
        testMode,
        disableAutoMd,
        configOverrides,
      });

      if (typeof extraInit === 'function') {
        await page.addInitScript(extraInit);
      }

      const qs = buildQuery({ debug: 1 }, query);
      const target = `${baseURL || ''}/ws-test.html${qs ? `?${qs}` : ''}`;
      await page.goto(target);

      if (waitForReady) {
        await page.waitForFunction(() => !!window.__WS_READY__, null, { timeout: 15_000 });
      }
    };

    await use(loader);
  },

  ensureWsReady: async ({ page }, use) => {
    const ensure = async ({ timeoutMs = 20_000 } = {}) => {
      return await page.evaluate(async ({ timeoutMs }) => {
        const ws = window.__fairchem_ws__ || window.__WS_API__;
        if (!ws || typeof ws.ensureConnected !== 'function') return false;
        const safeTimeout = Math.max(1000, Math.min(timeoutMs || 20000, 8000));
        const connect = ws.ensureConnected({ timeoutMs: safeTimeout })
          .then(() => true)
          .catch(() => false);
        return await Promise.race([
          connect,
          new Promise((resolve) => setTimeout(() => resolve(false), safeTimeout)),
        ]);
      }, { timeoutMs });
    };

    await use(ensure);
  },
});

export { expect };
