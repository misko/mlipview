// Per-test setup: assert cached Ray Serve GPU health (populated in globalSetup) and log device/model load status.
// Removed GPU health assertions: tests should not depend on external server hardware state.
beforeEach(() => {});

afterEach(() => {
  try {
    if (global.window && global.window.__MLIPVIEW_DEBUG_API) {
      global.window.__MLIPVIEW_DEBUG_API = false;
      if (process.env.TEST_VERBOSE === '1') {
        // eslint-disable-next-line no-console
        console.log('[per-test-health] disabling API debug logging');
      }
    }
  } catch {}
});

afterAll(() => {
  try {
    const w = global.window;
    if (w && w.__XR_HUD_POLL_INTERVAL) {
      clearInterval(w.__XR_HUD_POLL_INTERVAL);
    }
    const cleanups = (w && w.__MLIPVIEW_CLEANUP) || global.__MLIPVIEW_CLEANUP;
    if (Array.isArray(cleanups)) {
      for (const fn of cleanups) {
        try {
          fn();
        } catch {}
      }
    }
  } catch {}
});
