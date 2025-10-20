/** @jest-environment jsdom */
import { initNewViewer } from '../public/index.js';

// This test runs in fast-jsdom mode (no real servers). We stub global.fetch to observe URLs & payloads
// and return minimal shape responses with cache_key so the client stores and reuses it.

describe.skip('cache reuse flow chooses *_from_cache endpoints [cache removed]', () => {
  const calls = [];
  const makeJsonResponse = (obj, { status = 200 } = {}) => ({
    ok: status >= 200 && status < 300,
    status,
    json: async () => obj,
    text: async () => JSON.stringify(obj),
  });
  beforeAll(() => {
    // Ensure jsdom-like window and canvas
    global.window = global.window || {};
    global.document = global.document || {
      getElementById: (id) => ({
        id,
        textContent: '',
        width: 300,
        height: 100,
        getContext: () => ({
          clearRect() {},
          beginPath() {},
          moveTo() {},
          lineTo() {},
          stroke() {},
          arc() {},
          fill() {},
          fillText() {},
        }),
      }),
      createElement: (tag) => ({ tagName: tag.toUpperCase(), getContext: () => ({}) }),
    };
    // Minimal canvas for scene.js
    const canvas = {
      getContext: () => ({}),
      addEventListener: () => {},
      getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
    };
    // Attach energy canvas element lookup
    document.getElementById = (id) => {
      if (id === 'energyCanvas')
        return {
          width: 320,
          height: 80,
          getContext: () => ({
            clearRect() {},
            beginPath() {},
            moveTo() {},
            lineTo() {},
            stroke() {},
            arc() {},
            fill() {},
            fillText() {},
          }),
        };
      if (id === 'energyLabel') return { textContent: '' };
      if (id === 'btnMDRun') return { textContent: 'run' };
      if (id === 'status') return { textContent: '' };
      return { id, getContext: () => ({}) };
    };
    global.performance = global.performance || { now: () => Date.now() };
    // Stub fetch: respond based on URL path
    global.fetch = jest.fn(async (url, init) => {
      const body = init && init.body ? JSON.parse(init.body) : null;
      calls.push({ url: String(url), body });
      if (/\/serve\/simple$/.test(url)) {
        return makeJsonResponse({
          results: {
            energy: -1.23,
            forces: [
              [0, 0, 0],
              [0, 0, 0],
            ],
          },
          cache_key: 'K1',
        });
      }
      if (/\/serve\/relax_from_cache$/.test(url)) {
        // Return new cache key each call to ensure client updates it
        return makeJsonResponse({
          initial_energy: -1.23,
          final_energy: -1.5,
          positions: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          forces: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          steps_completed: 1,
          calculator: 'uma',
          cache_key: 'K2',
        });
      }
      if (/\/serve\/relax$/.test(url)) {
        return makeJsonResponse({
          initial_energy: -1.23,
          final_energy: -1.4,
          positions: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          forces: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          steps_completed: 1,
          calculator: 'uma',
          cache_key: 'K2',
        });
      }
      if (/\/serve\/simple_from_cache$/.test(url)) {
        return makeJsonResponse({
          results: {
            energy: -1.6,
            forces: [
              [0, 0, 0],
              [0, 0, 0],
            ],
          },
          cache_key: 'K3',
        });
      }
      if (/\/serve\/md_from_cache$/.test(url)) {
        return makeJsonResponse({
          initial_energy: -1.6,
          final_energy: -1.55,
          positions: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          velocities: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          forces: [
            [0, 0, 0],
            [0, 0, 0],
          ],
          steps_completed: 1,
          temperature: 300,
          calculator: 'uma',
          cache_key: 'K4',
        });
      }
      // Default OK empty
      return makeJsonResponse({ ok: true });
    });
    // Force fast-jsdom code path in viewer init
    global.window.__MLIPVIEW_TEST_MODE = true;
    global.window.__MLIPVIEW_NO_AUTO_MD = true; // avoid auto-run
  });

  beforeEach(() => {
    calls.length = 0;
  });

  test('uses from_cache for relax and simple when no user interaction between calls', async () => {
    // Initialize viewer
    const canvas = {
      getContext: () => ({}),
      addEventListener: () => {},
      getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
    };
    const api = await initNewViewer(canvas, {
      elements: ['O', 'H'],
      positions: [
        [0, 0, 0],
        [1, 0, 0],
      ],
      bonds: [],
    });
    expect(api).toBeTruthy();

    // 1) baselineEnergy triggers /serve/simple
    expect(calls.some((c) => /\/serve\/simple$/.test(c.url))).toBe(true);

    // 2) Perform relax step; since no user interaction since baseline, it should choose relax_from_cache
    calls.length = 0;
    const r1 = await api.relaxStep();
    expect(r1 && r1.applied).toBe(true);
    expect(calls.some((c) => /\/serve\/relax_from_cache$/.test(c.url))).toBe(true);

    // 3) Now fetch forces again (computeForces). Mark local force cache stale to force a network call; should prefer simple_from_cache
    calls.length = 0;
    api.state.forceCache && (api.state.forceCache.stale = true);
    await api.ff.computeForces({ sync: true });
    expect(calls.some((c) => /\/serve\/simple_from_cache$/.test(c.url))).toBe(true);

    // 4) MD step should also use md_from_cache in this no-interaction case
    calls.length = 0;
    const md = await api.mdStep({ steps: 1, temperature: 300 });
    expect(md && md.applied).toBe(true);
    expect(calls.some((c) => /\/serve\/md_from_cache$/.test(c.url))).toBe(true);
  });
});
