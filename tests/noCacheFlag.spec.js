/**
 * @jest-environment jsdom
 */
/**
 * Verifies that the debug flag ?noCache=1 disables *_from_cache usage.
 */

describe('noCache debug flag forces legacy endpoints', () => {
  let fetchCalls;

  beforeEach(() => {
    // Simulate debug flag via global to avoid URL mutation
    window.__MLIP_NO_CACHE = true;
    fetchCalls = [];
    global.fetch = jest.fn(async (url, opts) => {
      fetchCalls.push(String(url));
      // Return a response with a cache_key to mimic server behavior
      const body = { results: { energy: -1.0, forces: [[0, 0, 0]] }, cache_key: 'ckey-123' };
      if (/\/serve\/relax/.test(String(url))) {
        return {
          ok: true,
          json: async () => ({
            initial_energy: -1.0,
            final_energy: -1.1,
            positions: [[0, 0, 0]],
            forces: [[0, 0, 0]],
            steps_completed: 1,
            cache_key: 'ckey-234',
          }),
        };
      }
      if (/\/serve\/md/.test(String(url))) {
        return {
          ok: true,
          json: async () => ({
            initial_energy: -1.1,
            final_energy: -1.2,
            positions: [[0, 0, 0]],
            velocities: [[0, 0, 0]],
            forces: [[0, 0, 0]],
            steps_completed: 1,
            temperature: 300,
            cache_key: 'ckey-345',
          }),
        };
      }
      return { ok: true, json: async () => body };
    });
  });

  afterEach(() => {
    delete window.__MLIP_NO_CACHE;
    jest.resetModules();
    jest.clearAllMocks();
  });

  test('forces, relax, md use non-cache endpoints', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const canvas = {
      getContext: () => ({}),
      addEventListener: () => {},
      getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
    };
    const api = await initNewViewer(canvas, { elements: ['H'], positions: [[0, 0, 0]], bonds: [] });
    expect(api).toBeTruthy();

    // Baseline energy happens during init; should include /serve/simple only
    expect(fetchCalls.some((u) => u.includes('/serve/simple'))).toBe(true);
    expect(fetchCalls.some((u) => u.includes('/serve/simple_from_cache'))).toBe(false);

    // Without any interaction, relax should still use /serve/relax (not from_cache)
    fetchCalls.length = 0;
    const r = await api.relaxStep();
    expect(r && r.applied).toBe(true);
    expect(fetchCalls.some((u) => u.includes('/serve/relax'))).toBe(true);
    expect(fetchCalls.some((u) => u.includes('/serve/relax_from_cache'))).toBe(false);

    // MD should also use /serve/md
    fetchCalls.length = 0;
    const m = await api.mdStep({ steps: 1, temperature: 300 });
    expect(m && m.applied).toBe(true);
    expect(fetchCalls.some((u) => u.includes('/serve/md'))).toBe(true);
    expect(fetchCalls.some((u) => u.includes('/serve/md_from_cache'))).toBe(false);
  });
});
