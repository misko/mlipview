// Purpose: Validate safe sphere clamping keeps atoms within a selectable radius when positions overshoot.
import { test, expect } from './fixtures.js';

test('Safe sphere clamps manual and streamed runaway positions', async ({
    page,
    loadViewerPage,
    ensureWsReady,
}) => {
    test.setTimeout(60_000);

    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await loadViewerPage({
        query: { mol: 'molecules/water.xyz' },
        configOverrides: { safeSphereRadius: 4 },
    });

    expect(await ensureWsReady()).toBeTruthy();

    const manualClamp = await page.evaluate(() => {
        const radius = window.__MLIP_CONFIG.safeSphereRadius;
        const state = window.viewerApi.state;
        const idx = 0;
        const target = radius * 10;
        state.positions[idx].x = target;
        state.positions[idx].y = target * 0.5;
        state.positions[idx].z = -target * 0.25;
        const before = Math.hypot(state.positions[idx].x, state.positions[idx].y, state.positions[idx].z);
        state.markPositionsChanged();
        const after = Math.hypot(state.positions[idx].x, state.positions[idx].y, state.positions[idx].z);
        return { radius, before, after };
    });
    expect(manualClamp.before).toBeGreaterThan(manualClamp.radius);
    expect(Math.abs(manualClamp.after - manualClamp.radius)).toBeLessThan(1e-6);

    const streamClamp = await page.evaluate(async () => {
        const radius = window.__MLIP_CONFIG.safeSphereRadius;
        const ws = window.__fairchem_ws__;
        const state = window.viewerApi.state;
        const N = state.positions.length;
        const far = radius * 25;
        const positions = Array.from({ length: N }, (_, i) => [far, far * (i + 1), far * -0.5]);
        ws.injectTestResult({ seq: 999, positions, energy: 0, userInteractionCount: 0 });
        await new Promise((resolve) => setTimeout(resolve, 100));
        const dists = state.positions.map((p) => Math.hypot(p.x, p.y, p.z));
        return { radius, max: Math.max(...dists), min: Math.min(...dists) };
    });
    expect(streamClamp.max).toBeLessThanOrEqual(streamClamp.radius + 1e-6);
    expect(streamClamp.min).toBeGreaterThan(0);
});
