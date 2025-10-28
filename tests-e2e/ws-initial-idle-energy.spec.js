// Purpose: Validate that on page load with ?autoMD=0 the client
// sends an initial USER_INTERACTION and the server responds with
// an idle energy frame.
import { test, expect } from './fixtures.js';

async function waitForEnergy(page, timeoutMs = 12000) {
    return await page.evaluate(
        async ({ timeoutMs }) => {
            const ws = window.__fairchem_ws__;
            if (!ws || typeof ws.waitForEnergy !== 'function') return null;
            try {
                await ws.ensureConnected();
            } catch { }
            try {
                const res = await ws.waitForEnergy({ timeoutMs });
                return res || null;
            } catch {
                return null;
            }
        },
        { timeoutMs }
    );
}

test('initial load (autoMD=0) produces idle energy frame', async ({ page, loadViewerPage }) => {
    test.setTimeout(45000);
    // Collect only error-level logs for assertions (stdout mirroring handled by fixtures)
    const consoleErrors = [];
    page.on('console', (msg) => {
        const text = msg.text();
        if (msg.type() === 'error') consoleErrors.push(text);
    });
    page.on('pageerror', (err) => {
        const text = (err && (err.message || String(err))) || 'unknown pageerror';
        consoleErrors.push(text);
    });

    await loadViewerPage({ query: { autoMD: 0, wsDebug: 1 }, testMode: false });

    // Wait for the first energy-bearing frame (idle compute from init USER_INTERACTION)
    const res = await waitForEnergy(page, 12000);
    expect(res).toBeTruthy();
    if (res) {
        expect(typeof res.energy).toBe('number');
        // Positions are typically omitted for idle computes; assert non-presence when provided
        if (Array.isArray(res.positions)) expect(res.positions.length).toBeLessThanOrEqual(0);
    }

    // Ensure no protobuf init or WS decode errors
    const errBlob = consoleErrors.join('\n');
    expect(errBlob).not.toMatch(/protobuf-missing|Failed to parse message/);
});
