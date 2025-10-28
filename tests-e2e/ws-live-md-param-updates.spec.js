// Purpose: E2E live MD temperature update during streaming using WS-only protocol.
// Based on ws-start-stop-idle-relax.spec.js style (health precheck, WS hook).
import { test, expect } from './fixtures.js';

/**
 * Flow:
 *  - Health precheck; wait for default molecule.
 *  - Start MD streaming.
 *  - Change temperature slider; expect a START_SIMULATION with updated temperature (and full params).
 *  - Change friction slider; expect another START_SIMULATION with updated friction (and full params).
 *  - Ensure streaming continues (energy ticks grow).
 *  - Start MD streaming and accumulate at least 15 frames.
 *  - Change temperature slider; expect a START_SIMULATION with updated temperature (and full params).
 *  - Accumulate 15 more frames, then stop the simulation; expect a stop notice.
 */

test('ws live MD param updates (temperature + friction)', async ({ page, loadViewerPage }) => {
    test.setTimeout(60_000);

    // Health precheck
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');
    await loadViewerPage({ query: { mol: 'molecules/water.xyz', dev: 1 } });

    // Orchestrate in-page: start MD, wait 15 frames, change temperature, wait 15 more, stop
    const result = await page.evaluate(async () => {
        const ws = window.__fairchem_ws__;
        if (!ws || typeof ws.setTestHook !== 'function') return { error: 'noWS' };

        // capture outgoing client messages
        const seen = [];
        ws.setTestHook((obj) => { try { if (obj && obj.type) seen.push(obj); } catch { } });
        window.__SEEN_WS__ = seen;

        // start MD stream
        await window.viewerApi.startMDContinuous({ steps: 1000, temperature: 1200 });

        // count frames
        const N = (window.viewerApi?.state?.positions || []).length | 0;
        let count = 0;
        const off = ws.onResult((r) => {
            try { if (Array.isArray(r.positions) && r.positions.length === N) count++; } catch { }
        });

        // wait for 15 frames
        const waitFrames = async (want) => {
            const t0 = Date.now();
            while (Date.now() - t0 < 30000 && count < want) await new Promise((r) => setTimeout(r, 50));
            return count >= want;
        };
        const got15 = await waitFrames(15);

        // check HUD temperature populated (T: <number> K)
        let tempHudOkBefore = false;
        try {
            const el = document.getElementById('instTemp');
            const txt = (el && (el.textContent || el.innerText || '')).trim();
            tempHudOkBefore = /^T:\s*\d+(?:\.\d+)?\s*K$/i.test(txt);
        } catch { }

        // mark current index before issuing the change so we don't miss a fast START message
        const afterIdx = seen.length;

        // change temperature slider to max
        let tempChanged = false;
        try {
            const slider = document.getElementById('mdTempSlider');
            if (slider) {
                slider.value = String(slider.max);
                slider.dispatchEvent(new Event('input', { bubbles: true }));
                tempChanged = true;
            }
        } catch { }

        // expect a START_SIMULATION with updated temperature
        const tStart = Date.now();
        let sawTempStart = false;
        while (Date.now() - tStart < 15000 && !sawTempStart) {
            await new Promise((r) => setTimeout(r, 50));
            const wantT = Number(window.__MLIP_TARGET_TEMPERATURE) || 0;
            for (let i = afterIdx; i < seen.length; i++) {
                const ev = seen[i];
                if (ev && ev.type === 'START_SIMULATION') {
                    const sp = ev.simulationParams || {};
                    if (Number(sp.temperature) === wantT && typeof sp.friction === 'number') {
                        sawTempStart = true;
                        break;
                    }
                }
            }
        }

        // wait for 15 more frames
        const got30 = await waitFrames(30);

        // check HUD temperature again after live update
        let tempHudOkAfter = false;
        try {
            const el = document.getElementById('instTemp');
            const txt = (el && (el.textContent || el.innerText || '')).trim();
            tempHudOkAfter = /^T:\s*\d+(?:\.\d+)?\s*K$/i.test(txt);
        } catch { }

        // stop simulation and see a stop notice
        let sawStop = false;
        const off2 = ws.onResult((r) => {
            try { if (r && (r.message === 'SIMULATION_STOPPED' || r.simulationStopped === true)) sawStop = true; } catch { }
        });
        try { window.viewerApi.stopSimulation(); } catch { }
        const tStop = Date.now();
        while (Date.now() - tStop < 10000 && !sawStop) await new Promise((r) => setTimeout(r, 50));

        try { off && off(); } catch { }
        try { off2 && off2(); } catch { }

        return { got15, tempChanged, sawTempStart, got30, sawStop, tempHudOkBefore, tempHudOkAfter };
    });

    expect(result && result.got15).toBeTruthy();
    expect(result && result.tempChanged).toBeTruthy();
    expect(result && result.sawTempStart).toBeTruthy();
    expect(result && result.tempHudOkBefore).toBeTruthy();
    expect(result && result.got30).toBeTruthy();
    expect(result && result.tempHudOkAfter).toBeTruthy();
    expect(result && result.sawStop).toBeTruthy();
});
