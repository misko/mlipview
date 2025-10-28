// Purpose: End-to-end: when a cell is present, idle frames include energy and may include stress.
// Note: Server includes stress only if calculator provides it; this test asserts presence OR graceful absence.

import { test, expect } from './fixtures.js';

test.describe('WS protocol: cell and optional stress', () => {
  test('idle compute with cell includes energy and may include stress', async ({
    page,
    loadViewerPage,
  }) => {
    test.setTimeout(45_000);
    await loadViewerPage({ query: { autoMD: 0, wsDebug: 1 } });

    const result = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      const api = window.viewerApi;
      await ws.ensureConnected();
      const cell = [
        [5.0, 0, 0],
        [0, 5.0, 0],
        [0, 0, 5.0],
      ];
      // Correlate with a distinctive userInteractionCount to avoid stray frames
      const cur = 222333444;
      ws.setCounters({ userInteractionCount: cur });
      // Send positions and cell to trigger idle compute
      const pos = api.state.positions.map((p) => [p.x, p.y, p.z]);
      return await new Promise((resolve) => {
        const off = ws.onResult((r) => {
          try {
            if (r && typeof r.energy === 'number' && (r.userInteractionCount | 0) === cur) {
              const hasPositions = Array.isArray(r.positions) && r.positions.length > 0;
              const out = { hasStress: Array.isArray(r.stress), hasPositions, energy: r.energy };
              try {
                off && off();
              } catch { }
              resolve(out);
            }
          } catch { }
        });
        ws.userInteraction({ positions: pos, cell });
        setTimeout(() => {
          try {
            off && off();
          } catch { }
          resolve(null);
        }, 30000);
      });
    });

    expect(result).toBeTruthy();
    expect(result.energy).toBeDefined();
    if (result.hasPositions) {
      test
        .info()
        .annotations.push({
          type: 'note',
          description: 'Idle frame included positions (expected under protobuf migration).',
        });
    }
    // Allow either presence or absence depending on calculator support; log when absent to encourage server support
    if (!result.hasStress) {
      test
        .info()
        .annotations.push({
          type: 'note',
          description:
            'No stress in idle frame; ensure backend calculator populates stress when available.',
        });
    }
  });
});
