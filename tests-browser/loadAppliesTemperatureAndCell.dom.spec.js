/** @jest-environment jsdom */

// Verifies that applying parsed XYZ with temperature sets global and UI state,
// and that non-monoclinic cells trigger an error banner and are dropped.

import { applyParsedToViewer } from '../public/util/moleculeLoader.js';
import { buildCellFromParameters } from '../public/util/pbc.js';

function makeViewer() {
  const listeners = {};
  const bus = {
    on: (ev, fn) => {
      (listeners[ev] || (listeners[ev] = [])).push(fn);
    },
    emit: (ev) => {
      (listeners[ev] || []).forEach((fn) => fn());
    },
  };
  const state = {
    bus,
    elements: [],
    positions: [],
    bonds: [],
    cell: {
      a: { x: 1, y: 0, z: 0 },
      b: { x: 0, y: 1, z: 0 },
      c: { x: 0, y: 0, z: 1 },
      enabled: false,
      originOffset: { x: 0, y: 0, z: 0 },
    },
    showCell: true,
    markCellChanged() {
      bus.emit('cellChanged');
    },
    markPositionsChanged() {
      bus.emit('positionsChanged');
    },
    markBondsChanged() {
      bus.emit('bondsChanged');
    },
    dynamics: {},
  };
  const viewerApi = { state, recomputeBonds: () => {} };
  return viewerApi;
}

// Minimal error banner host
beforeEach(() => {
  document.body.innerHTML = '<div id="app"></div>';
});

describe('applyParsedToViewer temperature and cell handling', () => {
  test('sets global and state temperature', () => {
    const viewer = makeViewer();
    const parsed = { elements: ['H'], positions: [{ x: 0, y: 0, z: 0 }], temperature: 425 };
    applyParsedToViewer(viewer, parsed);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(425);
    expect(viewer.state.dynamics.targetTemperature).toBe(425);
  });

  test('accepts monoclinic cell params and drops non-monoclinic with error', () => {
    const viewer = makeViewer();
    const mono = buildCellFromParameters({ a: 10, b: 12, c: 15, alpha: 90, beta: 110, gamma: 90 });
    const parsed1 = { elements: ['H'], positions: [{ x: 0, y: 0, z: 0 }], cell: mono };
    applyParsedToViewer(viewer, parsed1);
    expect(viewer.state.cell && viewer.state.cell.enabled).toBe(true);
    // Non-monoclinic: alpha=100, beta=110, gamma=100 not convertible by permutation
    const tri = buildCellFromParameters({ a: 10, b: 11, c: 12, alpha: 100, beta: 110, gamma: 100 });
    const parsed2 = { elements: ['H'], positions: [{ x: 0, y: 0, z: 0 }], cell: tri };
    applyParsedToViewer(viewer, parsed2);
    // Should load but state.cell should be disabled or replaced; we drop cell (enabled false or missing vectors)
    expect(viewer.state.cell && viewer.state.cell.enabled).toBe(false);
    // Error banner element should be present
    const banner = document.getElementById('mlip-top-error-banner');
    expect(banner).toBeTruthy();
    expect(banner.textContent).toMatch(/not monoclinic/i);
  });
});
