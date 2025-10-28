/** @jest-environment jsdom */

import { applyParsedToViewer } from '../public/util/moleculeLoader.js';

function makeViewer() {
  const listeners = {};
  const bus = {
    on(ev, fn) {
      (listeners[ev] || (listeners[ev] = [])).push(fn);
    },
    emit(ev) {
      for (const fn of listeners[ev] || []) fn();
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
    showCell: false,
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
  return { state, recomputeBonds: () => {} };
}

describe('x-temperature-slider-prefill', () => {
  beforeEach(() => {
    document.body.innerHTML = '<div id="app"></div>';
    delete window.__MLIP_TARGET_TEMPERATURE;
  });

  test('desktop panel slider pre-fills from previously applied XYZ temperature', async () => {
    const viewer = makeViewer();
    applyParsedToViewer(viewer, {
      elements: ['H'],
      positions: [{ x: 0, y: 0, z: 0 }],
      temperature: 425,
    });
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(425);

    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app'), getViewer: () => viewer });

    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    expect(slider).toBeTruthy();
    expect(label).toBeTruthy();
    const shown = Number(label.textContent.replace(/[^0-9]/g, ''));
    expect(Math.abs(shown - 425)).toBeLessThanOrEqual(60);
  });
});

