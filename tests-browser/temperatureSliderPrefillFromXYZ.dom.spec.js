/** @jest-environment jsdom */

// Verifies that when a temperature is set before building the desktop panel/slider,
// the UI initializes the slider to the nearest tick and label shows that temperature.

import { applyParsedToViewer } from '../public/util/moleculeLoader.js';

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
  const viewerApi = { state, recomputeBonds: () => {} };
  return viewerApi;
}

beforeEach(() => {
  document.body.innerHTML = '<div id="app"></div>';
  delete window.__MLIP_TARGET_TEMPERATURE;
});

test('temperature from XYZ pre-fills slider on build', async () => {
  const viewer = makeViewer();
  // Simulate XYZ parsed temperature
  const parsed = { elements: ['H'], positions: [{ x: 0, y: 0, z: 0 }], temperature: 425 };
  applyParsedToViewer(viewer, parsed);
  expect(window.__MLIP_TARGET_TEMPERATURE).toBe(425);
  // Now build the desktop panel which installs the slider
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const panel = buildDesktopPanel({ attachTo: document.getElementById('app') });
  // Slider should exist and label reflect nearest tick to 425
  const slider = document.getElementById('mdTempSlider');
  const label = document.getElementById('tempLabel');
  expect(slider).toBeTruthy();
  expect(label.textContent).toMatch(/T=\d+K/);
  const shown = Number(label.textContent.replace(/[^0-9]/g, ''));
  expect(Math.abs(shown - 425)).toBeLessThanOrEqual(60); // nearest of 30 ticks across 0..3000 is within ~60K
});
