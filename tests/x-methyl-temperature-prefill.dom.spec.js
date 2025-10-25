/** @jest-environment jsdom */

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';
import { initTemperatureSlider } from '../public/ui/temperatureSlider.js';

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

function readFixture() {
  return `24
temp=500K (methyl radicals, reactive)
C    0.100   0.000   0.000
H    0.100   1.048   0.300
H    1.009  -0.524  -0.300
H   -0.809  -0.524   0.300
C    2.600   0.000   0.000
H    2.600   1.048  -0.300
H    3.509  -0.524   0.300
H    1.691  -0.524  -0.300
C    1.350   2.100   0.000
H    1.350   3.148   0.300
H    2.259   1.576  -0.300
H    0.441   1.576   0.300
C    1.350  -2.100   0.000
H    1.350  -1.052  -0.300
H    2.259  -2.624   0.300
H    0.441  -2.624  -0.300
C    1.350   0.000   2.100
H    1.350   1.048   2.400
H    2.259  -0.524   1.800
H    0.441  -0.524   2.400
C    1.350   0.000  -2.100
H    1.350   1.048  -2.400
H    2.259  -0.524  -1.800
H    0.441  -0.524  -2.400
`;
}

describe('x-methyl-temperature-prefill', () => {
  beforeEach(() => {
    document.body.innerHTML = '<div id="hud"></div>';
    delete window.__MLIP_TARGET_TEMPERATURE;
  });

  test('loading methyl radicals sets slider, state, and globals to 500K', () => {
    const viewer = makeViewer();
    const hud = document.getElementById('hud');
    const slider = initTemperatureSlider({ hudEl: hud, getViewer: () => viewer });

    const parsed = parseXYZ(readFixture());
    expect(parsed.temperature).toBe(500);
    applyParsedToViewer(viewer, parsed);

    const label = hud.querySelector('#tempLabel');
    const input = hud.querySelector('#mdTempSlider');
    expect(label && input).toBeTruthy();
    expect(label.textContent).toMatch(/500\s*K/);
    expect(slider.getTemperature()).toBe(500);
    expect(viewer.state.dynamics.targetTemperature).toBe(500);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(500);
  });
});

