import { computeRadialDelta } from '../public/vr/spherical-drag-core.js';

const baseState = {
  initialRadius: 10,
  controllerDist0: 0.5,
  ctrlVec0: [0.5, 0, 0],
  initialCtrlForward: [1, 0, 0],
};

const nearInputs = {
  controllerPos: [1, 0, 0],
  controllerForward: [1, 0, 0],
  cameraPos: [0, 0, 0],
};

const pushedFar = {
  controllerPos: [2, 0, 0],
  controllerForward: [1, 0, 0],
  cameraPos: [0, 0, 0],
};

const pulledClose = {
  controllerPos: [0.25, 0, 0],
  controllerForward: [1, 0, 0],
  cameraPos: [0, 0, 0],
};

describe('x-spherical-radial-modes', () => {
  test('distance mode grows deltaR when controller moves further', () => {
    const a = computeRadialDelta(baseState, nearInputs, { mode: 'distance' });
    const b = computeRadialDelta(baseState, pushedFar, { mode: 'distance' });
    expect(b.deltaR).toBeGreaterThan(a.deltaR);
  });

  test('projection mode equals projection component when aligned', () => {
    const res = computeRadialDelta(baseState, pushedFar, { mode: 'projection' });
    expect(Math.abs(res.deltaR - res.components.projection)).toBeLessThan(1e-6);
  });

  test('hybrid mode averages distance and projection', () => {
    const dist = computeRadialDelta(baseState, pushedFar, { mode: 'distance' });
    const proj = computeRadialDelta(baseState, pushedFar, { mode: 'projection' });
    const hybrid = computeRadialDelta(baseState, pushedFar, { mode: 'hybrid' });
    const expected = (dist.deltaR + proj.deltaR) * 0.5;
    expect(Math.abs(hybrid.deltaR - expected)).toBeLessThan(1e-6);
  });

  test('adaptive mode push yields positive deltaR', () => {
    const adaptive = computeRadialDelta(baseState, pushedFar, {
      mode: 'adaptive',
      adaptiveGain: 7,
    });
    expect(adaptive.deltaR).toBeGreaterThan(0);
  });

  test('adaptive mode pull yields negative deltaR', () => {
    const adaptive = computeRadialDelta(baseState, pulledClose, {
      mode: 'adaptive',
      adaptiveGain: 7,
    });
    expect(adaptive.deltaR).toBeLessThan(0);
  });

  test('forward mode responds to controller forward motion', () => {
    const res = computeRadialDelta(baseState, pushedFar, { mode: 'forward', forwardGain: 2 });
    expect(res.deltaR).toBeGreaterThan(0);
  });

  test('exp mode grows faster for larger expK/gain', () => {
    const a = computeRadialDelta(baseState, pushedFar, { mode: 'exp', expK: 1.1, expGain: 0.5 });
    const b = computeRadialDelta(baseState, pushedFar, { mode: 'exp', expK: 2.0, expGain: 1.0 });
    expect(b.deltaR).toBeGreaterThan(a.deltaR);
  });

  test('adaptive mode respects max fraction clamp', () => {
    const extreme = {
      controllerPos: [5, 0, 0],
      controllerForward: [1, 0, 0],
      cameraPos: [0, 0, 0],
    };
    const capped = computeRadialDelta(baseState, extreme, {
      mode: 'adaptive',
      adaptiveGain: 7,
      adaptiveMaxFrac: 1.2,
    });
    expect(Math.abs(capped.deltaR)).toBeLessThanOrEqual(baseState.initialRadius * 1.2 + 1e-6);
  });
});
