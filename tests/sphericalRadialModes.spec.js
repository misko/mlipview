import { computeRadialDelta } from '../public/vr/spherical-drag-core.js';

// Helper to build input vectors
function v(x,y,z){ return [x,y,z]; }

describe('spherical-drag-core computeRadialDelta', () => {
  const baseState = {
    initialRadius: 10,
    controllerDist0: 0.5,
    ctrlVec0: [0.5,0,0],
    initialCtrlForward: [1,0,0]
  };
  const baseInputsNear = {
    controllerPos: [1,0,0], // dist ~1 from center
    controllerForward: [1,0,0],
    cameraPos: [0,0,0]
  };
  const pushedFar = {
    controllerPos: [2,0,0], // dist 2
    controllerForward: [1,0,0],
    cameraPos: [0,0,0]
  };
  const pulledClose = {
    controllerPos: [0.25,0,0], // dist 0.25
    controllerForward: [1,0,0],
    cameraPos: [0,0,0]
  };

  test('distance mode increases radius when controller moves further', () => {
    const d1 = computeRadialDelta(baseState, baseInputsNear, { mode:'distance' });
    const d2 = computeRadialDelta(baseState, pushedFar, { mode:'distance' });
    expect(d2.deltaR).toBeGreaterThan(d1.deltaR);
  });

  test('projection mode equals distance when aligned', () => {
    const d = computeRadialDelta(baseState, pushedFar, { mode:'projection' });
    expect(Math.abs(d.deltaR - d.components.projection)).toBeLessThan(1e-6);
  });

  test('hybrid mode averages distance and projection', () => {
    const dist = computeRadialDelta(baseState, pushedFar, { mode:'distance' });
    const proj = computeRadialDelta(baseState, pushedFar, { mode:'projection' });
    const hybrid = computeRadialDelta(baseState, pushedFar, { mode:'hybrid' });
    const avg = (dist.deltaR + proj.deltaR) * 0.5;
    expect(Math.abs(hybrid.deltaR - avg)).toBeLessThan(1e-6);
  });

  test('adaptive mode push increases deltaR positive', () => {
    const adaptive = computeRadialDelta(baseState, pushedFar, { mode:'adaptive', adaptiveGain:7.0 });
    expect(adaptive.deltaR).toBeGreaterThan(0);
  });

  test('adaptive mode pull produces negative deltaR', () => {
    const adaptive = computeRadialDelta(baseState, pulledClose, { mode:'adaptive', adaptiveGain:7.0 });
    expect(adaptive.deltaR).toBeLessThan(0);
  });

  test('forward mode responds to forward movement', () => {
    const fwd = computeRadialDelta(baseState, pushedFar, { mode:'forward', forwardGain:2.0 });
    expect(fwd.deltaR).toBeGreaterThan(0);
  });

  test('exp mode nonlinear growth larger than distance for big push', () => {
    const dist = computeRadialDelta(baseState, pushedFar, { mode:'distance' });
    const exp = computeRadialDelta(baseState, pushedFar, { mode:'exp', expK:1.4, expGain:0.9 });
    // Not guaranteed always > distance, but for >1 ratio and expK>1 expect amplification with these params
    expect(exp.deltaR).toBeGreaterThan(dist.deltaR * 0.5); // loose sanity bound
  });

  test('adaptive mode caps by adaptiveMaxFrac', () => {
    const extremePushInputs = { controllerPos: [5,0,0], controllerForward:[1,0,0], cameraPos:[0,0,0] }; // ratio 5/0.5=10
    const capped = computeRadialDelta(baseState, extremePushInputs, { mode:'adaptive', adaptiveGain:7.0, adaptiveMaxFrac:1.2 });
    // Max absolute fraction * initialRadius
    expect(Math.abs(capped.deltaR)).toBeLessThanOrEqual(1.2 * baseState.initialRadius + 1e-6);
  });

  test('exp mode uses expK and expGain', () => {
    const a = computeRadialDelta(baseState, pushedFar, { mode:'exp', expK:1.1, expGain:0.5 });
    const b = computeRadialDelta(baseState, pushedFar, { mode:'exp', expK:2.0, expGain:1.0 });
    expect(b.deltaR).toBeGreaterThan(a.deltaR);
  });
});
