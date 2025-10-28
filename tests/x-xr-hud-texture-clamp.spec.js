import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

function makeViewer() {
  return {
    state: { dynamics: { energy: -5 }, bus: { on: () => {} } },
  };
}

describe('x-xr-hud-texture-clamp', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('clamps texture dimensions when scaling exceeds maxHUDTextureDim', () => {
    const viewer = makeViewer();
    window._viewer = viewer;
    const scene = makeHudScene();
    const hud = ensureWorldHUD({
      scene,
      getViewer: () => viewer,
      topEnergyScaleMult: 8,
      maxHUDTextureDim: 1024,
    });
    runHudFrames(scene, 1);

    const md = hud.planeTop.metadata?.hudTexSize;
    expect(md).toBeTruthy();
    expect(md.rawW).toBeGreaterThan(1024);
    expect(md.rawH).toBeGreaterThan(1024);
    expect(md.clamped).toBe(true);
    expect(md.w).toBeLessThanOrEqual(1024);
    expect(md.h).toBeLessThanOrEqual(1024);
    expect(md.max).toBe(1024);
  });
});

