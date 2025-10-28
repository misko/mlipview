import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

function makeViewer() {
  return {
    state: { dynamics: { energy: -4 }, bus: { on: () => {} } },
  };
}

describe('x-xr-hud-energy-scale', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('topEnergyScaleMult scales width and height proportionally', () => {
    const viewer = makeViewer();
    window._viewer = viewer;

    const scene1 = makeHudScene();
    const hud1 = ensureWorldHUD({
      scene: scene1,
      getViewer: () => viewer,
      topEnergyScaleMult: 1,
    });
    runHudFrames(scene1, 2);

    resetXRHudSingleton();

    const scene2 = makeHudScene();
    const hud2 = ensureWorldHUD({
      scene: scene2,
      getViewer: () => viewer,
      topEnergyScaleMult: 2,
    });
    runHudFrames(scene2, 2);

    const widthRatio = hud2.planeTop.opts.width / hud1.planeTop.opts.width;
    const heightRatio = hud2.planeTop.opts.height / hud1.planeTop.opts.height;
    expect(widthRatio).toBeGreaterThan(1.9);
    expect(widthRatio).toBeLessThan(2.1);
    expect(heightRatio).toBeGreaterThan(1.9);
    expect(heightRatio).toBeLessThan(2.1);
  });
});

