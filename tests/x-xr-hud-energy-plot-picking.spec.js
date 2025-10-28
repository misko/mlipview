import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

function makeViewer() {
  return {
    state: {
      dynamics: { energy: -3.2 },
      bus: { on: () => {} },
    },
  };
}

describe('x-xr-hud-energy-plot-picking', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('default layout keeps top plane non-pickable and farther away', () => {
    const viewer = makeViewer();
    window._viewer = viewer;
    const scene = makeHudScene();
    const hud = ensureWorldHUD({ scene, getViewer: () => viewer });
    expect(hud).toBeTruthy();
    runHudFrames(scene, 2);

    const bottomZ = hud.plane.position.z;
    const topZ = hud.planeTop.position.z;
    expect(bottomZ).toBeGreaterThan(1.0);
    expect(bottomZ).toBeLessThan(1.2);
    expect(topZ).toBeGreaterThan(2.0);
    expect(topZ).toBeLessThan(2.4);
    expect(hud.planeTop.isPickable).toBe(false);

    const widthRatio = hud.planeTop.opts.width / hud.plane.opts.width;
    expect(widthRatio).toBeGreaterThan(5.3);
    expect(widthRatio).toBeLessThan(5.5);
  });

  test('overrides can restore base distance and enable picking', () => {
    const viewer = makeViewer();
    window._viewer = viewer;

    const scene = makeHudScene();
    const hud = ensureWorldHUD({
      scene,
      getViewer: () => viewer,
      topEnergyDistanceMult: 1,
      topEnergyScaleMult: 1,
      topEnergyIsPickable: true,
    });
    runHudFrames(scene, 2);

    const dz = Math.abs(hud.planeTop.position.z - hud.plane.position.z);
    expect(dz).toBeLessThan(0.1);
    expect(hud.planeTop.isPickable).toBe(true);

    const widthRatio = hud.planeTop.opts.width / hud.plane.opts.width;
    expect(widthRatio).toBeGreaterThan(2.6);
    expect(widthRatio).toBeLessThan(2.8);
  });
});

