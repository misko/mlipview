import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

function makeViewer() {
  const listener = () => {};
  return {
    state: { dynamics: { energy: -1 }, bus: { on: () => listener } },
  };
}

describe('x-xr-hud-energy-depth', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('topEnergyDistanceMult pushes top panel forward without large Y drift', () => {
    const viewer = makeViewer();
    window._viewer = viewer;

    const scene1 = makeHudScene();
    const baseHud = ensureWorldHUD({ scene: scene1, getViewer: () => viewer, topEnergyDistanceMult: 1 });
    runHudFrames(scene1, 2);

    resetXRHudSingleton();

    const scene2 = makeHudScene();
    const farHud = ensureWorldHUD({ scene: scene2, getViewer: () => viewer, topEnergyDistanceMult: 3 });
    runHudFrames(scene2, 2);

    expect(farHud.planeTop.position.z).toBeGreaterThan(baseHud.planeTop.position.z * 2.5);
    const baseY = Math.abs(baseHud.planeTop.position.y);
    const farY = Math.abs(farHud.planeTop.position.y);
    expect(farY).toBeLessThan(baseY * 1.2 + 1e-6);
  });
});

