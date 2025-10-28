import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

function makeViewer() {
  return {
    state: { dynamics: { energy: -2 }, bus: { on: () => {} } },
  };
}

describe('x-xr-hud-position', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('top panel hovers near center and bottom panel slightly below', () => {
    const viewer = makeViewer();
    window._viewer = viewer;
    const scene = makeHudScene();
    const hud = ensureWorldHUD({ scene, getViewer: () => viewer });
    runHudFrames(scene, 1);

    const halfHeight = Math.tan((scene.activeCamera.fov || Math.PI / 2) / 2) * 1.1;
    const topY = hud.planeTop.position.y;
    const bottomY = hud.plane.position.y;

    const topMin = halfHeight * 0.1;
    const topMax = halfHeight * 0.3;
    expect(topY).toBeGreaterThanOrEqual(topMin);
    expect(topY).toBeLessThanOrEqual(topMax);

    const bottomAbs = Math.abs(bottomY);
    const bottomMin = halfHeight * 0.35;
    const bottomMax = halfHeight * 0.55;
    expect(bottomAbs).toBeGreaterThanOrEqual(bottomMin);
    expect(bottomAbs).toBeLessThanOrEqual(bottomMax);
  });
});

