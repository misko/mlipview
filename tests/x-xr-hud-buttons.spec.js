/**
 * Validates XR HUD button styling and highlight behaviour using shared HUD harness.
 */

import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
  findGuiNode,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

describe('x-xr-hud-buttons', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('forces and PBC buttons highlight when toggled', async () => {
    const scene = makeHudScene();
    const viewer = {
      state: {
        showForces: false,
        showCell: false,
        cell: { enabled: false },
        toggleForceVectorsVisibility() {
          this.showForces = !this.showForces;
        },
        toggleCellVisibility() {
          this.showCell = !this.showCell;
          this.cell.enabled = this.showCell;
        },
        toggleGhostCells: () => {},
      },
    };
    window._viewer = viewer;

    const hud = ensureWorldHUD({ scene, getViewer: () => viewer });
    expect(hud).toBeTruthy();

    expect(hud.plane.opts.height).toBeGreaterThan(0.25);
    expect(hud.plane.opts.height).toBeLessThan(0.28);

    const root = window.__XR_HUD_FALLBACK?.tex?._rootContainer;
    expect(root).toBeTruthy();

    const forcesBtn = findGuiNode(root, (node) => node?.text === 'Forces' || node?.textBlock?.text === 'Forces');
    const pbcBtn = findGuiNode(root, (node) => node?.text === 'PBC' || node?.textBlock?.text === 'PBC');
    expect(forcesBtn).toBeTruthy();
    expect(pbcBtn).toBeTruthy();

    const cm = await (async () => {
      for (let i = 0; i < 12; i += 1) {
        const model = hud.getControlsModel?.();
        if (model) return model;
        await new Promise((resolve) => setTimeout(resolve, 0));
      }
      return null;
    })();
    expect(cm).toBeTruthy();

    const beforeForces = forcesBtn.background;
    cm.toggleForces();
    await new Promise((resolve) => setTimeout(resolve, 0));
    runHudFrames(scene, 2);
    const afterForces = forcesBtn.background;
    expect(afterForces).not.toBe(beforeForces);
    expect(afterForces).toBe('rgba(40,140,100,0.9)');

    const beforePbc = pbcBtn.background;
    cm.togglePBC();
    await new Promise((resolve) => setTimeout(resolve, 0));
    runHudFrames(scene, 2);
    const afterPbc = pbcBtn.background;
    expect(afterPbc).not.toBe(beforePbc);
    expect(afterPbc).toBe('rgba(40,140,100,0.9)');
  });
});

