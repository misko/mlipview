/** @jest-environment jsdom */

import {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
  findGuiNode,
} from './utils/xrHudTestHarness.js';
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

async function waitForText(node, substring, attempts = 20) {
  for (let i = 0; i < attempts; i += 1) {
    if (node?.text?.includes(substring)) return;
    await new Promise((resolve) => setTimeout(resolve, 10));
  }
  throw new Error('Timed out waiting for text containing ' + substring);
}

describe('x-xr-energy-hud', () => {
  beforeEach(() => {
    setupXRHudTestEnv();
    resetXRHudSingleton();
  });

  test('energy panel builds and responds to forcesChanged events', async () => {
    const listeners = new Map();
    const bus = {
      on(evt, fn) {
        const arr = listeners.get(evt) || [];
        arr.push(fn);
        listeners.set(evt, arr);
      },
      emit(evt) {
        const arr = listeners.get(evt) || [];
        for (const fn of arr) {
          try {
            fn();
          } catch {}
        }
      },
    };

    const viewer = {
      state: {
        dynamics: { energy: -12.345 },
        bus,
      },
    };

    window.viewerApi = viewer;
    window._viewer = viewer;

    const scene = makeHudScene();
    const hud = ensureWorldHUD({ scene, getViewer: () => viewer });

    expect(hud).toBeTruthy();
    expect(window.__XR_HUD_FALLBACK).toBe(hud);
    expect(hud.energy).toBe(true);
    expect(hud.buttons).toEqual(expect.arrayContaining(['Relax', 'MD', 'Off', 'Forces', 'Reset']));

    runHudFrames(scene, 2);

    const root = hud.texTop._rootContainer;
    const energyImg = findGuiNode(root, (node) => node.id === 'xrWorldEnergyImg_top');
    const energyText = findGuiNode(root, (node) => node.id === 'xrWorldEnergyValue_top');

    expect(energyImg).toBeTruthy();
    expect(energyText).toBeTruthy();
    await waitForText(energyText, '-12.345');

    viewer.state.dynamics.energy = -10.001;
    bus.emit('forcesChanged');
    await waitForText(energyText, '-10.001');
  });
});
