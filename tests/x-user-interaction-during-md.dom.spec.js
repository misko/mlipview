/**
 * Ported browser-style test: verifies USER_INTERACTION frames are sent while MD is running
 * and that drag lock messages still arrive.
 */

import { jest } from '@jest/globals';
import { stubWebSocketAndHook } from './utils/wsTestStub.js';
import { getWS } from '../public/fairchem_ws_client.js';

describe('x-user-interaction-during-md', () => {
  beforeEach(() => {
    if (!global.window) global.window = global;
  });

  test('dragging during MD emits USER_INTERACTION frames', async () => {
    window.__MLIPVIEW_SERVER = 'ws://127.0.0.1:8000';
    window.location = { protocol: 'http:', host: '127.0.0.1:4000' };
    const { sent, emit } = stubWebSocketAndHook();
    const ws = getWS();
    await ws.ensureConnected();

    const atoms = [1, 1];
    const initialPositions = [
      [0, 0, 0],
      [1, 0, 0],
    ];
    ws.initSystem({ atomic_numbers: atoms, positions: initialPositions });
    ws.startSimulation({
      type: 'md',
      params: { calculator: 'lj', temperature: 0, timestep_fs: 1, friction: 0 },
    });

    emit({
      seq: 1,
      positions: [
        [0.01, 0, 0],
        [1.01, 0, 0],
      ],
      simStep: 1,
    });

    const draggedIndex = 1;
    const dragPositions = [
      [0, 0, 0],
      [1.2, 0.1, 0],
    ];
    ws.userInteraction({ positions: dragPositions, dragLockIndex: draggedIndex });

    const firstUI = sent.find((msg) => msg && (msg.type === 'USER_INTERACTION' || msg.type === 1));
    expect(firstUI).toBeTruthy();
    expect(firstUI.positionsCount).toBe(2);

    emit({
      seq: 2,
      positions: [
        [0.02, 0, 0],
        [1.05, 0, 0],
      ],
      simStep: 2,
    });

    ws.userInteraction({ positions: dragPositions, dragLockIndex: draggedIndex });
    const uiCount = sent.filter((msg) => msg && (msg.type === 'USER_INTERACTION' || msg.type === 1))
      .length;
    expect(uiCount).toBeGreaterThanOrEqual(2);
  });
});
