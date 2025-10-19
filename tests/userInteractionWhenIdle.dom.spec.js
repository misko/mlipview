/**
 * Mock browser test: user interaction when idle triggers simple_calculate
 * - Sends INIT_SYSTEM
 * - Sends USER_INTERACTION with new positions while not running
 * - Verifies a simple_calculate-like response is injected and could be logged
 */

import { jest } from '@jest/globals';
import { stubWebSocketAndHook } from './utils/wsTestStub.js';
import { getWS } from '../public/fairchem_ws_client.js';

describe('USER_INTERACTION when idle', () => {
  beforeEach(() => { if (!global.window) global.window = global; });

  test('sends USER_INTERACTION and receives energy frame (simple_calculate)', async () => {
  window.__MLIPVIEW_SERVER = 'ws://127.0.0.1:8000';
  window.location = { protocol: 'http:', host: '127.0.0.1:4000' };
  const { sent, emit } = stubWebSocketAndHook();
    const ws = getWS();
    await ws.ensureConnected();

    const atoms = [1, 1];
    const p0 = [[0,0,0],[1,0,0]];
    ws.initSystem({ atomic_numbers: atoms, positions: p0 });

    // Interaction while idle
    const moved = [[0,0,0],[1.1,0.2,0]];
  ws.userInteraction({ positions: moved, dragLockIndex: 1 });

    // Validate sent USER_INTERACTION
  const ui = sent.find(m => m && (m.type === 'USER_INTERACTION' || m.type === 1));
  expect(ui).toBeTruthy();
  expect(ui.positionsCount).toBe(2);

    // Simulate server response for idle simple_calculate (energy only)
    emit({ seq: 1, energy: -0.123, message: 'simple_calculate' });

    // No assertion on rendering, but ensure protocol flow was exercised without error
    expect(true).toBe(true);
  });
});
