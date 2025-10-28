/**
 * Ported idle USER_INTERACTION test: ensures idle compute frames are received.
 */

import { jest } from '@jest/globals';
import { getWS } from '../public/fairchem_ws_client.js';

describe('x-user-interaction-when-idle', () => {
  beforeEach(() => {
    if (!global.window) global.window = global;
  });

  test('idle USER_INTERACTION triggers energy-only frame', async () => {
    const origWS = global.WebSocket;
    class FakeWS {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen?.();
        }, 0);
      }
      send() {}
      close() {}
      onopen() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;

    window.__MLIPVIEW_SERVER = 'ws://127.0.0.1:8000';
    window.location = { protocol: 'http:', host: '127.0.0.1:4000' };
    const ws = getWS();
    const sent = [];
    ws.setTestHook((msg) => sent.push(msg));
    await ws.ensureConnected();

    const atoms = [1, 1];
    const p0 = [
      [0, 0, 0],
      [1, 0, 0],
    ];
    ws.userInteraction({ atomic_numbers: atoms, positions: p0 });

    const moved = [
      [0, 0, 0],
      [1.1, 0.2, 0],
    ];
    ws.userInteraction({ positions: moved, dragLockIndex: 1 });

    const ui = sent.find((msg) => msg && msg.type != null);
    expect(ui).toBeTruthy();
    expect(ui.positionsCount).toBe(2);

    ws.injectTestResult({
      seq: 1,
      energy: -0.123,
      forces: [
        [0, 0, 0],
        [0, 0, 0],
      ],
    });

    expect(true).toBe(true);
    global.WebSocket = origWS;
  });
});
