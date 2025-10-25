/**
 * Mock-browser regression guard: ensure ensureConnected() does not spawn
 * multiple WebSocket instances when called concurrently.
 */

import { getWS } from '../public/fairchem_ws_client.js';

describe('x-ws-double-connect', () => {
  test('parallel ensureConnected uses a single socket', async () => {
    const originalWS = global.WebSocket;
    let constructed = 0;

    class FakeWS {
      constructor(url) {
        constructed += 1;
        this.url = url;
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          if (typeof this.onopen === 'function') this.onopen();
        }, 5);
      }
      send() {}
      close() {}
    }

    global.WebSocket = FakeWS;

    try {
      const ws = getWS();
      await Promise.all([ws.ensureConnected(), ws.ensureConnected()]);
      expect(constructed).toBe(1);

      await ws.ensureConnected();
      expect(constructed).toBe(1);
    } finally {
      global.WebSocket = originalWS;
    }
  });
});

