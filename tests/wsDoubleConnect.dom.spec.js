/**
 * Mock-browser test: ensure we don't create two WS connections when ensureConnected
 * is called twice in parallel.
 */

import { getWS } from '../public/fairchem_ws_client.js';

describe('[mock-browser] WS double connect guard', () => {
  test('ensureConnected() called twice only creates one socket', async () => {
    // Spy on global WebSocket constructor
    const origWS = global.WebSocket;
    let constructed = 0;
    class FakeWS {
      constructor(url) {
        constructed++;
        this.url = url;
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen && this.onopen();
        }, 5);
      }
      send() {}
      close() {}
      onopen() {}
      onmessage() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;

    try {
      const ws = getWS();
      // Fire two connects in parallel
      await Promise.all([ws.ensureConnected(), ws.ensureConnected()]);
      expect(constructed).toBe(1);
      // Subsequent call should not re-create
      await ws.ensureConnected();
      expect(constructed).toBe(1);
    } finally {
      global.WebSocket = origWS;
    }
  });
});
