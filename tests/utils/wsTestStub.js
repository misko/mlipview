// Reusable WebSocket + WS test hook stub for jsdom tests.
export function stubWebSocketAndHook() {
  const sent = [];
  let wsOnMessage = null;
  global.WebSocket = class {
    constructor() {
      this.readyState = 1;
      setTimeout(() => this.onopen && this.onopen(), 0);
    }
    set binaryType(_) {}
    send(_) {}
    close() {}
    onopen() {}
    onerror() {}
    onclose() {}
    set onmessage(fn) {
      wsOnMessage = fn;
    }
    get onmessage() {
      return wsOnMessage;
    }
  };
  window.__WS_TEST_HOOK__ = (msg) => {
    sent.push(msg);
  };
  return {
    sent,
    emit: (obj) => {
      if (typeof window.__ON_WS_RESULT__ === 'function') window.__ON_WS_RESULT__(obj);
    },
  };
}
