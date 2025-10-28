// Lightweight event bus (typed via JSDoc for editor help)
export function createEventBus() {
  const listeners = new Map(); // event -> Set(fn)
  return {
    on(evt, fn) {
      if (!listeners.has(evt)) listeners.set(evt, new Set());
      listeners.get(evt).add(fn);
      return () => listeners.get(evt)?.delete(fn);
    },
    once(evt, fn) {
      const off = this.on(evt, (...args) => {
        off();
        fn(...args);
      });
      return off;
    },
    emit(evt, payload) {
      const set = listeners.get(evt);
      if (!set) return;
      for (const fn of Array.from(set)) {
        try {
          fn(payload);
        } catch (e) {
          console.warn('[eventBus] listener error', evt, e);
        }
      }
    },
    clear() {
      listeners.clear();
    },
  };
}
