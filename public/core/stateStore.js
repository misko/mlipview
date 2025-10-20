// Centralized shared state for WS-only frontend
// Manages user interaction counter, dragging/latched atoms, and frame acceptance logic

let __singleton = null;

export function createSharedStateStore() {
  const dragging = new Set();
  const latchedUntil = new Map(); // atomIndex -> expiresAt (ms since performance.now or Date.now)
  let userInteractionCount = 0;
  let totalInteractionCount = 0;

  function nowMs() {
    try {
      return typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now();
    } catch {
      return Date.now();
    }
  }

  return {
    // Counters
    getUserInteractionCount() {
      return userInteractionCount | 0;
    },
    getTotalInteractionCount() {
      return totalInteractionCount | 0;
    },
    setUserInteractionCount(n) {
      if (Number.isFinite(n)) userInteractionCount = n | 0;
      return userInteractionCount;
    },
    setTotalInteractionCount(n) {
      if (Number.isFinite(n)) totalInteractionCount = n | 0;
      return totalInteractionCount;
    },
    bumpUserInteraction(reason) {
      userInteractionCount = (userInteractionCount | 0) + 1;
      totalInteractionCount = (totalInteractionCount | 0) + 1;
      return { userInteractionCount, totalInteractionCount, reason };
    },
    bumpSimulation(reason) {
      totalInteractionCount = (totalInteractionCount | 0) + 1;
      return { totalInteractionCount, reason };
    },

    // Dragging set
    addDragging(i) {
      if (Number.isInteger(i) && i >= 0) dragging.add(i);
      return dragging;
    },
    removeDragging(i) {
      dragging.delete(i);
      return dragging;
    },
    clearDragging() {
      dragging.clear();
      return dragging;
    },
    getDraggingSet() {
      return dragging;
    },

    // Latch helpers (to avoid immediate snapback right after release)
    latchAtom(i, ms = 600) {
      if (!Number.isInteger(i) || i < 0) return;
      const t = nowMs() + Math.max(0, Number(ms) || 0);
      latchedUntil.set(i, t);
    },
    collectExcludeSet() {
      const out = new Set();
      for (const i of dragging) out.add(i);
      const tnow = nowMs();
      for (const [i, t] of latchedUntil) {
        if (t > tnow) out.add(i);
        else latchedUntil.delete(i);
      }
      return out;
    },

    // Acceptance logic for incoming frames
    shouldAcceptFrame(respUserInteractionCount) {
      // Accept if server didn't include the field
      if (respUserInteractionCount == null) return true;
      const cur = userInteractionCount | 0;
      // Discard if server's count is behind current front-end interaction counter
      return (respUserInteractionCount | 0) >= cur;
    },

    // Testing/reset utilities
    __reset() {
      dragging.clear();
      latchedUntil.clear();
      userInteractionCount = 0;
      totalInteractionCount = 0;
    },
  };
}

export function getSharedStateStore() {
  if (!__singleton) __singleton = createSharedStateStore();
  return __singleton;
}

export default { getSharedStateStore, createSharedStateStore };
