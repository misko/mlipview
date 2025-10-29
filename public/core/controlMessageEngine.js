// Control message engine: resolves timeline control messages to per-frame actions.

function clone(obj) {
  if (obj == null) return obj;
  try {
    return JSON.parse(JSON.stringify(obj));
  } catch {
    return obj;
  }
}

function ensureArray(val) {
  if (Array.isArray(val)) return val;
  return [];
}

function sanitizeAction(action) {
  if (!action || typeof action !== 'object') return null;
  const type = action.type;
  if (typeof type !== 'string') return null;
  switch (type) {
    case 'timeline.playbackSpeed': {
      const payload = {};
      if (Number.isFinite(action.fps)) payload.fps = Number(action.fps);
      if (Number.isFinite(action.speedMultiplier) && action.speedMultiplier > 0) {
        payload.speedMultiplier = Number(action.speedMultiplier);
      }
      if (payload.fps == null && payload.speedMultiplier == null) return null;
      if (Number.isFinite(action.transitionMs)) payload.transitionMs = Math.max(0, Math.floor(action.transitionMs));
      if (typeof action.easing === 'string' && action.easing) payload.easing = action.easing;
      return { type, ...payload };
    }
    case 'overlay.callout': {
      const payload = {
        text: typeof action.text === 'string' ? action.text : '',
        anchor: action.anchor && typeof action.anchor === 'object' ? clone(action.anchor) : null,
      };
      if (!payload.anchor) return null;
      if (Number.isFinite(action.textSize)) payload.textSize = Number(action.textSize);
      if (action.panelSize && typeof action.panelSize === 'object') {
        const { width, height } = action.panelSize;
        payload.panelSize = {};
        if (Number.isFinite(width)) payload.panelSize.width = Number(width);
        if (Number.isFinite(height)) payload.panelSize.height = Number(height);
      }
      if (action.style && typeof action.style === 'object') {
        payload.style = clone(action.style);
      }
      if (action.offset && Array.isArray(action.offset)) payload.offset = action.offset.map((v) => Number(v) || 0);
      return { type, ...payload };
    }
    case 'visual.opacityFocus': {
      const payload = {};
      if (action.focus && typeof action.focus === 'object') {
        payload.focus = clone(action.focus);
      }
      if (!payload.focus) return null;
      if (action.focusOpacity && typeof action.focusOpacity === 'object') {
        payload.focusOpacity = clone(action.focusOpacity);
      }
      if (action.backgroundOpacity && typeof action.backgroundOpacity === 'object') {
        payload.backgroundOpacity = clone(action.backgroundOpacity);
      }
      if (Number.isFinite(action.transitionMs)) payload.transitionMs = Math.max(0, Math.floor(action.transitionMs));
      return { type, ...payload };
    }
    default:
      return null;
  }
}

export function createControlMessageEngine(resolvers = {}) {
  const resolver = {
    resolveFrameIndexById: resolvers.resolveFrameIndexById || (() => null),
    offsetToIndex: resolvers.offsetToIndex || (() => null),
    getFrameCount: resolvers.getFrameCount || (() => 0),
  };

  let rawMessages = [];
  let resolvedMessages = [];
  let frameCount = 0;

  function setResolvers(next = {}) {
    if (typeof next.resolveFrameIndexById === 'function') resolver.resolveFrameIndexById = next.resolveFrameIndexById;
    if (typeof next.offsetToIndex === 'function') resolver.offsetToIndex = next.offsetToIndex;
    if (typeof next.getFrameCount === 'function') resolver.getFrameCount = next.getFrameCount;
    refreshResolved();
  }

  function updateFrameCount(count) {
    frameCount = Math.max(0, Number(count) || 0);
    refreshResolved();
  }

  function resolveFrameIndexFromRef(ref) {
    if (!ref || typeof ref !== 'object') return null;
    if (typeof ref.frameId === 'string' && ref.frameId) {
      const idxById = resolver.resolveFrameIndexById(ref.frameId);
      if (Number.isInteger(idxById)) return clampIndex(idxById);
    }
    if (Number.isFinite(ref.frameIndex)) {
      return clampIndex(Math.floor(ref.frameIndex));
    }
    if (Number.isFinite(ref.offset)) {
      const idx = resolver.offsetToIndex(ref.offset);
      if (Number.isInteger(idx)) return clampIndex(idx);
    }
    return null;
  }

  function clampIndex(idx) {
    if (!Number.isInteger(idx)) return null;
    if (frameCount <= 0) return null;
    if (idx < 0) return 0;
    if (idx >= frameCount) return frameCount - 1;
    return idx;
  }

  function resolveRange(range = {}) {
    if (frameCount <= 0) return null;
    const resolved = {};
    const startRef = range.start || {};
    const endRef = range.end || {};
    let startIndex = resolveFrameIndexFromRef(startRef);
    let endIndex = resolveFrameIndexFromRef(endRef);

    if (startIndex == null) startIndex = 0;
    if (endIndex == null) endIndex = frameCount - 1;
    if (endRef && typeof endRef.inclusive === 'boolean' && !endRef.inclusive) {
      endIndex -= 1;
    }
    startIndex = clampIndex(startIndex);
    endIndex = clampIndex(endIndex);
    if (startIndex == null || endIndex == null) return null;
    if (startIndex > endIndex) return null;
    resolved.startIndex = startIndex;
    resolved.endIndex = endIndex;
    return resolved;
  }

  function refreshResolved() {
    frameCount = Math.max(0, resolver.getFrameCount() || frameCount || 0);
    resolvedMessages = [];
    if (!frameCount || !rawMessages.length) return;
    rawMessages.forEach((msg, idx) => {
      const range = resolveRange(msg.range || {});
      if (!range) return;
      const actions = ensureArray(msg.actions)
        .map((action) => sanitizeAction(action))
        .filter(Boolean);
      if (!actions.length) return;
      resolvedMessages.push({
        id: msg.id || `control-${idx}`,
        priority: Number(msg.priority) || 0,
        range,
        actions,
        order: idx,
        raw: msg,
      });
    });
  }

  function setMessages(messages = []) {
    rawMessages = ensureArray(messages)
      .map((msg, idx) => {
        if (!msg || typeof msg !== 'object') return null;
        const actions = ensureArray(msg.actions).map((action) => sanitizeAction(action)).filter(Boolean);
        if (!actions.length) return null;
        const range = msg.range && typeof msg.range === 'object' ? clone(msg.range) : {};
        const entry = {
          id: typeof msg.id === 'string' && msg.id ? msg.id : `control-${idx}`,
          priority: Number(msg.priority) || 0,
          range,
          actions,
        };
        return entry;
      })
      .filter(Boolean);
    refreshResolved();
  }

  function evaluate(frameIndex) {
    if (!Number.isInteger(frameIndex) || frameCount <= 0 || !resolvedMessages.length) {
      return {
        speed: null,
        callout: null,
        opacity: null,
        activeMessages: [],
      };
    }
    const active = resolvedMessages.filter(
      (msg) => frameIndex >= msg.range.startIndex && frameIndex <= msg.range.endIndex,
    );
    if (!active.length) {
      return {
        speed: null,
        callout: null,
        opacity: null,
        activeMessages: [],
      };
    }
    const sorted = active
      .slice()
      .sort((a, b) => {
        const prio = (b.priority || 0) - (a.priority || 0);
        if (prio !== 0) return prio;
        return (a.order || 0) - (b.order || 0);
      });
    const result = {
      speed: null,
      callout: null,
      opacity: null,
      activeMessages: sorted.map((msg) => ({
        id: msg.id,
        priority: msg.priority,
        startIndex: msg.range.startIndex,
        endIndex: msg.range.endIndex,
      })),
    };
    for (const msg of sorted) {
      for (const action of msg.actions) {
        if (action.type === 'timeline.playbackSpeed' && !result.speed) {
          result.speed = { ...action, sourceId: msg.id, priority: msg.priority };
        } else if (action.type === 'overlay.callout' && !result.callout) {
          result.callout = { ...action, sourceId: msg.id, priority: msg.priority };
        } else if (action.type === 'visual.opacityFocus' && !result.opacity) {
          result.opacity = { ...action, sourceId: msg.id, priority: msg.priority };
        }
      }
    }
    return result;
  }

  return {
    setMessages,
    applySnapshot: (messages) => setMessages(messages),
    getSnapshot: () => rawMessages.map((msg) => ({
      id: msg.id,
      priority: msg.priority,
      range: clone(msg.range),
      actions: msg.actions.map((a) => clone(a)),
    })),
    evaluate,
    setResolvers,
    updateFrameCount,
    refresh: refreshResolved,
    clear: () => {
      rawMessages = [];
      resolvedMessages = [];
    },
  };
}

export default { createControlMessageEngine };
