// Timeline playback controller manages cadence, looping, and overrides.
// It remains agnostic of frame storage; the host supplies callbacks for stepping and state updates.

const DEFAULT_FPS = 20;

function clampFps(value) {
  const n = Number(value);
  if (!Number.isFinite(n) || n <= 0) return DEFAULT_FPS;
  return Math.max(1, Math.round(n));
}

function sanitizeStartFrame(ref = {}) {
  if (!ref || typeof ref !== 'object') return null;
  const out = {};
  if (typeof ref.frameId === 'string' && ref.frameId) out.frameId = ref.frameId;
  if (Number.isFinite(ref.offset)) out.offset = Math.floor(ref.offset);
  if (Number.isFinite(ref.frameIndex)) out.frameIndex = Math.max(0, Math.floor(ref.frameIndex));
  return Object.keys(out).length ? out : null;
}

function sanitizeLoopRange(range = {}) {
  if (!range || typeof range !== 'object') return null;
  const out = {};
  if (typeof range.startFrameId === 'string' && range.startFrameId) out.startFrameId = range.startFrameId;
  if (typeof range.endFrameId === 'string' && range.endFrameId) out.endFrameId = range.endFrameId;
  if (range.start && typeof range.start === 'object') {
    const start = {};
    if (typeof range.start.frameId === 'string' && range.start.frameId) start.frameId = range.start.frameId;
    if (Number.isFinite(range.start.offset)) start.offset = Math.floor(range.start.offset);
    if (Number.isFinite(range.start.frameIndex)) start.frameIndex = Math.max(0, Math.floor(range.start.frameIndex));
    if (Object.keys(start).length) out.start = start;
  }
  if (range.end && typeof range.end === 'object') {
    const end = {};
    if (typeof range.end.frameId === 'string' && range.end.frameId) end.frameId = range.end.frameId;
    if (Number.isFinite(range.end.offset)) end.offset = Math.floor(range.end.offset);
    if (Number.isFinite(range.end.frameIndex)) end.frameIndex = Math.max(0, Math.floor(range.end.frameIndex));
    if (typeof range.end.inclusive === 'boolean') end.inclusive = range.end.inclusive;
    if (Object.keys(end).length) out.end = end;
  }
  return Object.keys(out).length ? out : null;
}

export function createTimelinePlaybackController({
  defaultFps = DEFAULT_FPS,
  onStep = () => false,
  onPlaybackStateChange = () => {},
  getNow = () => (typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now()),
  schedule = (fn, ms) => setTimeout(fn, ms),
  cancel = (id) => clearTimeout(id),
} = {}) {
  let playbackStateListener = typeof onPlaybackStateChange === 'function' ? onPlaybackStateChange : () => {};
  const baseConfig = {
    defaultFps: clampFps(defaultFps),
    autoPlay: false,
    loop: false,
    loopRange: null,
    startFrame: null,
  };
  const overrides = {
    speed: null,
  };
  const state = {
    playing: false,
    timerId: null,
    lastTickAt: null,
  };

  function notifyPlaybackState() {
    try {
      playbackStateListener({
        playing: state.playing,
        effectiveFps: getEffectiveFps(),
      });
    } catch { /* noop */ }
  }

  function clearTimer() {
    if (state.timerId != null) {
      try { cancel(state.timerId); } catch { }
      state.timerId = null;
    }
  }

  function getEffectiveFps() {
    if (overrides.speed && overrides.speed.fps) {
      return clampFps(overrides.speed.fps);
    }
    return clampFps(baseConfig.defaultFps);
  }

  function tick() {
    state.timerId = null;
    state.lastTickAt = getNow();
    const keepPlaying = !!onStep({ fps: getEffectiveFps() });
    if (!state.playing) return;
    if (!keepPlaying) {
      stop();
      return;
    }
    scheduleNext();
  }

  function scheduleNext() {
    clearTimer();
    if (!state.playing) return;
    const fps = getEffectiveFps();
    const intervalMs = Math.max(1, Math.round(1000 / fps));
    state.timerId = schedule(tick, intervalMs);
  }

  function start({ immediate = false } = {}) {
    if (state.playing) {
      scheduleNext();
      return;
    }
    state.playing = true;
    state.lastTickAt = getNow();
    notifyPlaybackState();
    if (immediate) {
      tick();
      return;
    }
    scheduleNext();
  }

  function stop() {
    if (!state.playing) return false;
    clearTimer();
    state.playing = false;
    notifyPlaybackState();
    return true;
  }

  function setSpeedOverride(payload) {
    if (!payload || typeof payload !== 'object') {
      overrides.speed = null;
      scheduleNext();
      return;
    }
    if (payload.fps != null) {
      overrides.speed = {
        sourceId: payload.sourceId || payload.id || 'control-message',
        fps: clampFps(payload.fps),
        priority: payload.priority ?? 0,
      };
    } else if (payload.speedMultiplier != null && payload.speedMultiplier > 0) {
      overrides.speed = {
        sourceId: payload.sourceId || payload.id || 'control-message',
        fps: clampFps(baseConfig.defaultFps * Number(payload.speedMultiplier)),
        priority: payload.priority ?? 0,
      };
    } else {
      overrides.speed = null;
    }
    scheduleNext();
  }

  function clearOverrides() {
    overrides.speed = null;
    scheduleNext();
  }

  function setBaseConfig(cfg = {}) {
    if (!cfg || typeof cfg !== 'object') {
      return getSnapshot();
    }
    if (cfg.defaultFps != null) baseConfig.defaultFps = clampFps(cfg.defaultFps);
    if (typeof cfg.autoPlay === 'boolean') baseConfig.autoPlay = cfg.autoPlay;
    if (typeof cfg.loop === 'boolean') baseConfig.loop = cfg.loop;
    if (cfg.loopRange !== undefined) baseConfig.loopRange = sanitizeLoopRange(cfg.loopRange);
    if (cfg.startFrame !== undefined) baseConfig.startFrame = sanitizeStartFrame(cfg.startFrame);
    return getSnapshot();
  }

  function getSnapshot() {
    return {
      defaultFps: baseConfig.defaultFps,
      autoPlay: !!baseConfig.autoPlay,
      loop: !!baseConfig.loop,
      loopRange: baseConfig.loopRange ? { ...baseConfig.loopRange } : null,
      startFrame: baseConfig.startFrame ? { ...baseConfig.startFrame } : null,
    };
  }

  return {
    start,
    stop,
    isPlaying: () => state.playing,
    setSpeedOverride,
    clearOverrides,
    getEffectiveFps,
    setBaseConfig,
    applySnapshot: (cfg) => setBaseConfig(cfg),
    getSnapshot,
    shouldAutoPlay: () => !!baseConfig.autoPlay,
    getBaseConfig: () => getSnapshot(),
    getLoopConfig: () => ({
      loop: !!baseConfig.loop,
      loopRange: baseConfig.loopRange ? { ...baseConfig.loopRange } : null,
    }),
    getStartFrame: () => (baseConfig.startFrame ? { ...baseConfig.startFrame } : null),
    onPlaybackStateChange: (fn) => {
      if (typeof fn === 'function') {
        playbackStateListener = fn;
        notifyPlaybackState();
      }
    },
    updateTimer: () => scheduleNext(),
  };
}

export default { createTimelinePlaybackController };
