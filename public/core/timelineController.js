// Timeline controller — orchestrates playback vs live state.
// Depends on a timelineStore (ring buffer) and host-provided mutators.

function getRAF() {
  if (typeof requestAnimationFrame === 'function') {
    return requestAnimationFrame.bind(globalThis || window);
  }
  return (fn) => setTimeout(() => fn(Date.now()), 16);
}

function getCAF() {
  if (typeof cancelAnimationFrame === 'function') {
    return cancelAnimationFrame.bind(globalThis || window);
  }
  return (token) => clearTimeout(token);
}

export function createTimelineController({
  store,
  applyFrame,
  setInteractionsEnabled,
} = {}) {
  if (!store) throw new Error('timelineController requires a store');
  if (typeof applyFrame !== 'function') {
    throw new Error('timelineController requires applyFrame(frame, meta)');
  }

  const raf = getRAF();
  const caf = getCAF();

  const listeners = new Map(); // type -> Set
  const resumeMeta = {
    lastRunningKind: null,
    lastStartArgs: null,
  };

  let bufferSnapshot = store.getSnapshot();
  let mode = 'live';
  let playbackIndex = bufferSnapshot.liveIndex;
  let paused = false;

  let playbackStartIndex = -1;
  let playbackCurrentIndex = -1;
  let playbackToken = null;
  let playbackLastTick = null;

  function emit(type, detail) {
    const set = listeners.get(type);
    if (!set) return;
    for (const fn of set) {
      try {
        fn(detail);
      } catch {
        /* ignore */
      }
    }
  }

  function setMode(next) {
    if (mode === next) return;
    mode = next;
    const enableInteractions = mode === 'live';
    try {
      setInteractionsEnabled?.(enableInteractions);
    } catch {
      /* ignore */
    }
    emit('mode', { mode });
  }

  function stopPlayback({ forceLive = false } = {}) {
    if (playbackToken != null) {
      caf(playbackToken);
      playbackToken = null;
    }
    playbackLastTick = null;
    playbackStartIndex = -1;
    playbackCurrentIndex = -1;
    if (forceLive) setMode('live');
  }

  function handleStoreUpdate(snapshot) {
    bufferSnapshot = snapshot;
    if (mode === 'live' && playbackIndex !== snapshot.liveIndex) {
      playbackIndex = snapshot.liveIndex;
      emit('index', { index: playbackIndex });
    }
    emit('buffer', snapshot);
  }

  const unsubscribe = store.subscribe((snapshot) => {
    handleStoreUpdate(snapshot);
  });

  function applyFrameAt(index, { reason } = {}) {
    const frame = store.get(index);
    if (!frame) return false;
    applyFrame(frame, { reason });
    playbackIndex = index;
    emit('index', { index: playbackIndex });
    return true;
  }

  function schedulePlaybackTick() {
    if (playbackToken != null) caf(playbackToken);
    playbackToken = raf(playbackTick);
  }

  function playbackTick(ts) {
    if (mode !== 'playback') {
      playbackToken = null;
      return;
    }
    const snapshot = bufferSnapshot;
    const liveIdx = snapshot.liveIndex;
    if (playbackCurrentIndex >= liveIdx || liveIdx < 0) {
      stopPlayback();
      setMode('live');
      const resumeInfo = resumeMeta.lastRunningKind ? { ...resumeMeta } : null;
      emit('playback-complete', {
        resume: resumeInfo,
      });
      if (resumeInfo) {
        resumeMeta.lastRunningKind = null;
        resumeMeta.lastStartArgs = null;
      }
      paused = false;
      return;
    }
    const nextIndex = playbackCurrentIndex + 1;
    const nextFrame = store.get(nextIndex);
    if (!nextFrame) {
      // Underflow (evicted). Jump live.
      stopPlayback({ forceLive: true });
      const resumeInfo = resumeMeta.lastRunningKind ? { ...resumeMeta } : null;
      emit('playback-complete', {
        resume: resumeInfo,
        evicted: true,
      });
      if (resumeInfo) {
        resumeMeta.lastRunningKind = null;
        resumeMeta.lastStartArgs = null;
      }
      paused = false;
      return;
    }
    const targetDelta = 1000 / 20; // 20 fps playback
    if (playbackLastTick == null) playbackLastTick = ts - targetDelta;
    if (ts - playbackLastTick >= targetDelta) {
      playbackCurrentIndex = nextIndex;
      playbackLastTick = ts;
      applyFrame(nextFrame, { reason: 'timelinePlayback' });
      playbackIndex = playbackCurrentIndex;
      emit('index', { index: playbackIndex });
    }
    schedulePlaybackTick();
  }

  function scrubTo(index, { apply = true } = {}) {
    const snapshot = bufferSnapshot;
    if (!snapshot || snapshot.size === 0) return { applied: false };
    const clamped = Math.max(0, Math.min(index, snapshot.size - 1));
    stopPlayback();
    playbackStartIndex = -1;
    playbackCurrentIndex = clamped;
    const liveIdx = snapshot.liveIndex;
    const isLive = clamped === liveIdx;
    setMode(isLive ? 'live' : 'scrub');
    if (apply) {
      applyFrameAt(clamped, { reason: 'timelineScrub' });
    } else {
      playbackIndex = clamped;
      emit('index', { index: playbackIndex });
    }
    return { applied: apply, isLive };
  }

  function play() {
    const snapshot = bufferSnapshot;
    if (!snapshot || snapshot.size === 0) return { started: false };
    const startIndex = Math.max(0, Math.min(playbackIndex, snapshot.size - 1));
    const liveIdx = snapshot.liveIndex;
    if (startIndex >= liveIdx) {
      setMode('live');
      playbackIndex = liveIdx;
      emit('index', { index: playbackIndex });
      return { started: false, reason: 'already-live' };
    }
    paused = false;
    playbackStartIndex = startIndex;
    playbackCurrentIndex = startIndex;
    playbackLastTick = null;
    setMode('playback');
    applyFrameAt(startIndex, { reason: 'timelinePlaybackStart' });
    schedulePlaybackTick();
    return { started: true };
  }

  function beginPausedState({ runningKind, startArgs } = {}) {
    paused = true;
    resumeMeta.lastRunningKind = runningKind || null;
    resumeMeta.lastStartArgs = startArgs ? { ...startArgs } : null;
    stopPlayback();
    setMode('paused');
    emit('paused', { resume: { ...resumeMeta } });
  }

  function goLive() {
    const snapshot = bufferSnapshot;
    if (!snapshot || snapshot.size === 0) return { applied: false };
    const liveIdx = snapshot.liveIndex;
    stopPlayback({ forceLive: true });
    playbackIndex = liveIdx;
    paused = false;
    applyFrameAt(liveIdx, { reason: 'timelineLive' });
    setMode('live');
    const resumeInfo = resumeMeta.lastRunningKind ? { ...resumeMeta } : null;
    emit('playback-complete', { resume: resumeInfo, jumped: true });
    if (resumeInfo) {
      resumeMeta.lastRunningKind = null;
      resumeMeta.lastStartArgs = null;
    }
    return { applied: true };
  }

  function setPlaybackIndex(index) {
    playbackIndex = index;
    emit('index', { index: playbackIndex });
  }

  function debug() {
    const snapshot = bufferSnapshot;
    return {
      mode,
      playbackIndex,
      liveIndex: snapshot?.liveIndex ?? -1,
      size: snapshot?.size ?? 0,
      paused,
      resumeMeta: { ...resumeMeta },
    };
  }

  function on(type, fn) {
    if (typeof fn !== 'function') return () => {};
    let set = listeners.get(type);
    if (!set) {
      set = new Set();
      listeners.set(type, set);
    }
    set.add(fn);
    return () => set.delete(fn);
  }

  function destroy() {
    stopPlayback();
    unsubscribe?.();
    listeners.clear();
  }

  return {
    on,
    destroy,
    debug,
    getMode: () => mode,
    getPlaybackIndex: () => playbackIndex,
    setPlaybackIndex,
    scrubTo,
    play,
    goLive,
    beginPausedState,
  };
}

export default { createTimelineController };
