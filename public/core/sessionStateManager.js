// SessionStateManager centralises capture/restore of viewer, timeline, and websocket state.
// Dependencies are injected so the manager stays agnostic of Babylon.js internals.

const SCHEMA_VERSION = 4;

function migrateSnapshotSchema(snapshot) {
  if (!snapshot || typeof snapshot !== 'object') {
    throw new Error('Invalid session snapshot');
  }
  if (snapshot.schemaVersion === SCHEMA_VERSION) {
    return snapshot;
  }
  if (snapshot.schemaVersion === 3) {
    const clone = deepClone(snapshot);
    clone.schemaVersion = SCHEMA_VERSION;
    if (clone.timeline && Array.isArray(clone.timeline.controlMessages)) {
      clone.timeline.controlMessages = clone.timeline.controlMessages.map((msg, idx) => {
        if (!msg || typeof msg !== 'object') return null;
        const upgraded = { ...msg };
        if (typeof upgraded.label !== 'string' || !upgraded.label) {
          upgraded.label = typeof upgraded.id === 'string' ? upgraded.id : `control-${idx}`;
        }
        if (typeof upgraded.notes !== 'string') {
          delete upgraded.notes;
        }
        return upgraded;
      }).filter(Boolean);
    }
    return clone;
  }
  throw new Error(`Unsupported session snapshot schema: ${snapshot.schemaVersion}`);
}

function deepClone(obj) {
  try {
    return obj == null ? obj : JSON.parse(JSON.stringify(obj));
  } catch {
    return obj;
  }
}

function normalizeSource(meta = {}) {
  const kind = typeof meta.kind === 'string' ? meta.kind : 'unknown';
  const label = typeof meta.label === 'string' ? meta.label : null;
  const params = meta.params && typeof meta.params === 'object' ? { ...meta.params } : undefined;
  return { kind, ...(label ? { label } : {}), ...(params ? { params } : {}) };
}

function safeNumber(value, fallback = 0) {
  const n = Number(value);
  return Number.isFinite(n) ? n : fallback;
}

function normalizePlaybackConfig(playback = {}) {
  if (!playback || typeof playback !== 'object') return null;
  const cfg = {};
  if (typeof playback.startFrame === 'object' && playback.startFrame !== null) {
    const { frameId, offset, frameIndex } = playback.startFrame;
    cfg.startFrame = {};
    if (typeof frameId === 'string' && frameId) cfg.startFrame.frameId = frameId;
    if (Number.isFinite(offset)) cfg.startFrame.offset = Math.floor(offset);
    if (Number.isFinite(frameIndex)) cfg.startFrame.frameIndex = Math.max(0, Math.floor(frameIndex));
    if (!Object.keys(cfg.startFrame).length) delete cfg.startFrame;
  }
  if (typeof playback.autoPlay === 'boolean') cfg.autoPlay = playback.autoPlay;
  if (typeof playback.loop === 'boolean') cfg.loop = playback.loop;
  if (playback.loopRange && typeof playback.loopRange === 'object') {
    const { startFrameId, endFrameId, start, end } = playback.loopRange;
    cfg.loopRange = {};
    if (typeof startFrameId === 'string' && startFrameId) cfg.loopRange.startFrameId = startFrameId;
    if (typeof endFrameId === 'string' && endFrameId) cfg.loopRange.endFrameId = endFrameId;
    if (start && typeof start === 'object') {
      cfg.loopRange.start = {};
      if (Number.isFinite(start.offset)) cfg.loopRange.start.offset = Math.floor(start.offset);
      if (Number.isFinite(start.frameIndex)) cfg.loopRange.start.frameIndex = Math.max(0, Math.floor(start.frameIndex));
      if (typeof start.frameId === 'string' && start.frameId) cfg.loopRange.start.frameId = start.frameId;
      if (!Object.keys(cfg.loopRange.start).length) delete cfg.loopRange.start;
    }
    if (end && typeof end === 'object') {
      cfg.loopRange.end = {};
      if (Number.isFinite(end.offset)) cfg.loopRange.end.offset = Math.floor(end.offset);
      if (Number.isFinite(end.frameIndex)) cfg.loopRange.end.frameIndex = Math.max(0, Math.floor(end.frameIndex));
      if (typeof end.frameId === 'string' && end.frameId) cfg.loopRange.end.frameId = end.frameId;
      if (typeof end.inclusive === 'boolean') cfg.loopRange.end.inclusive = end.inclusive;
      if (!Object.keys(cfg.loopRange.end).length) delete cfg.loopRange.end;
    }
    if (!Object.keys(cfg.loopRange).length) delete cfg.loopRange;
  }
  if (Number.isFinite(playback.defaultFps)) cfg.defaultFps = Math.max(1, Math.round(playback.defaultFps));
  return Object.keys(cfg).length ? cfg : null;
}

export function createSessionStateManager(deps) {
  const {
    getViewerState,
    applyFullSnapshot,
    normalizeElement,
    energyPlot,
    frameBuffer,
    captureBaselineFromState,
    installBaseline,
    getInteractionCounters,
    setInteractionCounters,
    resetWsInitState,
    getWsClient,
    ensureWsInit,
    posToTriples,
    zOf,
    stateCellToArray,
    getVelocitiesForSnapshot,
    setMode,
    Mode,
    timelineState,
    timelineUiRef,
    clearTimelinePlayback,
    setTimelineOverlayVisible,
    setTimelineInteractionLock,
    applyTimelineFrame,
    refreshTimelineUi,
    setTimelineUiMode,
    lastContinuousOpts,
    rememberResume,
    stopSimulation,
    controlMessageEngine,
    timelinePlayback,
    getTimelinePlaybackSnapshot,
    getTimelineControlSnapshot,
    applyTimelinePlaybackSnapshot,
    applyControlMessageSnapshot,
    onSnapshotApplied,
  } = deps || {};

  let lastSource = { kind: 'init', label: 'initial' };
  let baselineSnapshot = null;

  function viewerToSnapshot(state = getViewerState()) {
    const elements = Array.isArray(state?.elements) ? state.elements.map(normalizeElement) : [];
    const positions = Array.isArray(state?.positions)
      ? state.positions.map((p) => [Number(p.x) || 0, Number(p.y) || 0, Number(p.z) || 0])
      : [];
    const velocities = Array.isArray(state?.dynamics?.velocities)
      ? state.dynamics.velocities.map((v) => Array.isArray(v)
        ? [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0]
        : [Number(v?.x) || 0, Number(v?.y) || 0, Number(v?.z) || 0])
      : undefined;
    const cell = stateCellToArray ? stateCellToArray(state?.cell) : null;
    return {
      elements,
      positions,
      velocities,
      cell,
      showCell: !!state?.showCell,
      showGhostCells: !!state?.showGhostCells,
    };
  }

  function captureSnapshot(meta = {}) {
    const state = getViewerState();
    const viewer = viewerToSnapshot(state);
    const energy = energyPlot?.exportSeries ? energyPlot.exportSeries() : { series: [], markerIndex: null };
    const frames = frameBuffer?.exportFrames ? frameBuffer.exportFrames() : [];
    const fbStats = frameBuffer?.stats ? frameBuffer.stats() : { capacity: frames.length };
    const counters = getInteractionCounters ? getInteractionCounters() : { user: 0, total: 0, lastApplied: 0 };
    let wsInfo = {};
    try {
      const ws = getWsClient?.();
      wsInfo = ws?.getState?.() || {};
    } catch { }
    const mode = deps?.getMode ? deps.getMode() : undefined;
    const runningMode = mode === Mode?.MD ? 'md' : mode === Mode?.Relax ? 'relax' : 'idle';
    const playbackSnapshot = getTimelinePlaybackSnapshot?.() ?? timelinePlayback?.getSnapshot?.();
    const controlSnapshot = getTimelineControlSnapshot?.() ?? controlMessageEngine?.getSnapshot?.();
    const snapshot = {
      schemaVersion: SCHEMA_VERSION,
      savedAt: new Date().toISOString(),
      source: normalizeSource(meta.kind ? meta : lastSource),
      viewer,
      energyPlot: energy,
      timeline: {
        capacity: fbStats.capacity,
        frames,
        lastLiveMode: runningMode,
        wasRunning: runningMode === 'md' || runningMode === 'relax',
        pendingSimParams:
          runningMode === 'md'
            ? deepClone(lastContinuousOpts?.md || {})
            : runningMode === 'relax'
              ? deepClone(lastContinuousOpts?.relax || {})
              : undefined,
        playback: playbackSnapshot ? normalizePlaybackConfig(playbackSnapshot) : null,
        controlMessages: Array.isArray(controlSnapshot) ? deepClone(controlSnapshot) : [],
      },
      websocket: {
        seq: Number(wsInfo.seq) || 0,
        nextSeq: Number(wsInfo.seq) || 0,
        clientAck: Number(wsInfo.clientAck) || Number(wsInfo.client_ack) || 0,
        userInteractionCount: Number(counters.user) || 0,
        totalInteractionCount: Number(counters.total) || Number(counters.user) || 0,
        simStep: Number(wsInfo.simStep ?? wsInfo.sim_step ?? 0) || 0,
      },
    };
    snapshot.websocket.nextSeq = (snapshot.websocket.seq | 0) + 1;
    return snapshot;
  }

  function setBaseline(snapshot, { reason, includeTimeline = true } = {}) {
    if (!snapshot) return;
    const clone = deepClone(snapshot);
    if (!includeTimeline && clone.timeline) {
      clone.timeline.frames = [];
      clone.timeline.wasRunning = false;
      clone.timeline.pendingSimParams = undefined;
    }
    baselineSnapshot = clone;
    try {
      const baseline = {
        elements: clone.viewer?.elements || [],
        positions: clone.viewer?.positions || [],
        velocities: clone.viewer?.velocities || null,
        cell: clone.viewer?.cell || null,
        showCell: !!clone.viewer?.showCell,
        showGhostCells: !!clone.viewer?.showGhostCells,
      };
      installBaseline?.(baseline, { reason: reason || clone.source?.kind || 'baseline' });
    } catch { }
  }

  function setBaselineFromState(meta, opts) {
    const snapshot = captureSnapshot(meta || lastSource);
    setBaseline(snapshot, opts);
    return snapshot;
  }

  function applyTimelineImport(snapshot) {
    const frames = snapshot?.timeline?.frames || [];
    frameBuffer?.importFrames(frames);
    energyPlot?.setMarker?.(null);
    if (frameBuffer && typeof frameBuffer.size === 'function' && frameBuffer.size() > 0) {
      const latest = frameBuffer.getByOffset(-1);
      if (latest && typeof latest.energyIndex === 'number') {
        energyPlot?.setMarker?.(latest.energyIndex);
      }
    }
    if (timelineUiRef?.value?.refresh) {
      timelineUiRef.value.refresh();
    } else {
      refreshTimelineUi?.();
    }
  }

  function deactivateTimelineUi() {
    clearTimelinePlayback?.();
    if (timelineState) {
      timelineState.active = false;
      timelineState.playing = false;
      timelineState.offset = -1;
      timelineState.suppressEnergy = false;
      timelineState.resumeMode = null;
    }
    setTimelineOverlayVisible?.(false);
    setTimelineInteractionLock?.(false);
    if (timelineUiRef?.value) {
      timelineUiRef.value.setMode?.('live');
      timelineUiRef.value.setActiveOffset?.(-1);
      timelineUiRef.value.refresh?.();
    } else {
      setTimelineUiMode?.('live');
      refreshTimelineUi?.();
    }
  }

  async function loadSnapshot(snapshot, { skipBaseline = false } = {}) {
    const normalized = migrateSnapshotSchema(deepClone(snapshot));
    const snap = normalized;

    lastSource = normalizeSource(snap.source || lastSource);
    stopSimulation?.();
    deactivateTimelineUi();

    const viewerPayload = {
      elements: snap.viewer?.elements,
      positions: snap.viewer?.positions,
      velocities: snap.viewer?.velocities,
      cell: snap.viewer?.cell,
    };
    await applyFullSnapshot?.(viewerPayload, { updateBaseline: false });

    const state = getViewerState();
    if (state) {
      state.showCell = !!snap.viewer?.showCell;
      state.showGhostCells = !!snap.viewer?.showGhostCells;
      state.markCellChanged?.();
    }

    energyPlot?.importSeries?.(snap.energyPlot || {});
    applyTimelineImport(snap);

    const counters = snap.websocket || {};
    const restoreUser = safeNumber(counters.userInteractionCount, 0);
    const restoreTotal = safeNumber(counters.totalInteractionCount, restoreUser);
    const restoreLastApplied = safeNumber(counters.lastApplied ?? counters.userInteractionCount, restoreUser);
    setInteractionCounters?.({
      user: restoreUser,
      total: restoreTotal,
      lastApplied: restoreLastApplied,
    });
    try {
      console.log('[SessionManager][loadSnapshot] counters restored', {
        user: restoreUser,
        total: restoreTotal,
        lastApplied: restoreLastApplied,
      });
    } catch { }

    if (!skipBaseline) {
      setBaseline(snap, { reason: snap.source?.kind || 'jsonLoad', includeTimeline: true });
    } else {
      setBaseline(snap, { reason: 'reset', includeTimeline: true });
    }

    applyTimelinePlaybackSnapshot?.(snap.timeline?.playback || null);
    resetWsInitState?.();

    try {
      const ws = getWsClient?.();
      if (ws?.seedSequencing) {
        const nextSeq = safeNumber(counters.nextSeq, counters.seq != null ? counters.seq + 1 : 1);
        const lastAck = safeNumber(counters.lastAck ?? counters.clientAck, 0);
        const simStep = safeNumber(counters.simStep, 0);
        ws.seedSequencing({
          nextSeq,
          ack: lastAck,
          userInteractionCount: restoreUser,
          simStep,
        });
        try {
          console.log('[SessionManager][loadSnapshot] seedSequencing applied', {
            nextSeq,
            ack: lastAck,
            userInteractionCount: restoreUser,
            simStep,
          });
        } catch { }
      }
      if (ws?.userInteraction) {
        const atomicNumbers = Array.isArray(state?.elements) ? state.elements.map(zOf) : [];
        const positions = posToTriples ? posToTriples(state) : [];
        const velocities = getVelocitiesForSnapshot ? getVelocitiesForSnapshot() : undefined;
        const cell = stateCellToArray ? stateCellToArray(state?.cell) : null;
        ws.userInteraction({
          atomic_numbers: atomicNumbers,
          positions,
          velocities,
          cell,
          natoms: atomicNumbers.length,
          full_update: true,
        });
        try {
          console.log('[SessionManager][loadSnapshot] sent full_update snapshot', {
            natoms: atomicNumbers.length,
          });
        } catch { }
      }
    } catch { }

    try {
      await ensureWsInit?.({ allowOffline: false });
      try {
        console.log('[SessionManager][loadSnapshot] ensureWsInit complete');
      } catch { }
    } catch { }

    if (frameBuffer && typeof frameBuffer.size === 'function' && frameBuffer.size() > 0) {
      if (timelineState) {
        timelineState.active = true;
        timelineState.suppressEnergy = true;
        timelineState.offset = -1;
        timelineState.resumeMode = snap.timeline?.wasRunning
          ? (snap.timeline?.lastLiveMode === 'relax' ? Mode?.Relax : Mode?.MD)
          : null;
      }
      setMode?.(Mode?.Timeline || 'timeline');
      setTimelineOverlayVisible?.(true);
      setTimelineInteractionLock?.(true);
      applyTimelineFrame?.(-1);
      if (timelineUiRef?.value) {
        timelineUiRef.value.setMode?.('paused');
        timelineUiRef.value.setActiveOffset?.(-1);
        timelineUiRef.value.refresh?.();
      } else {
        setTimelineUiMode?.('paused');
        refreshTimelineUi?.();
      }
    } else {
      setMode?.(Mode?.Idle || 'idle');
    }

    if (snap.timeline?.wasRunning) {
      const kind = snap.timeline.lastLiveMode === 'relax' ? 'relax' : 'md';
      if (kind === 'md') {
        lastContinuousOpts && (lastContinuousOpts.md = deepClone(snap.timeline.pendingSimParams || lastContinuousOpts?.md || {}));
      } else {
        lastContinuousOpts && (lastContinuousOpts.relax = deepClone(snap.timeline.pendingSimParams || lastContinuousOpts?.relax || {}));
      }
      rememberResume?.(kind, deepClone(snap.timeline.pendingSimParams || {}));
    }

    applyControlMessageSnapshot?.(snap.timeline?.controlMessages || []);
    onSnapshotApplied?.({
      source: snap.source,
      playback: snap.timeline?.playback || null,
      controlMessages: snap.timeline?.controlMessages || [],
    });

    return snap;
  }

  async function resetToLastSource() {
    if (!baselineSnapshot) {
      setBaselineFromState(lastSource);
    }
    if (!baselineSnapshot) return false;
    await loadSnapshot(deepClone(baselineSnapshot), { skipBaseline: true });
    return true;
  }

  async function loadFromFile(file) {
    let text;
    if (typeof file === 'string') {
      text = file;
    } else if (file && typeof file.text === 'function') {
      text = await file.text();
    } else {
      throw new Error('Unsupported file handle');
    }
    const parsed = JSON.parse(text);
    const applied = await loadSnapshot(parsed);
    setBaseline(applied, { includeTimeline: true, reason: applied.source?.kind || 'jsonLoad' });
    return applied;
  }

  function saveToFile({ filename, pretty = 2 } = {}) {
    const snapshot = captureSnapshot(lastSource);
    const json = JSON.stringify(snapshot, null, pretty);
    if (typeof document !== 'undefined') {
      const blob = new Blob([json], { type: 'application/json' });
      const link = document.createElement('a');
      const name = filename || `mlipview-session-${Date.now()}.json`;
      link.href = URL.createObjectURL(blob);
      link.download = name;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(() => URL.revokeObjectURL(link.href), 200);
    }
    return snapshot;
  }

  return {
    captureSnapshot,
    loadSnapshot,
    saveToFile,
    loadFromFile,
    resetToLastSource,
    noteSource(meta) {
      lastSource = normalizeSource(meta);
      setBaselineFromState(lastSource, { includeTimeline: false });
      return lastSource;
    },
    setBaselineFromState,
    setBaseline: (snapshot, opts) => setBaseline(snapshot, opts),
    getBaseline: () => deepClone(baselineSnapshot),
    getLastSource: () => ({ ...lastSource }),
  };
}

export default { createSessionStateManager };
