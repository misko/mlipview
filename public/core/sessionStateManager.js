// SessionStateManager centralises capture/restore of viewer, timeline, and websocket state.
// Dependencies are injected so the manager stays agnostic of Babylon.js internals.

const SCHEMA_VERSION = 1;

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
    if (!snapshot || snapshot.schemaVersion !== SCHEMA_VERSION) {
      throw new Error('Unsupported session snapshot schema');
    }

    lastSource = normalizeSource(snapshot.source || lastSource);
    stopSimulation?.();
    deactivateTimelineUi();

    const viewerPayload = {
      elements: snapshot.viewer?.elements,
      positions: snapshot.viewer?.positions,
      velocities: snapshot.viewer?.velocities,
      cell: snapshot.viewer?.cell,
    };
    await applyFullSnapshot?.(viewerPayload, { updateBaseline: false });

    const state = getViewerState();
    if (state) {
      state.showCell = !!snapshot.viewer?.showCell;
      state.showGhostCells = !!snapshot.viewer?.showGhostCells;
      state.markCellChanged?.();
    }

    energyPlot?.importSeries?.(snapshot.energyPlot || {});
    applyTimelineImport(snapshot);

    const counters = snapshot.websocket || {};
    const timelineFrames = snapshot?.timeline?.frames || [];
    const lastFrame = timelineFrames.length ? timelineFrames[timelineFrames.length - 1] : null;
    const userCount = Number(counters.userInteractionCount);
    const totalCount = Number(counters.totalInteractionCount);
    const fallbackUser = Number.isFinite(userCount) ? userCount : (Number.isFinite(lastFrame?.userInteractionCount) ? lastFrame.userInteractionCount : 0);
    const lastAppliedCount = Number.isFinite(lastFrame?.userInteractionCount)
      ? lastFrame.userInteractionCount
      : fallbackUser;
    setInteractionCounters?.({
      user: fallbackUser,
      total: Number.isFinite(totalCount) ? totalCount : fallbackUser,
      lastApplied: lastAppliedCount,
    });

    if (!skipBaseline) {
      setBaseline(snapshot, { reason: snapshot.source?.kind || 'jsonLoad', includeTimeline: true });
    } else {
      setBaseline(snapshot, { reason: 'reset', includeTimeline: true });
    }

    resetWsInitState?.();

    try {
      const ws = getWsClient?.();
      if (ws?.seedSequencing) {
        ws.seedSequencing({
          nextSeq: Number(counters.seq ?? counters.nextSeq) || 0,
          ack: Number(counters.clientAck) || 0,
          userInteractionCount: Number(counters.userInteractionCount) || 0,
          simStep: Number(counters.simStep) || 0,
        });
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
      }
    } catch { }

    try {
      await ensureWsInit?.({ allowOffline: false });
    } catch { }

    if (frameBuffer && typeof frameBuffer.size === 'function' && frameBuffer.size() > 0) {
      if (timelineState) {
        timelineState.active = true;
        timelineState.suppressEnergy = true;
        timelineState.offset = -1;
        timelineState.resumeMode = snapshot.timeline?.wasRunning
          ? (snapshot.timeline?.lastLiveMode === 'relax' ? Mode?.Relax : Mode?.MD)
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

    if (snapshot.timeline?.wasRunning) {
      const kind = snapshot.timeline.lastLiveMode === 'relax' ? 'relax' : 'md';
      if (kind === 'md') {
        lastContinuousOpts && (lastContinuousOpts.md = deepClone(snapshot.timeline.pendingSimParams || lastContinuousOpts?.md || {}));
      } else {
        lastContinuousOpts && (lastContinuousOpts.relax = deepClone(snapshot.timeline.pendingSimParams || lastContinuousOpts?.relax || {}));
      }
      rememberResume?.(kind, deepClone(snapshot.timeline.pendingSimParams || {}));
    }

    return snapshot;
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
    const snapshot = JSON.parse(text);
    await loadSnapshot(snapshot);
    setBaseline(snapshot, { includeTimeline: true, reason: snapshot.source?.kind || 'jsonLoad' });
    return snapshot;
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
