// SessionStateManager centralises capture/restore of viewer, timeline, and websocket state.
// Dependencies are injected so the manager stays agnostic of Babylon.js internals.

const SCHEMA_VERSION = 6;
const VALID_MESH_MODES = new Set(['solid', 'soft']);

type MeshMode = 'solid' | 'soft';

type MeshModeArray = MeshMode[];

interface MeshModeSnapshot {
  default?: MeshModeArray;
  current?: MeshModeArray;
}

interface MeshModesSnapshot {
  atoms?: MeshModeSnapshot;
  bonds?: MeshModeSnapshot;
}

interface PlaybackConfigSnapshot {
  autoPlay?: boolean;
  defaultFps?: number;
  loop?: boolean;
  loopRange?: unknown;
  startFrame?: unknown;
}

interface SnapshotMeta {
  kind?: string;
  label?: string;
  params?: Record<string, unknown>;
}

interface SourceMeta {
  kind: string;
  label?: string;
  params?: Record<string, unknown>;
}

interface BaselineOptions {
  reason?: string;
  includeTimeline?: boolean;
}

interface SaveSnapshotOptions {
  filename?: string;
  pretty?: number;
}

interface LoadSnapshotOptions {
  skipBaseline?: boolean;
}

function normalizeScalar(v, round = true) {
  const num = Number(v);
  if (!Number.isFinite(num)) return 0;
  return round ? Math.round(num) : num;
}

function normalizeVec3(value, round = true) {
  if (Array.isArray(value) && value.length >= 3) {
    return [
      normalizeScalar(value[0], round),
      normalizeScalar(value[1], round),
      normalizeScalar(value[2], round),
    ];
  }
  if (value && typeof value === 'object') {
    const { x, y, z } = value;
    return [normalizeScalar(x, round), normalizeScalar(y, round), normalizeScalar(z, round)];
  }
  return [0, 0, 0];
}

function serializeBonds(bonds) {
  if (!Array.isArray(bonds) || bonds.length === 0) {
    return {
      atomIndexA: [],
      atomIndexB: [],
      cellOffsetA: [],
      cellOffsetB: [],
      length: [],
      weight: [],
      opacity: [],
      flags: { inRing: [], crossing: [] },
      imageDelta: [],
    };
  }

  const payload = {
    atomIndexA: [],
    atomIndexB: [],
    cellOffsetA: [],
    cellOffsetB: [],
    length: [],
    weight: [],
    opacity: [],
    flags: { inRing: [], crossing: [] },
    imageDelta: [],
  };

  for (const bond of bonds) {
    const i = Number(bond?.i);
    const j = Number(bond?.j);
    payload.atomIndexA.push(Number.isInteger(i) ? i : 0);
    payload.atomIndexB.push(Number.isInteger(j) ? j : 0);
    payload.opacity.push(
      typeof bond?.opacity === 'number' && Number.isFinite(bond.opacity) ? bond.opacity : 1,
    );
    payload.length.push(
      typeof bond?.length === 'number' && Number.isFinite(bond.length) ? bond.length : null,
    );
    payload.weight.push(
      typeof bond?.weight === 'number' && Number.isFinite(bond.weight) ? bond.weight : null,
    );
    payload.flags.inRing.push(bond?.inRing ? 1 : 0);
    payload.flags.crossing.push(bond?.crossing ? 1 : 0);
    payload.cellOffsetA.push(normalizeVec3(bond?.cellOffsetA));
    payload.cellOffsetB.push(normalizeVec3(bond?.cellOffsetB));
    payload.imageDelta.push(normalizeVec3(bond?.imageDelta));
  }

  return payload;
}

function deserializeBonds(snapshotBonds) {
  if (!snapshotBonds) return [];

  if (Array.isArray(snapshotBonds)) {
    return snapshotBonds
      .map((bond) => {
        if (!bond || typeof bond !== 'object') return null;
        const i = Number(bond.i);
        const j = Number(bond.j);
        if (!Number.isInteger(i) || !Number.isInteger(j)) return null;
        return {
          i,
          j,
          length:
            typeof bond.length === 'number' && Number.isFinite(bond.length) ? bond.length : null,
          weight:
            typeof bond.weight === 'number' && Number.isFinite(bond.weight) ? bond.weight : null,
          opacity:
            typeof bond.opacity === 'number' && Number.isFinite(bond.opacity) ? bond.opacity : 1,
          inRing: !!bond.inRing,
          crossing: !!bond.crossing,
          imageDelta: normalizeVec3(bond.imageDelta),
          cellOffsetA: normalizeVec3(bond.cellOffsetA),
          cellOffsetB: normalizeVec3(bond.cellOffsetB),
        };
      })
      .filter(Boolean);
  }

  if (typeof snapshotBonds !== 'object' || snapshotBonds === null) return [];

  const {
    atomIndexA = [],
    atomIndexB = [],
    cellOffsetA = [],
    cellOffsetB = [],
    length = [],
    weight = [],
    opacity = [],
    flags = {},
    imageDelta = [],
  } = snapshotBonds;

  const inRing = Array.isArray(flags?.inRing) ? flags.inRing : [];
  const crossing = Array.isArray(flags?.crossing) ? flags.crossing : [];
  const count = Math.max(atomIndexA.length, atomIndexB.length);
  const result = [];

  for (let idx = 0; idx < count; idx++) {
    const i = Number(atomIndexA[idx]);
    const j = Number(atomIndexB[idx]);
    if (!Number.isInteger(i) || !Number.isInteger(j)) continue;
    result.push({
      i,
      j,
      length:
        typeof length[idx] === 'number' && Number.isFinite(length[idx]) ? Number(length[idx]) : null,
      weight:
        typeof weight[idx] === 'number' && Number.isFinite(weight[idx]) ? Number(weight[idx]) : null,
      opacity:
        typeof opacity[idx] === 'number' && Number.isFinite(opacity[idx])
          ? Number(opacity[idx])
          : 1,
      inRing: !!inRing[idx],
      crossing: !!crossing[idx],
      imageDelta: normalizeVec3(imageDelta[idx]),
      cellOffsetA: normalizeVec3(cellOffsetA[idx]),
      cellOffsetB: normalizeVec3(cellOffsetB[idx]),
    });
  }

  return result;
}

function migrateSnapshotSchema(snapshot: any) {
  if (!snapshot || typeof snapshot !== 'object') {
    throw new Error('Invalid session snapshot');
  }
  let current = deepClone(snapshot);
  let version = Number(current.schemaVersion);
  if (!Number.isFinite(version)) version = 3;
  if (version > SCHEMA_VERSION) {
    throw new Error(
      `Unsupported session snapshot schema: ${version} (max supported ${SCHEMA_VERSION})`,
    );
  }

  while (version < SCHEMA_VERSION) {
    switch (version) {
      case 3: {
        if (current.timeline && Array.isArray(current.timeline.controlMessages)) {
          current.timeline.controlMessages = current.timeline.controlMessages
            .map((msg, idx) => {
              if (!msg || typeof msg !== 'object') return null;
              const upgraded = { ...msg };
              if (typeof upgraded.label !== 'string' || !upgraded.label) {
                upgraded.label = typeof upgraded.id === 'string' ? upgraded.id : `control-${idx}`;
              }
              if (typeof upgraded.notes !== 'string') {
                delete upgraded.notes;
              }
              return upgraded;
            })
            .filter(Boolean);
        }
        version = 4;
        current.schemaVersion = version;
        break;
      }
      case 4: {
        current.viewer = current.viewer || {};
        const atomCount = Array.isArray(current.viewer.positions) ? current.viewer.positions.length : 0;
        if (!current.viewer.meshAssignments || typeof current.viewer.meshAssignments !== 'object') {
          current.viewer.meshAssignments = {
            atoms: atomCount ? new Array(atomCount).fill('solid') : [],
            bonds: [],
          };
        } else {
          if (!Array.isArray(current.viewer.meshAssignments.atoms)) {
            current.viewer.meshAssignments.atoms = atomCount ? new Array(atomCount).fill('solid') : [];
          }
          if (!Array.isArray(current.viewer.meshAssignments.bonds)) {
            current.viewer.meshAssignments.bonds = [];
          }
        }
        current.render = current.render || {};
        if (!current.render.overrides || typeof current.render.overrides !== 'object') {
          current.render.overrides = {};
        }
        if (!current.render.overrides.atoms) current.render.overrides.atoms = {};
        if (!current.render.overrides.bonds) current.render.overrides.bonds = {};
        current.schemaVersion = 5;
        break;
      }
      case 5: {
        current.viewer = current.viewer || {};
        if (!Array.isArray(current.viewer.periodicBonds)) {
          current.viewer.periodicBonds = [];
        } else {
          current.viewer.periodicBonds = current.viewer.periodicBonds
            .map((bond) => {
              if (!bond || typeof bond !== 'object') return null;
              const i = Number(bond.i);
              const j = Number(bond.j);
              if (!Number.isInteger(i) || !Number.isInteger(j)) return null;
              const opacity =
                typeof bond.opacity === 'number' && Number.isFinite(bond.opacity)
                  ? bond.opacity
                  : 1;
              const imageDelta = Array.isArray(bond.imageDelta)
                ? bond.imageDelta.slice(0, 3).map((v) =>
                    Number.isFinite(Number(v)) ? Math.round(Number(v)) : 0
                  )
                : [0, 0, 0];
              return {
                i,
                j,
                opacity,
                crossing: !!bond.crossing,
                imageDelta,
              };
            })
            .filter(Boolean);
        }
        if (!Array.isArray(current.viewer.ghostAtoms)) {
          current.viewer.ghostAtoms = [];
        }
        if (!Array.isArray(current.viewer.ghostBondMeta)) {
          current.viewer.ghostBondMeta = [];
        }
        const legacyBonds = normaliseLegacyPeriodicBonds(current.viewer.periodicBonds);
        current.viewer.bonds = serializeBonds(legacyBonds);
        delete current.viewer.periodicBonds;
        version = 6;
        current.schemaVersion = version;
        break;
      }
      default:
        throw new Error(`Unsupported session snapshot schema: ${version}`);
    }
  }

  const normalizedBonds = (() => {
    if (current.viewer?.periodicBonds) {
      const legacy = normaliseLegacyPeriodicBonds(current.viewer.periodicBonds);
      delete current.viewer.periodicBonds;
      return legacy;
    }
    return deserializeBonds(current.viewer?.bonds);
  })();
  current.viewer = current.viewer || {};
  current.viewer.bonds = serializeBonds(normalizedBonds);
  current.schemaVersion = SCHEMA_VERSION;
  return current;
}

function deepClone(obj: any) {
  try {
    return obj == null ? obj : JSON.parse(JSON.stringify(obj));
  } catch {
    return obj;
  }
}

function normalizeSource(meta: SnapshotMeta = {} as SnapshotMeta): SourceMeta {
  const kind = typeof meta.kind === 'string' && meta.kind ? meta.kind : 'unknown';
  const label = typeof meta.label === 'string' && meta.label ? meta.label : undefined;
  const params =
    meta.params && typeof meta.params === 'object'
      ? { ...meta.params }
      : undefined;
  const source: SourceMeta = { kind };
  if (label) source.label = label;
  if (params) source.params = params;
  return source;
}

function safeNumber(value, fallback = 0) {
  const n = Number(value);
  return Number.isFinite(n) ? n : fallback;
}

function normalizePlaybackConfig(playback: any = {}) {
  if (!playback || typeof playback !== 'object') return null;
  const cfg: any = {};
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

function normaliseLegacyPeriodicBonds(periodicBonds) {
  if (!Array.isArray(periodicBonds)) return [];
  return periodicBonds
    .map((bond) => {
      if (!bond || typeof bond !== 'object') return null;
      const i = Number(bond.i);
      const j = Number(bond.j);
      if (!Number.isInteger(i) || !Number.isInteger(j)) return null;
      const imageDelta = Array.isArray(bond.imageDelta)
        ? bond.imageDelta.slice(0, 3).map((v) => {
            const num = Number(v);
            if (!Number.isFinite(num)) return 0;
            const rounded = Math.round(num);
            return Number.isFinite(rounded) ? rounded : 0;
          })
        : [0, 0, 0];
      return {
        i,
        j,
        length: typeof bond.length === 'number' ? bond.length : null,
        weight:
          typeof bond.weight === 'number' && Number.isFinite(bond.weight) ? bond.weight : null,
        opacity:
          typeof bond.opacity === 'number' && Number.isFinite(bond.opacity) ? bond.opacity : 1,
        inRing: !!bond.inRing,
        crossing: !!bond.crossing,
        imageDelta,
        cellOffsetA: [0, 0, 0],
        cellOffsetB: imageDelta,
      };
    })
    .filter(Boolean);
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
    getMeshModes,
    applyMeshModes,
    clearOpacityMask,
    controlMessageEngine,
    timelinePlayback,
    getTimelinePlaybackSnapshot,
    getTimelineControlSnapshot,
    applyTimelinePlaybackSnapshot,
    applyControlMessageSnapshot,
    onSnapshotApplied,
  } = deps || {};

  let lastSource: SourceMeta = { kind: 'init', label: 'initial' };
  let baselineSnapshot = null;

  function sanitizeModeArray(
    array: MeshModeArray | null | undefined,
    fallback: MeshMode = 'solid',
    expectedLength?: number,
  ): MeshModeArray | null {
    if (!Array.isArray(array)) return null;
    if (
      typeof expectedLength === 'number' &&
      expectedLength >= 0 &&
      array.length !== expectedLength
    ) {
      return null;
    }
    return array.map((value) => (VALID_MESH_MODES.has(value) ? (value as MeshMode) : fallback));
  }

  function buildOverrides(
    defaultModes: MeshModeArray = [],
    currentModes: MeshModeArray = [],
  ): Record<MeshMode, number[]> {
    const overrides: Record<MeshMode, number[]> = { soft: [], solid: [] };
    if (!Array.isArray(defaultModes)) return overrides;
    for (let i = 0; i < defaultModes.length; i++) {
      const defaultMode = VALID_MESH_MODES.has(defaultModes[i])
        ? (defaultModes[i] as MeshMode)
        : 'solid';
      const currentMode = VALID_MESH_MODES.has(currentModes[i])
        ? (currentModes[i] as MeshMode)
        : defaultMode;
      if (currentMode !== defaultMode) overrides[currentMode].push(i);
    }
    return overrides;
  }

  function deriveCurrentModes(
    defaultModes: MeshModeArray | null,
    overrides: Record<string, number[]> | undefined,
  ): MeshModeArray | null {
    if (!Array.isArray(defaultModes)) return null;
    const current = defaultModes.slice();
    if (!overrides || typeof overrides !== 'object') return current;
    for (const [mode, indices] of Object.entries(overrides)) {
      if (!VALID_MESH_MODES.has(mode)) continue;
      if (!Array.isArray(indices)) continue;
      for (const idx of indices) {
        const n = Number(idx);
        if (Number.isInteger(n) && n >= 0 && n < current.length) {
          current[n] = mode as MeshMode;
        }
      }
    }
    return current;
  }

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
    const bonds = serializeBonds(Array.isArray(state?.bonds) ? state.bonds : []);
    const ghostAtoms = Array.isArray(state?.ghostImages)
      ? state.ghostImages
          .map((ghost) => {
            if (!ghost || !Number.isInteger(ghost.atomIndex)) return null;
            const position = Array.isArray(ghost.position)
              ? ghost.position.slice(0, 3).map((v) => Number(v) || 0)
              : [
                  Number(ghost.position?.x) || 0,
                  Number(ghost.position?.y) || 0,
                  Number(ghost.position?.z) || 0,
                ];
            return {
              atomIndex: ghost.atomIndex,
              shift: Array.isArray(ghost.shift)
                ? ghost.shift.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              position,
            };
          })
          .filter(Boolean)
      : [];
    const ghostBondMeta = Array.isArray(state?.ghostBondMeta)
      ? state.ghostBondMeta
          .map((meta) => {
            if (!meta || !meta.base) return null;
            const i = Number(meta.base.i);
            const j = Number(meta.base.j);
            if (!Number.isInteger(i) || !Number.isInteger(j)) return null;
            return {
              base: { i, j },
              shiftA: Array.isArray(meta.shiftA)
                ? meta.shiftA.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              shiftB: Array.isArray(meta.shiftB)
                ? meta.shiftB.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              imageDelta: Array.isArray(meta.imageDelta)
                ? meta.imageDelta.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
            };
          })
          .filter(Boolean)
      : [];
    return {
      elements,
      positions,
      velocities,
      cell,
      showCell: !!state?.showCell,
      showGhostCells: !!state?.showGhostCells,
      bonds,
      ghostAtoms,
      ghostBondMeta,
    };
  }

  function captureSnapshot(meta?: SnapshotMeta) {
    const state: any = getViewerState();
    const viewer = viewerToSnapshot(state);
    const energy = energyPlot?.exportSeries ? energyPlot.exportSeries() : { series: [], markerIndex: null };
    const frames = frameBuffer?.exportFrames ? frameBuffer.exportFrames() : [];
    const fbStats = frameBuffer?.stats ? frameBuffer.stats() : { capacity: frames.length };
    const counters = getInteractionCounters ? getInteractionCounters() : { user: 0, total: 0, lastApplied: 0 };
    let wsInfo: any = {};
    try {
      const ws = getWsClient?.();
      wsInfo = ws?.getState?.() || {};
    } catch { }
    const mode = deps?.getMode ? deps.getMode() : undefined;
    const runningMode = mode === Mode?.MD ? 'md' : mode === Mode?.Relax ? 'relax' : 'idle';
    const playbackSnapshot =
      (getTimelinePlaybackSnapshot?.() ??
        timelinePlayback?.getSnapshot?.()) as PlaybackConfigSnapshot | null;
    const controlSnapshot = getTimelineControlSnapshot?.() ?? controlMessageEngine?.getSnapshot?.();
    const meshModes = (getMeshModes?.() || {}) as MeshModesSnapshot;
    const atomCount = viewer.positions.length;
    const bondCount = Array.isArray(state?.bonds) ? state.bonds.length : 0;
    const atomDefaults =
      sanitizeModeArray(meshModes.atoms?.default, 'solid', atomCount) ||
      (atomCount ? new Array(atomCount).fill('solid') : []);
    const atomCurrents =
      sanitizeModeArray(meshModes.atoms?.current, 'solid', atomCount) || atomDefaults.slice();
    const bondDefaults =
      sanitizeModeArray(meshModes.bonds?.default, 'solid', bondCount) ||
      (bondCount ? new Array(bondCount).fill('solid') : []);
    const bondCurrents =
      sanitizeModeArray(meshModes.bonds?.current, 'solid', bondCount) || bondDefaults.slice();

    const atomOverridesFull = buildOverrides(atomDefaults || [], atomCurrents || []);
    const bondOverridesFull = buildOverrides(bondDefaults || [], bondCurrents || []);

    const atomOverrides: Record<string, number[]> = {};
    for (const mode of Object.keys(atomOverridesFull)) {
      if (atomOverridesFull[mode].length) atomOverrides[mode] = atomOverridesFull[mode];
    }
    const bondOverrides: Record<string, number[]> = {};
    for (const mode of Object.keys(bondOverridesFull)) {
      if (bondOverridesFull[mode].length) bondOverrides[mode] = bondOverridesFull[mode];
    }

    const snapshot: any = {
      schemaVersion: SCHEMA_VERSION,
      savedAt: new Date().toISOString(),
      source: normalizeSource(meta?.kind ? meta : lastSource),
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
    snapshot.viewer.meshAssignments = {
      atoms: atomDefaults,
      bonds: bondDefaults,
    };
    if (!snapshot.render || typeof snapshot.render !== 'object') snapshot.render = {};
    const overridesPayload: { atoms?: Record<string, number[]>; bonds?: Record<string, number[]> } = {};
    if (Object.keys(atomOverrides).length) overridesPayload.atoms = atomOverrides;
    if (Object.keys(bondOverrides).length) overridesPayload.bonds = bondOverrides;
    if (Object.keys(overridesPayload).length) snapshot.render.overrides = overridesPayload;
    else if (snapshot.render.overrides) delete snapshot.render.overrides;
    if (snapshot.render && !Object.keys(snapshot.render).length) delete snapshot.render;
    return snapshot;
  }

  function setBaseline(snapshot: any, options: BaselineOptions = {}) {
    const { reason, includeTimeline = true } = options;
    if (!snapshot) return;
    const clone: any = deepClone(snapshot);
    if (!includeTimeline && clone.timeline) {
      clone.timeline.frames = [];
      clone.timeline.wasRunning = false;
      clone.timeline.pendingSimParams = undefined;
    }
    const normalizedBonds = deserializeBonds(clone.viewer?.bonds);
    clone.viewer.bonds = serializeBonds(normalizedBonds);
    baselineSnapshot = clone;
    try {
      const baseline = {
        elements: clone.viewer?.elements || [],
        positions: clone.viewer?.positions || [],
        velocities: clone.viewer?.velocities || null,
        cell: clone.viewer?.cell || null,
        showCell: !!clone.viewer?.showCell,
        showGhostCells: !!clone.viewer?.showGhostCells,
        bonds: normalizedBonds,
        ghostAtoms: Array.isArray(clone.viewer?.ghostAtoms)
          ? clone.viewer.ghostAtoms.map((ghost) => ({
              atomIndex: Number(ghost.atomIndex) || 0,
              shift: Array.isArray(ghost.shift)
                ? ghost.shift.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              position: Array.isArray(ghost.position)
                ? ghost.position.slice(0, 3).map((v) => Number(v) || 0)
                : [
                    Number(ghost.position?.x) || 0,
                    Number(ghost.position?.y) || 0,
                    Number(ghost.position?.z) || 0,
                  ],
            }))
          : [],
        ghostBondMeta: Array.isArray(clone.viewer?.ghostBondMeta)
          ? clone.viewer.ghostBondMeta.map((meta) => ({
              base: {
                i: Number(meta.base?.i) || 0,
                j: Number(meta.base?.j) || 0,
              },
              shiftA: Array.isArray(meta.shiftA)
                ? meta.shiftA.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              shiftB: Array.isArray(meta.shiftB)
                ? meta.shiftB.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
              imageDelta: Array.isArray(meta.imageDelta)
                ? meta.imageDelta.slice(0, 3).map((v) => Number(v) || 0)
                : [0, 0, 0],
            }))
          : [],
      };
      installBaseline?.(baseline, { reason: reason || clone.source?.kind || 'baseline' });
    } catch { }
  }

  function setBaselineFromState(meta?: SnapshotMeta, opts: BaselineOptions = {}) {
    const snapshot = captureSnapshot(meta ?? lastSource);
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

  async function loadSnapshot(snapshot: any, options: LoadSnapshotOptions = {}) {
    const { skipBaseline = false } = options;
    const normalized: any = migrateSnapshotSchema(deepClone(snapshot));
    const snap: any = normalized;
    const normalizedBonds = deserializeBonds(snap.viewer?.bonds);
    snap.viewer.bonds = serializeBonds(normalizedBonds);

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
      state.ghostImages = Array.isArray(snap.viewer?.ghostAtoms)
        ? snap.viewer.ghostAtoms.map((ghost) => ({
            atomIndex: Number(ghost.atomIndex) || 0,
            shift: Array.isArray(ghost.shift)
              ? ghost.shift.slice(0, 3).map((v) => Number(v) || 0)
              : [0, 0, 0],
            position: Array.isArray(ghost.position)
              ? ghost.position.slice(0, 3).map((v) => Number(v) || 0)
              : [
                  Number(ghost.position?.x) || 0,
                  Number(ghost.position?.y) || 0,
                  Number(ghost.position?.z) || 0,
                ],
          }))
        : [];
      state.ghostBondMeta = Array.isArray(snap.viewer?.ghostBondMeta)
        ? snap.viewer.ghostBondMeta.map((meta) => ({
            base: {
              i: Number(meta.base?.i) || 0,
              j: Number(meta.base?.j) || 0,
            },
            shiftA: Array.isArray(meta.shiftA)
              ? meta.shiftA.slice(0, 3).map((v) => Number(v) || 0)
              : [0, 0, 0],
            shiftB: Array.isArray(meta.shiftB)
              ? meta.shiftB.slice(0, 3).map((v) => Number(v) || 0)
              : [0, 0, 0],
            imageDelta: Array.isArray(meta.imageDelta)
              ? meta.imageDelta.slice(0, 3).map((v) => Number(v) || 0)
              : [0, 0, 0],
          }))
        : [];
      state.bonds = normalizedBonds;
      state.markBondsChanged?.();
    }

    clearOpacityMask?.({ refresh: true });
    const meshAssignments = snap.viewer?.meshAssignments || {};
    const overrides = snap.render?.overrides || {};
    const atomCount = Array.isArray(snap.viewer?.positions) ? snap.viewer.positions.length : 0;
    const bondCount = Array.isArray(meshAssignments.bonds) ? meshAssignments.bonds.length : 0;
    const atomDefaults =
      sanitizeModeArray(meshAssignments.atoms, 'solid', atomCount) ||
      (atomCount ? new Array(atomCount).fill('solid') : []);
    const bondDefaults =
      sanitizeModeArray(meshAssignments.bonds, 'solid', bondCount) ||
      (bondCount ? new Array(bondCount).fill('solid') : []);
    const atomCurrents = deriveCurrentModes(atomDefaults, overrides.atoms);
    const bondCurrents = deriveCurrentModes(bondDefaults, overrides.bonds);
    applyMeshModes?.({
      atomDefault: atomDefaults || undefined,
      atomCurrent: atomCurrents || undefined,
      bondDefault: bondDefaults || undefined,
      bondCurrent: bondCurrents || undefined,
    });

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
        const lastSeq = safeNumber(
          counters.seq,
          Number.isFinite(nextSeq) ? Math.max(0, Number(nextSeq) - 1) : 0,
        );
        const lastAck = safeNumber(counters.lastAck ?? counters.clientAck, 0);
        const simStep = safeNumber(counters.simStep, 0);
        ws.seedSequencing({
          nextSeq,
          lastSeq,
          ack: lastAck,
          userInteractionCount: restoreUser,
          simStep,
        });
        try {
          console.log('[SessionManager][loadSnapshot] seedSequencing applied', {
            nextSeq,
            lastSeq,
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

  async function loadFromFile(file: string | File) {
    let text: string;
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

  function saveToFile({ filename, pretty = 2 }: SaveSnapshotOptions = {}) {
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
    noteSource(meta: SnapshotMeta) {
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
