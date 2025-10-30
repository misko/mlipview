// index.js — streamlined viewer orchestration with WS/protobuf backend
// Refactored: unified streaming loop, idle listener, predictable throttle, lighter debounces.

import { createMoleculeState } from './domain/moleculeState.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { createScene } from './render/scene.js';
import { createMoleculeView } from './render/moleculeView.js';
import { createPickingService } from './core/pickingService.js';
import installTouchControls from './ui/touchControls.js';
import { createManipulationService } from './domain/manipulationService.js';
import { createVRSupport } from './vr/setup.js';
import { createVRPicker } from './vr/vr-picker.js';
import { __count } from './util/funcCount.js';
import { DEFAULT_MD_FRICTION, DEFAULT_MIN_STEP_INTERVAL_MS, MAX_STEP } from './util/constants.js';
import { getWS } from './fairchem_ws_client.js';
import { showErrorBanner, hideErrorBanner } from './ui/errorBanner.js';
import { SYMBOL_TO_Z, Z_TO_SYMBOL, OMOL25_ELEMENTS } from './data/periodicTable.js';
import { createFrameBuffer } from './core/frameBuffer.js';
import { installTimeline } from './ui/timeline.js';
import { createSessionStateManager } from './core/sessionStateManager.js';
import { createTimelinePlaybackController } from './core/timelinePlaybackController.js';
import { createControlMessageEngine } from './core/controlMessageEngine.js';
import { createCalloutLayer } from './render/calloutLayer.js';
import { createTimelineEditorPanel } from './ui/timelineEditorPanel.js';
import { createTimelineEditorStatus } from './ui/timelineEditorStatus.js';

/* ──────────────────────────────────────────────────────────────────────────
   Tunables & Defaults
   ────────────────────────────────────────────────────────────────────────── */

const DEFAULTS = {
  DRAG_THROTTLE_MS: 60,
  POS_ENERGY_DEBOUNCE_MS: 50,
  LATCH_AFTER_DRAG_MS: 600,
  LATCH_AFTER_VR_MS: 800,
  RPS_WINDOW_MS: 2000,
  CONTINUOUS_STEPS: 1000,
  SAFE_SPHERE_RADIUS: 20,
  ROTATION_LATCH_MS: 600,
  BOND_SCROLL_STEP_RAD: Math.PI / 36,
  WS_CONNECT_TIMEOUT_MS: 15000,
};

const TEMP_PLACEHOLDER = 'T: — K';

const wsStateLog = [];
const reconnectState = {
  attempts: 0,
  nextAttemptAt: 0,
  countdownTimer: null,
  bannerVisible: false,
};

const OMOL25_SET = new Set(Array.from(OMOL25_ELEMENTS || []));

const pendingResume = { kind: null, opts: null };
const lastContinuousOpts = { md: null, relax: null };

function cfgNumber(v, fallback) {
  const n = Number(v);
  return Number.isFinite(n) ? n : fallback;
}

const TUNABLES = {
  get DRAG_THROTTLE_MS() {
    if (typeof window === 'undefined') return DEFAULTS.DRAG_THROTTLE_MS;
    return cfgNumber(window.__MLIP_CONFIG?.dragThrottleMs, DEFAULTS.DRAG_THROTTLE_MS);
  },
  get POS_ENERGY_DEBOUNCE_MS() {
    if (typeof window === 'undefined') return DEFAULTS.POS_ENERGY_DEBOUNCE_MS;
    return cfgNumber(window.__MLIP_CONFIG?.posEnergyDebounceMs, DEFAULTS.POS_ENERGY_DEBOUNCE_MS);
  },
  get LATCH_AFTER_DRAG_MS() {
    if (typeof window === 'undefined') return DEFAULTS.LATCH_AFTER_DRAG_MS;
    return cfgNumber(window.__MLIP_CONFIG?.latchAfterDragMs, DEFAULTS.LATCH_AFTER_DRAG_MS);
  },
  get LATCH_AFTER_VR_MS() {
    if (typeof window === 'undefined') return DEFAULTS.LATCH_AFTER_VR_MS;
    return cfgNumber(window.__MLIP_CONFIG?.latchAfterVrMs, DEFAULTS.LATCH_AFTER_VR_MS);
  },
  get RPS_WINDOW_MS() {
    if (typeof window === 'undefined') return DEFAULTS.RPS_WINDOW_MS;
    return cfgNumber(window.__MLIP_CONFIG?.rpsWindowMs, DEFAULTS.RPS_WINDOW_MS);
  },
  get SAFE_SPHERE_RADIUS() {
    if (typeof window === 'undefined') return DEFAULTS.SAFE_SPHERE_RADIUS;
    return cfgNumber(window.__MLIP_CONFIG?.safeSphereRadius, DEFAULTS.SAFE_SPHERE_RADIUS);
  },
  get ROTATION_LATCH_MS() {
    if (typeof window === 'undefined') return DEFAULTS.ROTATION_LATCH_MS;
    return cfgNumber(window.__MLIP_CONFIG?.rotationLatchMs, DEFAULTS.ROTATION_LATCH_MS);
  },
  get BOND_SCROLL_STEP_RAD() {
    if (typeof window === 'undefined') return DEFAULTS.BOND_SCROLL_STEP_RAD;
    return cfgNumber(window.__MLIP_CONFIG?.bondScrollStepRad, DEFAULTS.BOND_SCROLL_STEP_RAD);
  },
  get WS_CONNECT_TIMEOUT_MS() {
    if (typeof window === 'undefined') return DEFAULTS.WS_CONNECT_TIMEOUT_MS;
    return cfgNumber(window.__MLIP_CONFIG?.wsConnectTimeoutMs, DEFAULTS.WS_CONNECT_TIMEOUT_MS);
  },
};

const wsStateListeners = new Set();
function addWsStateListener(fn) {
  if (typeof fn === 'function') wsStateListeners.add(fn);
  return () => wsStateListeners.delete(fn);
}

function dispatchWsState(evt) {
  const snapshot = { ...(evt || {}), timestamp: Date.now() };
  wsStateLog.push(snapshot);
  while (wsStateLog.length > 200) wsStateLog.shift();
  for (const fn of wsStateListeners) {
    try { fn(snapshot); } catch { }
  }
}

try {
  if (typeof window !== 'undefined' && !window.__MLIPVIEW_WS_STATE_BRIDGE__) {
    const previous = typeof window.__WS_ON_STATE__ === 'function' ? window.__WS_ON_STATE__ : null;
    const bridged = function wsStateBridge(evt) {
      dispatchWsState(evt || {});
      if (previous) {
        try { previous(evt); } catch { }
      }
    };
    bridged.__MLIPVIEW_BRIDGED__ = true;
    window.__WS_ON_STATE__ = bridged;
    window.__MLIPVIEW_WS_STATE_BRIDGE__ = true;
  }
} catch { }

try {
  console.log('[config:tunables]', {
    dragThrottleMs: TUNABLES.DRAG_THROTTLE_MS,
    posEnergyDebounceMs: TUNABLES.POS_ENERGY_DEBOUNCE_MS,
    latchAfterDragMs: TUNABLES.LATCH_AFTER_DRAG_MS,
    latchAfterVrMs: TUNABLES.LATCH_AFTER_VR_MS,
    rpsWindowMs: TUNABLES.RPS_WINDOW_MS,
    bondScrollStepRad: TUNABLES.BOND_SCROLL_STEP_RAD,
    wsConnectTimeoutMs: TUNABLES.WS_CONNECT_TIMEOUT_MS,
  });
} catch { }

/* ──────────────────────────────────────────────────────────────────────────
   Utils
   ────────────────────────────────────────────────────────────────────────── */

const env = {
  now: () => (performance?.now?.() ?? Date.now()),
  apiOn: () => (typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API),
};

const dbg = {
  log: (...a) => { try { console.log(...a); } catch { } },
  warn: (...a) => { try { console.warn(...a); } catch { } },
};

const qbool = (key) => {
  try {
    const q = new URLSearchParams(window.location?.search || '');
    const v = q.get(key);
    return v === '1' || String(v).toLowerCase() === 'true';
  } catch { return false; }
};

const EDIT_MODE = qbool('edit');

// Predictable throttle (leading + trailing)
const throttle = (fn, wait) => {
  let last = 0, timer = null, lastArgs = null;
  return (...args) => {
    const now = env.now();
    const remaining = wait - (now - last);
    lastArgs = args;
    if (remaining <= 0) {
      if (timer) { clearTimeout(timer); timer = null; }
      last = now; fn(...lastArgs);
    } else if (!timer) {
      timer = setTimeout(() => { last = env.now(); timer = null; fn(...lastArgs); }, remaining);
    }
  };
};

const zOf = (e) =>
  (typeof e === 'number') ? e :
    (typeof e === 'string') ? (SYMBOL_TO_Z[e] || 0) :
      (e?.Z || e?.atomicNumber || e?.z ||
        SYMBOL_TO_Z[e?.symbol] || SYMBOL_TO_Z[e?.sym] || SYMBOL_TO_Z[e?.S] || 0);

const posToTriples = (state) => state.positions.map(p => [p.x, p.y, p.z]);
const triplesForIndices = (state, indices) =>
  indices.map((idx) => {
    const p = state.positions[idx] || { x: 0, y: 0, z: 0 };
    return [p.x, p.y, p.z];
  });

const applyTriples = (state, triples, { exclude, forceAll = false } = {}) => {
  const N = Math.min(state.positions.length, triples?.length || 0);
  if (!N) return;
  if (!forceAll && exclude?.size) {
    for (let i = 0; i < N; i++) {
      if (exclude.has(i)) continue;
      const [x, y, z] = triples[i]; const tp = state.positions[i]; tp.x = x; tp.y = y; tp.z = z;
    }
  } else {
    for (let i = 0; i < N; i++) {
      const [x, y, z] = triples[i]; const tp = state.positions[i]; tp.x = x; tp.y = y; tp.z = z;
    }
  }
};

const toSymbol = (el) => {
  const str = typeof el === 'string' ? el.trim() : '';
  if (str && SYMBOL_TO_Z[str]) return str;
  const num = Number(el);
  if (Number.isFinite(num) && Z_TO_SYMBOL[num]) return Z_TO_SYMBOL[num];
  throw new Error(`Unknown element: ${el}`);
};

const assertAllowedSymbol = (sym) => {
  if (!OMOL25_SET.has(sym)) {
    throw new Error(`Element ${sym} not supported by OMol25`);
  }
};

const toVec3 = (value, fallback = { x: 0, y: 0, z: 0 }) => {
  if (!value) return { ...fallback };
  if (Array.isArray(value)) {
    return {
      x: Number(value[0]) || 0,
      y: Number(value[1]) || 0,
      z: Number(value[2]) || 0,
    };
  }
  if (typeof value === 'object') {
    return {
      x: Number(value.x) || 0,
      y: Number(value.y) || 0,
      z: Number(value.z) || 0,
    };
  }
  return { ...fallback };
};

function clampVectorToRadius(pos, radius) {
  if (!pos || typeof pos !== 'object') return false;
  const x = Number(pos.x) || 0;
  const y = Number(pos.y) || 0;
  const z = Number(pos.z) || 0;
  const r2 = x * x + y * y + z * z;
  const limit = radius * radius;
  if (r2 === 0 || r2 <= limit) return false;
  const r = Math.sqrt(r2);
  if (!Number.isFinite(r) || r === 0) return false;
  const scale = radius / r;
  pos.x = x * scale;
  pos.y = y * scale;
  pos.z = z * scale;
  return true;
}

function enforceSafeSphere(state) {
  const radius = TUNABLES.SAFE_SPHERE_RADIUS;
  if (!Number.isFinite(radius) || radius <= 0) return false;
  const pts = state?.positions;
  if (!Array.isArray(pts) || !pts.length) return false;
  let changed = false;
  for (const p of pts) {
    if (clampVectorToRadius(p, radius)) changed = true;
  }
  if (changed && env.apiOn()) dbg.log('[safeSphere][clamp]', { radius, count: pts.length });
  return changed;
}

const stateCellToArray = (c) => (c && c.enabled) ? [
  [c.a.x, c.a.y, c.a.z],
  [c.b.x, c.b.y, c.b.z],
  [c.c.x, c.c.y, c.c.z],
] : null;

const setText = (id, txt) => { try { const el = document.getElementById(id); if (el) el.textContent = txt; } catch { } };

const onWin = (fn) => { try { if (typeof window !== 'undefined') fn(window); } catch { } };

/* ──────────────────────────────────────────────────────────────────────────
   Runtime config
   ────────────────────────────────────────────────────────────────────────── */

if (typeof window !== 'undefined') {
  window.__MLIP_CONFIG ||= {};
  if (!Number.isFinite(window.__MLIP_CONFIG.minStepIntervalMs))
    window.__MLIP_CONFIG.minStepIntervalMs = DEFAULT_MIN_STEP_INTERVAL_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.mdFriction))
    window.__MLIP_CONFIG.mdFriction = DEFAULT_MD_FRICTION;

  if (!Number.isFinite(window.__MLIP_CONFIG.dragThrottleMs))
    window.__MLIP_CONFIG.dragThrottleMs = DEFAULTS.DRAG_THROTTLE_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.posEnergyDebounceMs))
    window.__MLIP_CONFIG.posEnergyDebounceMs = DEFAULTS.POS_ENERGY_DEBOUNCE_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.latchAfterDragMs))
    window.__MLIP_CONFIG.latchAfterDragMs = DEFAULTS.LATCH_AFTER_DRAG_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.latchAfterVrMs))
    window.__MLIP_CONFIG.latchAfterVrMs = DEFAULTS.LATCH_AFTER_VR_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.rpsWindowMs))
    window.__MLIP_CONFIG.rpsWindowMs = DEFAULTS.RPS_WINDOW_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.safeSphereRadius))
    window.__MLIP_CONFIG.safeSphereRadius = DEFAULTS.SAFE_SPHERE_RADIUS;

  const autoMD = new URLSearchParams(window.location?.search || '').get('autoMD');
  if (autoMD === '0' || String(autoMD).toLowerCase() === 'false') window.__MLIPVIEW_NO_AUTO_MD = true;
}

const cfg = {
  get minStepIntervalMs() {
    return (typeof window !== 'undefined') ? window.__MLIP_CONFIG.minStepIntervalMs : DEFAULT_MIN_STEP_INTERVAL_MS;
  },
  set minStepIntervalMs(v) {
    const n = Math.max(1, Math.round(Number(v) || 0));
    if (typeof window !== 'undefined') window.__MLIP_CONFIG.minStepIntervalMs = n;
  },
  get mdFriction() {
    return (typeof window !== 'undefined') ? window.__MLIP_CONFIG.mdFriction : DEFAULT_MD_FRICTION;
  },
  set mdFriction(v) {
    const n = Math.max(0, Math.min(5, Number(v) || 0));
    if (typeof window !== 'undefined') window.__MLIP_CONFIG.mdFriction = n;
  },
};

export function setMinStepInterval(ms) {
  cfg.minStepIntervalMs = ms;
  onWin(w => dbg.log('[config] set minStepIntervalMs =', w.__MLIP_CONFIG.minStepIntervalMs));
  return cfg.minStepIntervalMs;
}
export function setMdFriction(f) {
  cfg.mdFriction = f;
  onWin(w => dbg.log('[config] set mdFriction =', w.__MLIP_CONFIG.mdFriction));
  return cfg.mdFriction;
}

/* ──────────────────────────────────────────────────────────────────────────
   Feature flags
   ────────────────────────────────────────────────────────────────────────── */

const FEATURES = (typeof window !== 'undefined' && (window.__MLIP_FEATURES ||= {
  RELAX_LOOP: true,
  MD_LOOP: true,
  ENERGY_TRACE: true,
  FORCE_VECTORS: true,
})) || { RELAX_LOOP: true, MD_LOOP: true, ENERGY_TRACE: true, FORCE_VECTORS: true };

function enableFeatureFlag(name, value = true) {
  if (typeof window !== 'undefined') window.__MLIP_FEATURES[name] = value;
  else FEATURES[name] = value;
  return true;
}

/* ──────────────────────────────────────────────────────────────────────────
   Energy plot
   ────────────────────────────────────────────────────────────────────────── */

const energyPlot = (() => {
  let series = []; // {i,E,kind}
  let last = undefined;
  let markerIndex = null;
  let ctx = null, canvas = null, label = null;

  const init = () => {
    if (canvas) return;
    canvas = document.getElementById('energyCanvas');
    if (!canvas) return;
    ctx = canvas.getContext('2d');
    label = document.getElementById('energyLabel');
  };

  const draw = () => {
    init(); if (!ctx) return;
    const W = canvas.width, H = canvas.height;
    ctx.clearRect(0, 0, W, H);
    if (series.length === 0) {
      if (label) label.textContent = 'E steps=0';
      return;
    }
    let minE = Infinity, maxE = -Infinity;
    for (const p of series) { if (p.E < minE) minE = p.E; if (p.E > maxE) maxE = p.E; }
    if (maxE - minE < 1e-12) maxE = minE + 1e-12;
    if (series.length === 1) {
      const p = series[0];
      const x = W / 2;
      const y = H - 2 - ((p.E - minE) / (maxE - minE)) * (H - 4);
      ctx.fillStyle = '#6fc2ff';
      ctx.beginPath();
      ctx.arc(x, y, 3, 0, Math.PI * 2);
      ctx.fill();
    } else {
      ctx.strokeStyle = '#6fc2ff';
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let k = 0; k < series.length; k++) {
        const p = series[k];
        const x = (k / (series.length - 1)) * (W - 4) + 2;
        const y = H - 2 - ((p.E - minE) / (maxE - minE)) * (H - 4);
        if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
    if (markerIndex != null && series.length) {
      const idx = Math.max(0, Math.min(series.length - 1, markerIndex | 0));
      const marker = series[idx];
      const x = series.length === 1 ? W / 2 : (idx / (series.length - 1)) * (W - 4) + 2;
      const y = H - 2 - ((marker.E - minE) / (maxE - minE)) * (H - 4);
      ctx.save();
      ctx.strokeStyle = '#ffd95e';
      ctx.fillStyle = '#ffd95e';
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 2]);
      ctx.beginPath();
      ctx.moveTo(x, 2);
      ctx.lineTo(x, H - 2);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.beginPath();
      ctx.arc(x, y, 3, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    }
    if (label) {
      const range = (maxE - minE);
      const rangeLabel = range > 0 ? range.toExponential(2) : '0.00e+0';
      label.textContent = `E steps=${series.length} range=${rangeLabel}`;
    }
  };

  const push = (E, kind) => {
    if (typeof E !== 'number' || !isFinite(E)) return;
    if (E === last && !(typeof window !== 'undefined' && window.__ALLOW_DUPLICATE_ENERGY_TICKS)) return;
    const i = series.length; series.push({ i, E, kind }); last = E;
    setText('instEnergy', 'E: ' + E.toFixed(2));
    draw();
    return i;
  };

  const reset = () => { series = []; last = undefined; markerIndex = null; draw(); setText('instEnergy', 'E: —'); };
  const length = () => series.length;
  const setMarker = (idx) => {
    if (!Number.isFinite(idx)) {
      markerIndex = null;
    } else {
      markerIndex = Math.max(0, Math.min(series.length - 1, idx | 0));
    }
    draw();
  };
  const clearMarker = () => { markerIndex = null; draw(); };
  const getMarker = () => {
    if (markerIndex == null || !series.length) return null;
    const idx = Math.max(0, Math.min(series.length - 1, markerIndex | 0));
    const entry = series[idx] || {};
    return { index: idx, energy: entry.E, length: series.length };
  };
  const exportSeries = () => ({
    series: series.map(({ E, kind }) => ({ energy: E, kind })),
    markerIndex: markerIndex != null ? markerIndex : null,
  });
  const importSeries = ({ series: incoming, markerIndex: marker } = {}) => {
    series = [];
    last = undefined;
    markerIndex = null;
    if (Array.isArray(incoming)) {
      for (const entry of incoming) {
        if (!entry) continue;
        const energy = Number(entry.energy);
        if (!Number.isFinite(energy)) continue;
        const kind = entry.kind && typeof entry.kind === 'string' ? entry.kind : undefined;
        const idx = series.length;
        series.push({ i: idx, E: energy, kind });
        last = energy;
      }
    }
    if (Number.isFinite(marker) && series.length) {
      markerIndex = Math.max(0, Math.min(series.length - 1, marker | 0));
    } else {
      markerIndex = null;
    }
    if (series.length) {
      const tail = series[series.length - 1];
      setText('instEnergy', 'E: ' + tail.E.toFixed(2));
    } else {
      setText('instEnergy', 'E: —');
    }
    draw();
  };

  return { push, reset, length, setMarker, clearMarker, getMarker, exportSeries, importSeries };
})();

/* ──────────────────────────────────────────────────────────────────────────
   Viewer initialization
   ────────────────────────────────────────────────────────────────────────── */

export async function initNewViewer(canvas, { elements, positions, bonds }) {
  __count('index#initNewViewer');

  onWin(w => {
    if (w.__MLIPVIEW_DEBUG_API == null) w.__MLIPVIEW_DEBUG_API = qbool('debug');
    w.__MLIPVIEW_DEBUG_PICK = qbool('debugPick');
    w.__MLIPVIEW_DEBUG_SELECT = qbool('debugSelect');
    w.__MLIPVIEW_DEBUG_UI = qbool('debugUI');
  });

  const state = createMoleculeState({ elements, positions, bonds });
  try { enforceSafeSphere(state); } catch { }
  try { state.__initialPositions = (state.positions || []).map(p => ({ x: p.x, y: p.y, z: p.z })); } catch { }

  const clonePositionList = (list) => {
    if (!Array.isArray(list)) return [];
    return list.map((p) => toVec3(p, { x: 0, y: 0, z: 0 }));
  };
  const normalizeElement = (el) => {
    if (el && typeof el === 'object') {
      if (typeof el.symbol === 'string') {
        return toSymbol(el.symbol);
      }
      if (el.Z != null) {
        return toSymbol(el.Z);
      }
    }
    return toSymbol(el);
  };
  const cloneVelocityVector = (v) => {
    if (Array.isArray(v)) {
      return [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0];
    }
    if (v && typeof v === 'object') {
      return [Number(v.x) || 0, Number(v.y) || 0, Number(v.z) || 0];
    }
    return [0, 0, 0];
  };
  const cloneVelocityList = (list) => {
    if (!Array.isArray(list) || !list.length) return null;
    return list.map((v) => cloneVelocityVector(v));
  };
  const cloneCellSnapshot = (cell) => {
    if (!cell || !cell.a || !cell.b || !cell.c) return null;
    return {
      a: toVec3(cell.a, { x: 0, y: 0, z: 0 }),
      b: toVec3(cell.b, { x: 0, y: 0, z: 0 }),
      c: toVec3(cell.c, { x: 0, y: 0, z: 0 }),
      originOffset: toVec3(cell.originOffset || { x: 0, y: 0, z: 0 }, { x: 0, y: 0, z: 0 }),
      enabled: !!cell.enabled,
    };
  };
  const cellSnapshotToMatrix = (cell) => {
    if (!cell) return null;
    return [
      [Number(cell.a?.x) || 0, Number(cell.a?.y) || 0, Number(cell.a?.z) || 0],
      [Number(cell.b?.x) || 0, Number(cell.b?.y) || 0, Number(cell.b?.z) || 0],
      [Number(cell.c?.x) || 0, Number(cell.c?.y) || 0, Number(cell.c?.z) || 0],
    ];
  };
  const captureResetBaselineFromState = () => ({
    elements: Array.isArray(state.elements) ? state.elements.map(normalizeElement) : [],
    positions: clonePositionList(state.positions),
    velocities: state.dynamics?.velocities ? cloneVelocityList(state.dynamics.velocities) : null,
    cell: cloneCellSnapshot(state.cell),
    showCell: !!state.showCell,
    showGhostCells: !!state.showGhostCells,
  });
  const installResetBaseline = (baseline, meta = {}) => {
    if (!baseline) {
      state.__resetBaseline = null;
      return;
    }
    state.__resetBaseline = {
      elements: Array.isArray(baseline.elements) ? baseline.elements.map(normalizeElement) : [],
      positions: clonePositionList(baseline.positions),
      velocities: baseline.velocities ? cloneVelocityList(baseline.velocities) : null,
      cell: cloneCellSnapshot(baseline.cell),
      showCell: !!baseline.showCell,
      showGhostCells: !!baseline.showGhostCells,
    };
    if (env.apiOn()) {
      dbg.log('[resetBaseline] updated', {
        reason: meta?.reason || 'unknown',
        natoms: state.__resetBaseline.positions?.length || 0,
      });
    }
  };
  const seedResetBaseline = (reason = 'seed') => {
    installResetBaseline(captureResetBaselineFromState(), { reason });
    try {
      state.__initialCellSnapshot = cloneCellSnapshot(state.cell);
    } catch { }
  };

  try {
    state.elements = Array.isArray(state.elements) ? state.elements.map(normalizeElement) : [];
  } catch { }
  seedResetBaseline('init');

  // Versions & caches
  state.forceCache = { version: 0, energy: NaN, forces: [], stress: null, stale: true };
  let structureVersion = 0;
  let resetEpoch = 0;
  let userInteractionVersion = 0;
  let totalInteractionVersion = 0;
  let lastAppliedUIC = 0;
  let testHoldUIC = null;
  let currentStepBudget = null;
  const frameBuffer = createFrameBuffer({ capacity: 500 });
  let timelineUi = null;
  const timelineUiRef = { value: null };
  let sessionStateManager = null;
  const getInteractionCountersSnapshot = () => ({
    user: userInteractionVersion | 0,
    total: totalInteractionVersion | 0,
    lastApplied: lastAppliedUIC | 0,
  });
  const setInteractionCountersFromSnapshot = ({ user, total, lastApplied } = {}) => {
    if (Number.isFinite(user)) userInteractionVersion = user | 0;
    if (Number.isFinite(total)) totalInteractionVersion = Math.max(userInteractionVersion, total | 0);
    else totalInteractionVersion = Math.max(userInteractionVersion, totalInteractionVersion);
    if (Number.isFinite(lastApplied)) lastAppliedUIC = lastApplied | 0;
    else lastAppliedUIC = Math.max(lastAppliedUIC, userInteractionVersion);
  };
  let timelineOverlay = null;
  const timelineState = {
    active: false,
    playing: false,
    offset: -1,
    resumeMode: null,
    suppressEnergy: false,
  };
  const playbackRuntime = {
    startOffset: null,
    loopStart: null,
    loopEnd: null,
  };
  let calloutLayer = null;
  let timelineEditorPanel = null;
  let timelineEditorStatus = null;
  let editorStatusFrame = { isLive: true };
  let editorStatusSelection = null;

  function buildSelectionSnapshot() {
    const sel = state.selection;
    if (!sel || !sel.kind) return null;
    if (sel.kind === 'atom') {
      const idx = sel.data?.index;
      if (!Number.isInteger(idx) || idx < 0 || idx >= state.positions.length) return null;
      const pos = state.positions[idx] || { x: 0, y: 0, z: 0 };
      const element = state.elements?.[idx] || null;
      return {
        kind: 'atom',
        index: idx,
        element,
        position: [pos.x ?? 0, pos.y ?? 0, pos.z ?? 0],
      };
    }
    if (sel.kind === 'bond') {
      const atoms = sel.data?.atoms || sel.data?.pair || [];
      return {
        kind: 'bond',
        atoms: Array.isArray(atoms) ? atoms.slice(0, 2) : [],
        elements: Array.isArray(atoms)
          ? atoms.slice(0, 2).map((i) => state.elements?.[i] || null)
          : [],
      };
    }
    return null;
  }

  function updateEditorStatus() {
    if (!timelineEditorStatus) return;
    timelineEditorStatus.update({
      frame: editorStatusFrame,
      selection: editorStatusSelection,
    });
  }

  function refreshEditorSelection() {
    editorStatusSelection = buildSelectionSnapshot();
    updateEditorStatus();
  }

  function refreshEditorFrame(offset) {
    if (!timelineEditorStatus) return;
    if (!timelineState.active || !Number.isFinite(offset)) {
      editorStatusFrame = { isLive: true };
      updateEditorStatus();
      return;
    }
    const entry = frameBuffer.getByOffset(offset);
    if (!entry) {
      editorStatusFrame = { isLive: true };
      updateEditorStatus();
      return;
    }
    const frameIndex = offsetToFrameIndex(offset);
    editorStatusFrame = {
      isLive: false,
      frameId: entry.id,
      offset,
      frameIndex,
    };
    updateEditorStatus();
  }

  const controlEngine = createControlMessageEngine({
    resolveFrameIndexById: (frameId) => frameBuffer.resolveFrameIndex(frameId),
    offsetToIndex: (offset) => offsetToFrameIndex(offset),
    getFrameCount: () => frameBuffer.stats().size,
  });

  const timelinePlayback = createTimelinePlaybackController({
    defaultFps: 20,
    onStep: () => timelineStepForward(),
    onPlaybackStateChange: ({ playing }) => {
      timelineState.playing = !!playing;
      if (timelineUi) timelineUi.setMode(playing ? 'playing' : 'paused');
    },
  });

  const activeControlState = {
    speed: null,
    callouts: [],
    opacity: null,
  };
  const timelineLockClass = 'timeline-readonly-overlay';

  function timelineLocked() {
    return timelineState.active;
  }
  function ensureTimelineEditable(reason) {
    if (!timelineLocked()) return true;
    if (env.apiOn()) dbg.warn('[timeline][block]', reason || '?');
    return false;
  }

  function resetStepBudgetDueToInteraction(reason) {
    if (!currentStepBudget) return;
    currentStepBudget.remaining = currentStepBudget.limit;
    currentStepBudget.overshoot = 0;
    if (env.apiOn()) dbg.log('[stepBudget][reset]', { reason, limit: currentStepBudget.limit, kind: currentStepBudget.kind });
  }

  const modifiedByVersion = new Map();
  const draggingAtoms = new Set();
  const rotatingAtoms = new Map();
  const latchedUntil = new Map();
  const bondLatch = {
    atoms: new Set(),
    active: false,
    key: null,
    anchor: null,
    movingRoot: null,
  };
  const latchAtom = (i, ms = TUNABLES.LATCH_AFTER_DRAG_MS) => {
    const now = env.now(); latchedUntil.set(i, now + ms);
  };
  const markRotatingAtoms = (indices) => {
    if (!Array.isArray(indices) || !indices.length) return;
    const now = env.now();
    const expiry = now + TUNABLES.ROTATION_LATCH_MS;
    for (const idx of indices) {
      rotatingAtoms.set(idx, expiry);
      latchedUntil.set(idx, expiry);
    }
  };
  const currentDraggedAtom = { idx: null };
  function sendBondLatchUpdate({ bumpReason = null } = {}) {
    if (!bondLatch.active || bondLatch.atoms.size === 0) return;
    if (bumpReason) {
      bumpUser(bumpReason);
      modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set());
      const modSet = modifiedByVersion.get(userInteractionVersion);
      for (const idx of bondLatch.atoms) modSet.add(idx);
    }
    const indices = Array.from(bondLatch.atoms);
    const triples = triplesForIndices(state, indices);
    try {
      const ws = getWS();
      ws.setCounters?.({ userInteractionCount: userInteractionVersion });
      ws.userInteraction?.({ positions: { indices, triples } });
    } catch { }
  }
  function clearBondLatch({ sendFinal = false, bumpReason = 'bondRotateFinalize' } = {}) {
    if (!bondLatch.active) return;
    const latchedIndices = Array.from(bondLatch.atoms);
    if (sendFinal && latchedIndices.length) {
      sendBondLatchUpdate({ bumpReason });
    }
    bondLatch.atoms.clear();
    bondLatch.active = false;
    bondLatch.key = null;
    bondLatch.anchor = null;
    bondLatch.movingRoot = null;
    for (const idx of latchedIndices) {
      rotatingAtoms.delete(idx);
      latchedUntil.delete(idx);
    }
  }
  function setBondLatchFromGroup(group) {
    if (!group || !Array.isArray(group.sideAtoms) || group.sideAtoms.length === 0) {
      clearBondLatch({ sendFinal: true });
      return;
    }
    const newKey = `${group.i}|${group.j}|${group.orientation ?? 'null'}`;
    const same =
      bondLatch.active &&
      bondLatch.key === newKey &&
      bondLatch.anchor === group.anchor &&
      bondLatch.movingRoot === group.movingRoot &&
      bondLatch.atoms.size === group.sideAtoms.length &&
      group.sideAtoms.every((idx) => bondLatch.atoms.has(idx));
    if (same) return;
    clearBondLatch({ sendFinal: true });
    bondLatch.atoms = new Set(group.sideAtoms);
    bondLatch.active = bondLatch.atoms.size > 0;
    bondLatch.key = newKey;
    bondLatch.anchor = group.anchor;
    bondLatch.movingRoot = group.movingRoot;
  }

  const buildExcludeSet = () => {
    const out = new Set(draggingAtoms);
    const now = env.now();
    for (const [i, t] of latchedUntil) { if (t > now) out.add(i); else latchedUntil.delete(i); }
    for (const [i, t] of rotatingAtoms) { if (t > now) out.add(i); else rotatingAtoms.delete(i); }
    if (bondLatch.active) {
      for (const i of bondLatch.atoms) out.add(i);
    }
    return out;
  };

  const bumpUser = (reason) => {
    userInteractionVersion++; totalInteractionVersion++; structureVersion++;
    resetStepBudgetDueToInteraction(reason);
    if (state.forceCache) { state.forceCache.stale = true; state.forceCache.version = structureVersion; }
    if (env.apiOn()) dbg.log('[version][user]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
  };
  const bumpSim = (reason) => {
    totalInteractionVersion++;
    if (env.apiOn()) dbg.log('[version][sim]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
  };

  const updateForces = (forces, { reason } = {}) => {
    const DBG = (typeof window !== 'undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search || ''));
    if (Array.isArray(forces) && forces.length) {
      state.forces = forces;
      state.bus?.emit?.('forcesChanged');
      if (DBG) dbg.log('[Forces][update]', reason || '?', 'len=', forces.length);
    }
  };

  const updateEnergyForces = ({ energy, forces, stress = null, reason }) => {
    if (typeof energy === 'number') { state.dynamics ||= {}; state.dynamics.energy = energy; }
    const effectiveForces = Array.isArray(forces) ? forces : (state.forces || []);
    if (Array.isArray(forces)) updateForces(forces, { reason });
    state.forceCache = { version: structureVersion, energy: state.dynamics?.energy, forces: effectiveForces, stress, stale: false };
    if (env.apiOn()) dbg.log('[forces-cache]', reason, { n: effectiveForces.length, E: state.forceCache.energy });
  };

  const velocitiesToArray = (vels, count) => {
    if (!Array.isArray(vels) || count <= 0) return undefined;
    const out = new Array(count);
    for (let i = 0; i < count; i++) {
      const v = vels[i];
      if (!v) {
        out[i] = [0, 0, 0];
        continue;
      }
      if (Array.isArray(v)) {
        out[i] = [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0];
      } else {
        out[i] = [Number(v.x) || 0, Number(v.y) || 0, Number(v.z) || 0];
      }
    }
    return out;
  };

  const resetInteractionCachesAfterSnapshot = () => {
    modifiedByVersion.clear();
    draggingAtoms.clear();
    rotatingAtoms.clear();
    latchedUntil.clear();
    clearBondLatch({ sendFinal: false });
  };

  const buildCellObject = (cellMatrix) => {
    if (!Array.isArray(cellMatrix) || cellMatrix.length !== 3) return null;
    const safe = cellMatrix.map((row) => (Array.isArray(row) ? row : [0, 0, 0]));
    return {
      a: { x: Number(safe[0][0]) || 0, y: Number(safe[0][1]) || 0, z: Number(safe[0][2]) || 0 },
      b: { x: Number(safe[1][0]) || 0, y: Number(safe[1][1]) || 0, z: Number(safe[1][2]) || 0 },
      c: { x: Number(safe[2][0]) || 0, y: Number(safe[2][1]) || 0, z: Number(safe[2][2]) || 0 },
      enabled: true,
      originOffset: { x: 0, y: 0, z: 0 },
    };
  };

  const updateCellFromMatrix = (matrix) => {
    if (!Array.isArray(matrix) || matrix.length !== 3) return;
    const obj = buildCellObject(matrix);
    if (!obj) return;
    state.cell ||= obj;
    state.cell.a = obj.a;
    state.cell.b = obj.b;
    state.cell.c = obj.c;
    state.cell.enabled = obj.enabled;
    state.cell.originOffset = obj.originOffset;
    state.markCellChanged?.();
  };

  const applyFullSnapshot = async ({ elements, positions, velocities, cell } = {}, opts = {}) => {
    const { updateBaseline = true } = opts || {};
    if (!ensureTimelineEditable('applyFullSnapshot')) return false;
    const nextElements = Array.isArray(elements) && elements.length
      ? elements.map(normalizeElement)
      : (state.elements || []).map(normalizeElement);
    const nextPositions = Array.isArray(positions) && positions.length
      ? positions.map((p) => toVec3(p))
      : (state.positions || []).map((p) => ({ x: p.x, y: p.y, z: p.z }));

    if (nextElements.length !== nextPositions.length) {
      throw new Error('Elements and positions length mismatch');
    }

    const natoms = nextElements.length;
    const nextVelocities = velocitiesToArray(velocities, natoms);

    state.elements = nextElements.slice();
    state.positions = nextPositions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
    state.__initialPositions = state.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
    state.forceCache && (state.forceCache.stale = true);

    state.dynamics ||= {};
    if (nextVelocities) {
      state.dynamics.velocities = nextVelocities.map((v) => [v[0], v[1], v[2]]);
    } else if (state.dynamics.velocities) {
      delete state.dynamics.velocities;
    }

    if (Array.isArray(cell) && cell.length === 3) {
      updateCellFromMatrix(cell);
      state.__initialCellSnapshot = {
        a: { ...state.cell.a },
        b: { ...state.cell.b },
        c: { ...state.cell.c },
        originOffset: { ...(state.cell.originOffset || { x: 0, y: 0, z: 0 }) },
      };
    }

    state.selection = { kind: null, data: null };
    state.markSelectionChanged?.();
    resetInteractionCachesAfterSnapshot();

    bumpUser('fullSnapshotLocal');
    state.markPositionsChanged();
    if (updateBaseline) {
      seedResetBaseline('applyFullSnapshot');
    }

    __wsState.inited = false;
    const ok = await ensureWsInit({ allowOffline: false });
    if (!ok) throw new Error('WS connection unavailable for full snapshot');
    __wsState.inited = true;
    __wsState.lastAtomCount = natoms;
    __wsState.lastCellKey = cellKey();
    return natoms;
  };

  const addAtom = async ({ element = 'H', position, velocity } = {}) => {
    if (!ensureTimelineEditable('addAtom')) return false;
    const elements = (state.elements || []).map(normalizeElement);
    const positions = (state.positions || []).map((p) => ({ x: p.x, y: p.y, z: p.z }));
    const symbol = toSymbol(element);
    assertAllowedSymbol(symbol);
    elements.push(symbol);
    positions.push(toVec3(position, { x: 0, y: 0, z: 0 }));

    let velocities = undefined;
    if (state.dynamics?.velocities) {
      velocities = state.dynamics.velocities.map((v) => {
        if (Array.isArray(v)) {
          return [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0];
        }
        return [Number(v?.x) || 0, Number(v?.y) || 0, Number(v?.z) || 0];
      });
    }
    if (velocity) {
      velocities ||= Array.from({ length: elements.length }, () => [0, 0, 0]);
      velocities[velocities.length - 1] = Array.isArray(velocity)
        ? [Number(velocity[0]) || 0, Number(velocity[1]) || 0, Number(velocity[2]) || 0]
        : [Number(velocity.x) || 0, Number(velocity.y) || 0, Number(velocity.z) || 0];
    } else if (velocities) {
      velocities.push([0, 0, 0]);
    }

    const cellArray = stateCellToArray(state.cell);
    await applyFullSnapshot({ elements, positions, velocities, cell: cellArray }, { updateBaseline: false });
  };

  const addAtomAtPosition = async (element, position, velocity) => {
    if (!ensureTimelineEditable('addAtomAtPosition')) return false;
    const symbol = toSymbol(element);
    assertAllowedSymbol(symbol);
    return addAtom({ element: symbol, position, velocity });
  };

  const addAtomAtOrigin = async (element) => addAtomAtPosition(element, { x: 0, y: 0, z: 0 });

  const removeAtoms = async (indices = []) => {
    if (!ensureTimelineEditable('removeAtoms')) return false;
    const idx = Array.from(new Set(indices.map((n) => Number(n)))).filter((n) => Number.isInteger(n) && n >= 0);
    if (!idx.length) return;
    const sorted = idx.sort((a, b) => b - a);
    const elements = (state.elements || []).map(normalizeElement);
    const positions = (state.positions || []).map((p) => ({ x: p.x, y: p.y, z: p.z }));
    let velocities = state.dynamics?.velocities
      ? state.dynamics.velocities.map((v) => {
          if (Array.isArray(v)) {
            return [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0];
          }
          return [Number(v?.x) || 0, Number(v?.y) || 0, Number(v?.z) || 0];
        })
      : undefined;

    for (const i of sorted) {
      if (i < elements.length) {
        elements.splice(i, 1);
        positions.splice(i, 1);
        velocities && velocities.splice(i, 1);
      }
    }

    const cellArray = stateCellToArray(state.cell);
    await applyFullSnapshot({ elements, positions, velocities, cell: cellArray }, { updateBaseline: false });
  };

  const removeAtomByIndex = async (index) => {
    if (!ensureTimelineEditable('removeAtomByIndex')) return false;
    if (!Number.isInteger(index) || index < 0) return false;
    const before = state.elements?.length || 0;
    await removeAtoms([index]);
    const after = state.elements?.length || 0;
    return after !== before;
  };

  const getVelocitiesIfFresh = () => {
    const v = state.dynamics?.velocities;
    if (!Array.isArray(v) || v.length !== state.elements.length) return null;
    for (let i = 0; i < v.length; i++) {
      const r = v[i]; if (!r || r.length !== 3) return null;
      if (!Number.isFinite(r[0]) || !Number.isFinite(r[1]) || !Number.isFinite(r[2])) return null;
    }
    return v;
  };

  // WS session init guard
  const __wsState = { inited: false, lastAtomCount: 0, lastCellKey: 'off' };
  const cellKey = () => {
    try {
      const c = state.cell; const on = !!(state.showCell && c && c.enabled);
      if (!on || !c) return 'off';
      return [c.a.x, c.a.y, c.a.z, c.b.x, c.b.y, c.b.z, c.c.x, c.c.y, c.c.z].map(Number).map(x => x || 0).join(',');
    } catch { return 'off'; }
  };
  const resetWsInitState = () => {
    __wsState.inited = false;
    __wsState.lastAtomCount = 0;
    __wsState.lastCellKey = 'off';
  };

  async function ensureWsInit({ allowOffline = true, timeoutMs = TUNABLES.WS_CONNECT_TIMEOUT_MS } = {}) {
    const ws = getWS();
    let connected = false;
    try {
      connected = await ws.ensureConnected({ timeoutMs });
    } catch (err) {
      resetWsInitState();
      if (!allowOffline) throw err;
      if (env.apiOn()) dbg.warn('[WS][ensureInit][connect] failed', err?.message || err);
      return false;
    }
    const wsState = (() => {
      try { return typeof ws?.getState === 'function' ? ws.getState() : null; }
      catch { return null; }
    })();
    const ready =
      !!wsState?.connected ||
      (typeof ws?.readyState === 'number' && ws.readyState === 1) ||
      (!!connected);
    if (!connected || !ws || !ready) {
      resetWsInitState();
      if (!allowOffline) throw new Error('WS not connected');
      return false;
    }
    const nat = (state.elements || []).length;
    const ckey = cellKey();
    const need = !__wsState.inited || __wsState.lastAtomCount !== nat || __wsState.lastCellKey !== ckey;
    if (need) {
      const atomic_numbers = state.elements.map(zOf);
      const positionsTriples = posToTriples(state);
      const v = getVelocitiesIfFresh();
      const cell = stateCellToArray(state.cell);
      let seq = null;
      try {
        if (env.apiOn()) {
          dbg.log('[WS][ensureInit][send]', {
            nat,
            zCount: Array.isArray(atomic_numbers) ? atomic_numbers.length : 0,
            posCount: Array.isArray(positionsTriples) ? positionsTriples.length : 0,
            cell: !!cell,
          });
        }
        ws.setCounters?.({ userInteractionCount: userInteractionVersion });
        seq = ws.userInteraction({
          natoms: nat,
          atomic_numbers,
          positions: positionsTriples,
          velocities: v || undefined,
          cell,
          full_update: true,
        });
      } catch (err) {
        resetWsInitState();
        if (!allowOffline) throw err;
        if (env.apiOn()) dbg.warn('[WS][ensureInit][send] failed', err?.message || err);
        return false;
      }
      const waitMs = Math.max(4000, Number(timeoutMs) || 0);
      if (typeof ws.waitForClientSeq === 'function' && Number.isFinite(seq) && seq > 0) {
        try {
          await ws.waitForClientSeq(seq, { timeoutMs: waitMs });
        } catch (err) {
          resetWsInitState();
          if (!allowOffline) throw err;
          if (env.apiOn()) dbg.warn('[WS][ensureInit][wait] failed', err?.message || err);
          return true;
        }
      }
      __wsState.inited = true;
      __wsState.lastAtomCount = nat;
      __wsState.lastCellKey = ckey;
      if (env.apiOn()) dbg.log('[WS][ensureInit]', { nat, v: !!v, hasCell: !!cell });
    }
    return true;
  }

  // Force provider using one-shot WS
  let lastForceResult = { energy: NaN, forces: [] };
  let inFlight = false;

  async function fetchRemoteForces({ awaitResult = false } = {}) {
    __count('index#fetchRemoteForces');
    if (inFlight && !awaitResult) return;
    while (inFlight && awaitResult) await new Promise(r => setTimeout(r, 10));
    inFlight = true;
    try {
      if (!state.elements?.length) return lastForceResult;
      const ws = getWS();
      const ok = await ensureWsInit({ allowOffline: true });
      if (!ok) return lastForceResult;
      const positions = posToTriples(state);
      try { ws.setCounters({ userInteractionCount: userInteractionVersion }); } catch { }
      if (env.apiOn()) dbg.log('[WS][simple] USER_INTERACTION', { n: positions.length, uic: userInteractionVersion });
      ws.userInteraction({ positions });
      const t0 = env.now();
      const { energy, forces } = await ws.waitForEnergy({ timeoutMs: 5000 });
      const timingMs = Math.round((env.now() - t0) * 100) / 100;
      if (env.apiOn()) dbg.log('[WS][simple][recv]', { energy, forcesLen: forces?.length || 0, timingMs });
      if (typeof energy === 'number') {
        lastForceResult = { energy, forces: forces || [] };
        updateEnergyForces({ energy, forces, reason: 'fetchRemoteForcesWS' });
        energyPlot.push(energy, 'forces');
      }
    } catch { } finally { inFlight = false; }
    return lastForceResult;
  }

  const ff = {
    computeForces: ({ sync } = {}) => {
      if (sync) return fetchRemoteForces({ awaitResult: true });
      fetchRemoteForces();
      return lastForceResult;
    },
  };

  async function requestSimpleCalculateNow() {
    try {
      const ws = getWS();
      const ok = await ensureWsInit({ allowOffline: true });
      if (!ok) return energyPlot.length();
      ws.setCounters?.({ userInteractionCount: userInteractionVersion });
      ws.userInteraction?.({ positions: posToTriples(state) });
      const { energy, forces } = await ws.waitForEnergy({ timeoutMs: 5000 });
      if (typeof energy === 'number') {
        updateEnergyForces({ energy, forces, reason: 'requestSimpleCalculateNow' });
        energyPlot.push(energy, 'forces');
      }
      return energyPlot.length();
    } catch { return energyPlot.length(); }
  }

  // Bonds
  const __origMarkPositionsChanged = state.markPositionsChanged?.bind(state) || null;
  const bondService = createBondService(state);
  const recomputeBonds = (reason = 'manual') => {
    __count('index#recomputeBonds');
    const bonds = bondService.recomputeAndStore();
    try { view.rebuildBonds(bonds); } catch { }
    if (reason !== 'markPositionsChanged') recordInteraction('rebonds');
    if (env.apiOn()) dbg.log(`[bonds] recomputed after ${reason} (count=${bonds ? bonds.length : 0})`);
    return bonds;
  };
  const __origMarkCellChanged = state.markCellChanged?.bind(state) || null;
  state.markPositionsChanged = (...a) => {
    enforceSafeSphere(state);
    structureVersion++;
    if (state.forceCache) { state.forceCache.stale = true; state.forceCache.version = structureVersion; }
    const r = __origMarkPositionsChanged ? __origMarkPositionsChanged(...a) : undefined;
    try { recomputeBonds('markPositionsChanged'); } catch { }
    return r;
  };
  state.markCellChanged = (...a) => {
    const r = __origMarkCellChanged ? __origMarkCellChanged(...a) : undefined;
    try { recomputeBonds('cellChanged'); } catch { }
    return r;
  };
  if (!state.bonds?.length) { try { recomputeBonds(); } catch { } }

  // Scene & view
  const { engine, scene, camera } = await createScene(canvas);

  const cameraControlStats = { detachCalls: 0, attachCalls: 0 };
  if (typeof window !== 'undefined' && window.__MLIPVIEW_TEST_MODE && camera && !camera.__MLIPVIEW_CONTROL_HOOKED__) {
    const origDetach = typeof camera.detachControl === 'function' ? camera.detachControl.bind(camera) : null;
    const origAttach = typeof camera.attachControl === 'function' ? camera.attachControl.bind(camera) : null;
    camera.detachControl = function (...args) {
      cameraControlStats.detachCalls += 1;
      return origDetach ? origDetach(...args) : undefined;
    };
    camera.attachControl = function (...args) {
      cameraControlStats.attachCalls += 1;
      return origAttach ? origAttach(...args) : undefined;
    };
    Object.defineProperty(camera, '__MLIPVIEW_CONTROL_HOOKED__', { value: true, configurable: true });
  }
  onWin(w => {
    w.__MLIPVIEW_DEBUG_TOUCH = !!w.__MLIPVIEW_DEBUG_TOUCH;
    w.enableTouchDebug = (on) => { w.__MLIPVIEW_DEBUG_TOUCH = !!on; dbg.log('[touchDebug] set to', w.__MLIPVIEW_DEBUG_TOUCH); };
  });
  const view = createMoleculeView(scene, state);
  const manipulation = createManipulationService(state, { bondService });

  function transformLocalToWorld(vec) {
    if (!vec) return null;
    const root = state?.moleculeRoot;
    if (root?.getWorldMatrix && typeof BABYLON !== 'undefined' && BABYLON.Vector3) {
      try {
        const v = BABYLON.Vector3.TransformCoordinates(
          new BABYLON.Vector3(vec.x || 0, vec.y || 0, vec.z || 0),
          root.getWorldMatrix()
        );
        return { x: v.x, y: v.y, z: v.z };
      } catch { }
    }
    return { x: vec.x || 0, y: vec.y || 0, z: vec.z || 0 };
  }

  function getAtomWorldPosition(idx) {
    if (!Number.isInteger(idx) || idx < 0) return null;
    const p = state.positions?.[idx];
    if (!p) return null;
    return transformLocalToWorld(p);
  }

  function getBondWorldMidpoint(i, j) {
    const a = getAtomWorldPosition(i);
    const b = getAtomWorldPosition(j);
    if (!a || !b) return null;
    return {
      x: (a.x + b.x) / 2,
      y: (a.y + b.y) / 2,
      z: (a.z + b.z) / 2,
    };
  }

  function ensureCalloutLayer() {
    if (!calloutLayer && scene) {
      calloutLayer = createCalloutLayer({
        scene,
        babylon: typeof BABYLON !== 'undefined' ? BABYLON : null,
        getAtomPositionWorld: getAtomWorldPosition,
        getBondMidpointWorld: getBondWorldMidpoint,
      });
    }
    return calloutLayer;
  }

  // Picking & touch controls
  let pickingRef = null;
  const pickingProxy = { _debug: { get dragActive() { try { return !!(pickingRef?._debug?.dragActive); } catch { return false; } } } };
  const NO_TOUCH = qbool('noTouch');
  try {
    if (NO_TOUCH) { onWin(w => { w.__MLIPVIEW_NO_TOUCH = true; dbg.warn('[touchControls] skipped by ?noTouch=1'); }); }
    else { installTouchControls({ canvas, scene, camera, picking: pickingProxy }); dbg.log('[touchControls] installed'); }
  } catch (e) { dbg.warn('[touchControls] install failed', e?.message || e); }

  // Interactions log
  const interactions = [];
  let __suppressNextPosChange = false;
  const recordInteraction = (kind) => {
    const E = state.dynamics?.energy ?? 0;
    const i = interactions.length; interactions.push({ i, kind, E });
    if (kind && kind !== 'dragMove') {
      if (!state.forces && Array.isArray(window?.__RELAX_FORCES) && window.__RELAX_FORCES.length)
        updateForces(window.__RELAX_FORCES, { reason: 'lateAttachInteraction' });
      state.bus?.emit?.('forcesChanged');
    }
  };

  // Drag emission (throttled)
  const emitDuringDrag = throttle(() => {
    try {
      enforceSafeSphere(state);
      bumpUser('dragMove');
      const ws = getWS();
      ws.setCounters?.({ userInteractionCount: userInteractionVersion });
      // Send only the currently dragged atoms (sparse delta)
      const idxs = Array.from(draggingAtoms || []);
      if (idxs.length === 0) return;
      const triples = idxs.map(i => {
        const p = state.positions[i] || { x: 0, y: 0, z: 0 };
        return [p.x, p.y, p.z];
      });
      ws.userInteraction?.({ positions: { indices: idxs, triples } });
    } catch { }
  }, TUNABLES.DRAG_THROTTLE_MS);

  // Wrapped manipulation
  function selectedAtomIndex() {
    try { return state.selection?.kind === 'atom' ? state.selection.data.index : null; }
    catch { return null; }
  }
  const wrappedManipulation = {
    beginDrag: (...a) => {
      if (!ensureTimelineEditable('manipulation.beginDrag')) return false;
      const ok = manipulation.beginDrag?.(...a);
      if (ok) {
        const idx = selectedAtomIndex(); if (idx != null) { currentDraggedAtom.idx = idx; draggingAtoms.add(idx); }
      }
      return ok;
    },
    updateDrag: (...a) => {
      if (!ensureTimelineEditable('manipulation.updateDrag')) return false;
      const r = manipulation.updateDrag?.(...a);
      if (r) {
        bumpUser('dragMove');
        const idx = selectedAtomIndex(); if (idx != null) { modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set()); modifiedByVersion.get(userInteractionVersion).add(idx); currentDraggedAtom.idx = idx; draggingAtoms.add(idx); }
        emitDuringDrag(); recordInteraction('dragMove');
      }
      return r;
    },
    endDrag: (...a) => {
      if (!ensureTimelineEditable('manipulation.endDrag')) return false;
      const had = (currentDraggedAtom.idx != null) || draggingAtoms.size > 0;
      let dragSource = null; try { dragSource = manipulation._debug?.getDragState?.()?.source || null; } catch { }
      const r = manipulation.endDrag?.(...a);
      if (had) {
        bumpUser('dragEnd');
        const idx = (currentDraggedAtom.idx != null) ? currentDraggedAtom.idx : selectedAtomIndex();
        if (idx != null) {
          modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set());
          modifiedByVersion.get(userInteractionVersion).add(idx);
          draggingAtoms.delete(idx);
          if (dragSource === 'vr') latchAtom(idx, TUNABLES.LATCH_AFTER_VR_MS);
        }
        try { const ws = getWS(); ws.setCounters?.({ userInteractionCount: userInteractionVersion }); ws.userInteraction?.({ positions: posToTriples(state) }); } catch { }
        currentDraggedAtom.idx = null;
      }
      recordInteraction('dragEnd');
      return r;
    },
    setDragPlane: (...a) => {
      if (!ensureTimelineEditable('manipulation.setDragPlane')) return;
      return manipulation.setDragPlane?.(...a);
    },
    rotateBond: (...a) => {
      if (!ensureTimelineEditable('manipulation.rotateBond')) return false;
      const r = manipulation.rotateBond?.(...a);
      if (r) {
        let sideAtoms = null; try { sideAtoms = manipulation._debug?.getLastRotation?.()?.sideAtoms || null; } catch { }
        bumpUser('bondRotate');
        if (Array.isArray(sideAtoms) && sideAtoms.length) {
          modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set());
          for (const i of sideAtoms) modifiedByVersion.get(userInteractionVersion).add(i);
          markRotatingAtoms(sideAtoms);
        }
        if (running.kind === null) ff.computeForces();
        recordInteraction('bondRotate');
        __suppressNextPosChange = true;
        sendBondLatchUpdate();
      }
      return r;
    },
  };
  try { wrappedManipulation._debug = manipulation._debug; } catch { }
  try { Object.defineProperty(wrappedManipulation, '__isWrappedManipulation', { value: true }); } catch { wrappedManipulation.__isWrappedManipulation = true; }

  // Picking with energy hook
  const selectionService = createSelectionService(state);
  try {
    state.bus?.on?.('selectionChanged', (sel) => {
      try {
        if (sel && sel.kind === 'bond') {
          const group = sel.data?.rotationGroup;
          if (group && Array.isArray(group.sideAtoms) && group.sideAtoms.length) {
            setBondLatchFromGroup(group);
            return;
          }
        }
        clearBondLatch({ sendFinal: true });
      } catch { }
      if (EDIT_MODE) {
        try { refreshEditorSelection(); } catch { }
      }
    });
  } catch { }
  const picking = (pickingRef = createPickingService(
    scene,
    view,
    selectionService,
    {
      manipulation: new Proxy({}, { get: (_, k) => wrappedManipulation?.[k] }),
      camera,
      energyHook: ({ kind }) => {
        const dragging = draggingAtoms.size > 0 || currentDraggedAtom.idx != null;
        if (mode === Mode.Idle) { if (dragging) emitDuringDrag(); else ff.computeForces(); }
        recordInteraction(kind || 'drag');
      },
      bondScrollStep: TUNABLES.BOND_SCROLL_STEP_RAD,
    },
  ));

  try {
    timelineUi = installTimeline({
      host: typeof document !== 'undefined' ? (document.getElementById('app') || canvas?.parentElement || document.body) : null,
      capacity: 500,
      getOffsets: () => frameBuffer.listOffsets(),
      getActiveOffset: () => timelineState.offset,
      onRequestOffset: (offset) => handleTimelineOffsetRequest(offset),
      onRequestPlay: (offset) => handleTimelinePlayRequest(offset),
      onRequestPause: () => handleTimelinePauseRequest(),
      onRequestLive: () => handleTimelineLiveRequest(),
    });
    timelineUiRef.value = timelineUi;
    timelineUi?.refresh?.();
  } catch (err) {
    if (env.apiOn()) dbg.warn('[timeline] install failed', err?.message || err);
  }

  if (EDIT_MODE) {
    try {
      const attachTarget = typeof document !== 'undefined' ? (document.getElementById('app') || document.body) : null;
      if (!timelineEditorStatus) {
        timelineEditorStatus = createTimelineEditorStatus({ attachTo: attachTarget });
        updateEditorStatus();
      }
      if (!timelineEditorPanel) {
        timelineEditorPanel = createTimelineEditorPanel({
          attachTo: attachTarget,
          getMessages: () => controlEngine.getSnapshot?.() || [],
          setMessages: (messages) => controlEngine.setMessages(messages),
          refreshControlEngine: () => {
            controlEngine.refresh?.();
            recomputePlaybackRuntime();
          },
          getPlaybackConfig: () => (timelinePlayback.getSnapshot ? timelinePlayback.getSnapshot() : timelinePlayback.getBaseConfig?.() || {}),
          setPlaybackConfig: (cfg) => {
            timelinePlayback.applySnapshot(cfg || {});
            recomputePlaybackRuntime();
            timelineUi?.refresh?.();
          },
          getCurrentOffset: () => timelineState.offset,
          offsetToFrameId: (offset) => {
            const entry = frameBuffer.getByOffset(offset);
            return entry?.id || null;
          },
          offsetToIndex: (offset) => offsetToFrameIndex(offset),
          getSelectionSnapshot: () => buildSelectionSnapshot(),
          getCameraTarget: () => {
            try {
              if (!camera || !camera.target) return null;
              return [camera.target.x ?? 0, camera.target.y ?? 0, camera.target.z ?? 0];
            } catch {
              return null;
            }
          },
          onCommit: () => {
            if (timelineState.active && Number.isFinite(timelineState.offset)) {
              applyTimelineFrame(timelineState.offset);
            }
            timelineUi?.refresh?.();
            try {
              sessionStateManager?.setBaselineFromState({ kind: 'json', label: 'timeline-edit' }, { includeTimeline: true });
            } catch { /* noop */ }
            refreshEditorSelection();
            refreshEditorFrame(timelineState.offset);
          },
        });
      } else {
        timelineEditorPanel.refresh();
      }
      refreshEditorSelection();
      refreshEditorFrame(timelineState.offset);
    } catch (err) {
      if (env.apiOn()) dbg.warn('[timeline-editor] init failed', err?.message || err);
    }
  }

  // VR support
  let vrPicker = null; try { vrPicker = createVRPicker({ scene, view }); } catch (e) { dbg.warn('[VR] vrPicker init failed', e?.message || e); }
  const vr = createVRSupport(scene, {
    picking: {
      ...picking,
      view,
      vrPicker,
      selectionService,
      manipulation: new Proxy({}, { get: (_, k) => wrappedManipulation?.[k] }),
      molState: state,
    },
  });
  try { vr.init().then(res => dbg.log(res.supported ? '[VR] support initialized (auto)' : '[VR] not supported')); } catch (e) { dbg.warn('[VR] auto init failed', e?.message || e); }

  // Baseline energy
  async function baselineEnergy() {
    __count('index#baselineEnergy');
    interactions.length = 0; energyPlot.reset();
    let res; try { res = await ff.computeForces({ sync: !!window?.__FORCE_SYNC_REMOTE__ }); state.dynamics ||= {}; state.dynamics.energy = res.energy; }
    catch { res = { energy: NaN, forces: [] }; }
    window.__RELAX_FORCES = res.forces || [];
    if (res.forces?.length) updateForces(res.forces, { reason: 'baselineEnergy' });
  }
  try { baselineEnergy().catch(() => { }); } catch { }

  // Seed initial positions if missing
  try {
    const need = !state.__initialPositions?.length || state.__initialPositions.length !== (state.positions?.length || 0);
    if (need && state.positions?.length) {
      state.__initialPositions = state.positions.map(p => ({ x: p.x, y: p.y, z: p.z }));
      dbg.log('[Reset] baseline seeded after baselineEnergy', { count: state.__initialPositions.length });
      if (!state.__resetBaseline || !(state.__resetBaseline.positions?.length)) {
        seedResetBaseline('baselineEnergy');
      }
    }
  } catch { }

  // Capture initial cell snapshot
  try {
    if (!state.__initialCellSnapshot && state.cell?.a && state.cell?.b && state.cell?.c) {
      const c = state.cell;
      state.__initialCellSnapshot = {
        a: { x: c.a.x, y: c.a.y, z: c.a.z },
        b: { x: c.b.x, y: c.b.y, z: c.b.z },
        c: { x: c.c.x, y: c.c.y, z: c.c.z },
        originOffset: c.originOffset ? { x: c.originOffset.x || 0, y: c.originOffset.y || 0, z: c.originOffset.z || 0 } : { x: 0, y: 0, z: 0 },
        enabled: !!c.enabled,
      };
      dbg.log('[Reset] captured initial cell snapshot');
    }
  } catch { }

  // positionsChanged debounce
  let posTickScheduled = false;
  state.bus.on('positionsChanged', () => {
    const dragging = draggingAtoms.size > 0 || currentDraggedAtom.idx != null;
    if (dragging) { emitDuringDrag(); return; }
    if (posTickScheduled) return;
    posTickScheduled = true;
    setTimeout(() => {
      posTickScheduled = false;
      if (__suppressNextPosChange) { __suppressNextPosChange = false; return; }
      if (mode === Mode.Idle) { bumpUser('posChange'); ff.computeForces(); recordInteraction('posChange'); }
    }, TUNABLES.POS_ENERGY_DEBOUNCE_MS);
  });

  // RPS counter
  const rps = { samples: [], windowMs: TUNABLES.RPS_WINDOW_MS, value: 0 };
  const resetRPS = () => { rps.samples.length = 0; rps.value = 0; setText('rpsLabel', 'RPS: --'); };
  const noteReqDone = () => {
    const now = env.now();
    rps.samples.push(now);
    while (rps.samples.length && now - rps.samples[0] > rps.windowMs) rps.samples.shift();
    rps.value = (rps.samples.length >= 2) ? ((rps.samples.length - 1)) * 1000 / (rps.samples[rps.samples.length - 1] - rps.samples[0]) : 0;
    setText('rpsLabel', 'RPS: ' + rps.value.toFixed(1));
  };

  // Focus helpers (kept for completeness; not used here)
  const isAppFocused = () => {
    try {
      if (typeof window === 'undefined' || typeof document === 'undefined') return true;
      if (window.__MLIPVIEW_TEST_MODE) return true;
      if (document.hidden) return false;
      return typeof document.hasFocus === 'function' ? !!document.hasFocus() : true;
    } catch { return true; }
  };
  const waitForFocus = () => new Promise(resolve => {
    if (isAppFocused()) return resolve();
    const cleanup = () => { try { window.removeEventListener('focus', onFocus); document.removeEventListener('visibilitychange', onVis); } catch { } };
    const onFocus = () => { cleanup(); resolve(); };
    const onVis = () => { if (!document.hidden && isAppFocused()) { cleanup(); resolve(); } };
    try { window.addEventListener('focus', onFocus); document.addEventListener('visibilitychange', onVis); } catch { resolve(); }
  });

  /* ────────────────────────────────────────────────────────────────────────
     One-shot WS steps
     ──────────────────────────────────────────────────────────────────────── */

  async function doOneStepViaWS(kind, params) {
    const ws = getWS();
    await ensureWsInit({ allowOffline: false });
    try { ws.setCounters({ userInteractionCount: userInteractionVersion, simStep: 0 }); } catch { }
    try {
      const v = getVelocitiesIfFresh();
      const prepSeq = ws.userInteraction({ positions: posToTriples(state), velocities: v || undefined });
      if (typeof ws.waitForClientSeq === 'function' && Number.isFinite(prepSeq) && prepSeq > 0) {
        const waitMs = Math.max(4000, Number(TUNABLES.WS_CONNECT_TIMEOUT_MS) || 0);
        await ws.waitForClientSeq(prepSeq, { timeoutMs: waitMs });
      }
    } catch (err) {
      if (env.apiOn()) dbg.warn(`[${kind}Step][prep] failed`, err?.message || err);
      throw err;
    }
    const r = await ws.requestSingleStep({ type: kind, params });
    const epochAtSend = resetEpoch;
    if (epochAtSend !== resetEpoch) return { stale: true, staleReason: 'staleEpoch', epochAtSend, resetEpoch };
    if (Array.isArray(r.positions) && r.positions.length === state.positions.length) {
      const exclude = buildExcludeSet();
      applyTriples(state, r.positions, { exclude });
      __suppressNextPosChange = true; state.markPositionsChanged(); bumpSim(kind === 'md' ? 'mdStepApply' : 'relaxStepApply');
    }
    updateEnergyForces({ energy: r.energy, forces: r.forces, stress: r.stress || null, reason: kind });
    if (kind === 'md') {
      if (Array.isArray(r.velocities) && r.velocities.length === state.elements.length) {
        state.dynamics ||= {}; state.dynamics.velocities = r.velocities;
      }
      if (typeof r.temperature === 'number' && isFinite(r.temperature)) {
        state.dynamics ||= {}; state.dynamics.temperature = r.temperature;
        setText('instTemp', 'T: ' + r.temperature.toFixed(1) + ' K');
      }
    }
    return { applied: true, energy: r.energy, stepType: kind, netMs: 0, parseMs: 0 };
  }

  async function relaxStep() {
    __count('index#relaxStep');
    if (!ensureTimelineEditable('relaxStep')) return { blocked: true };
    try {
      const r = await doOneStepViaWS('relax', { calculator: 'uma', fmax: state?.optimizer?.fmax || 0.05, max_step: MAX_STEP, optimizer: 'bfgs' });
      if (r?.applied) { recordInteraction('relaxStep'); energyPlot.push(state.dynamics?.energy, 'relax'); }
      return r;
    } catch (e) { dbg.warn('[relaxStep] failed', e); return { error: e?.message || String(e) }; }
  }

  async function mdStep(opts = {}) {
    __count('index#mdStep');
    if (!ensureTimelineEditable('mdStep')) return { blocked: true };
    try {
      let temperature = opts.temperature;
      try { if ((temperature == null || !Number.isFinite(temperature)) && typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null) temperature = Number(window.__MLIP_TARGET_TEMPERATURE); } catch { }
      const r = await doOneStepViaWS('md', {
        calculator: opts.calculator || 'uma',
        temperature,
        timestep_fs: opts.timestep_fs,
        friction: (Number.isFinite(opts.friction) ? opts.friction : cfg.mdFriction),
      });
      if (r?.applied) { recordInteraction('mdStep'); energyPlot.push(state.dynamics?.energy, 'md'); }
      return r;
    } catch (e) { dbg.warn('[mdStep] failed', e); return { error: e?.message || String(e) }; }
  }

  /* ────────────────────────────────────────────────────────────────────────
     Streaming WS loops (unified)
     ──────────────────────────────────────────────────────────────────────── */

  function applyFramePayload(kind, frame, {
    ws = null,
    isLive = false,
    allowEnergyPlot = true,
    forceAll = false,
    advanceCounters = false,
    noteFrame = false,
    updateUIC = false,
  } = {}) {
    if (!frame) return { applied: false };
    const uic = Number(frame.userInteractionCount ?? frame.user_interaction_count ?? 0) | 0;
    const seq = Number(frame.seq) || 0;
    if (isLive && ws && typeof ws.ack === 'function' && seq > 0) {
      try { ws.ack(seq); } catch { }
    }
    const inTestMode = (typeof window !== 'undefined') && !!window.__MLIPVIEW_TEST_MODE;
    const isTestInjected = frame.__testInjected === true;
    if (inTestMode) {
      if (isTestInjected) {
        testHoldUIC = uic;
      } else if (testHoldUIC != null) {
        if (uic <= testHoldUIC) {
          if (env.apiOn()) dbg.log(`[${kind}WS][discard-test-hold]`, { server: uic, hold: testHoldUIC });
          return { applied: false, skipped: 'test-hold', uic };
        }
        testHoldUIC = null;
      }
    }
    if (isLive && uic < lastAppliedUIC) {
      if (env.apiOn()) dbg.log(`[${kind}WS][discard-stale-applied]`, { server: uic, lastAppliedUIC });
      return { applied: false, skipped: 'stale', uic };
    }

    const exclude = forceAll ? null : buildExcludeSet();
    let applied = false;
    let energyPlotIndex = null;
    if (Array.isArray(frame.positions) && frame.positions.length === state.positions.length) {
      applyTriples(state, frame.positions, { exclude, forceAll });
      __suppressNextPosChange = true;
      state.markPositionsChanged();
      if (advanceCounters) {
        bumpSim(kind === 'md' ? 'mdWS' : (kind === 'relax' ? 'relaxWS' : 'idleWS'));
      }
      applied = true;
    }

    const energyValue = (typeof frame.energy === 'number') ? frame.energy : state.dynamics?.energy;
    updateEnergyForces({
      energy: energyValue,
      forces: Array.isArray(frame.forces) && frame.forces.length ? frame.forces : undefined,
      stress: frame.stress || null,
      reason: isLive ? kind + 'WS' : kind + 'Replay',
    });

    if (kind === 'md') {
      if (Array.isArray(frame.velocities) && frame.velocities.length === state.elements.length) {
        state.dynamics ||= {};
        state.dynamics.velocities = frame.velocities;
      }
      if (typeof frame.temperature === 'number' && Number.isFinite(frame.temperature)) {
        state.dynamics ||= {};
        state.dynamics.temperature = frame.temperature;
        setText('instTemp', 'T: ' + frame.temperature.toFixed(1) + ' K');
      } else if (!isLive && typeof state.dynamics?.temperature === 'number') {
        setText('instTemp', 'T: ' + state.dynamics.temperature.toFixed(1) + ' K');
      }
      try {
        const tm = typeof window !== 'undefined' && !!window.__MLIPVIEW_TEST_MODE;
        const haveV = Array.isArray(state.dynamics?.velocities) && state.dynamics.velocities.length === state.elements.length;
        if (tm && !haveV && state.positions?.length === state.elements.length) {
          state.dynamics ||= {};
          state.dynamics.velocities = state.positions.map(() => [0, 0, 0]);
        }
      } catch { }
    }

    if (allowEnergyPlot && typeof energyValue === 'number' && Number.isFinite(energyValue)) {
      const pushedIndex = energyPlot.push(energyValue, kind);
      if (Number.isFinite(pushedIndex)) {
        energyPlotIndex = pushedIndex;
      } else if (energyPlot.length() > 0) {
        energyPlotIndex = energyPlot.length() - 1;
      }
    } else if (typeof energyValue === 'number' && Number.isFinite(energyValue)) {
      setText('instEnergy', 'E: ' + energyValue.toFixed(2));
    }

    if (noteFrame) noteReqDone();
    if (updateUIC && isLive && uic > lastAppliedUIC) lastAppliedUIC = uic;

    return { applied: applied || Number.isFinite(energyValue), energy: energyValue, uic, energyIndex: energyPlotIndex };
  }

  function handleStreamFrame(kind, ws, frame) {
    const allowEnergy = !timelineState.suppressEnergy;
    const result = applyFramePayload(kind, frame, {
      ws,
      isLive: true,
      allowEnergyPlot: allowEnergy,
      forceAll: false,
      advanceCounters: true,
      noteFrame: true,
      updateUIC: true,
    });
    if (result?.applied) {
      frameBuffer.record(kind, frame, { energyIndex: result.energyIndex });
      recomputePlaybackRuntime();
      if (timelineUi) timelineUi.refresh();
    }
  }

  function ensureTimelineOverlay() {
    if (typeof document === 'undefined') return;
    if (timelineOverlay) return;
    const host = document.getElementById('app') || canvas?.parentElement || document.body;
    if (!host) return;
    const overlay = document.createElement('div');
    overlay.className = timelineLockClass;
    overlay.dataset.testid = 'timeline-overlay';
    overlay.style.position = 'absolute';
    overlay.style.left = '0';
    overlay.style.right = '0';
    overlay.style.top = '0';
    overlay.style.bottom = '0';
    overlay.style.pointerEvents = 'none';
    overlay.style.background = 'rgba(0, 0, 0, 0.02)';
    overlay.style.display = 'none';
    overlay.style.zIndex = '26';
    overlay.innerHTML = '<div style="position:absolute;bottom:26px;right:26px;padding:6px 12px;border-radius:8px;background:rgba(10,16,24,0.9);color:#e4ecff;font:12px system-ui,sans-serif;pointer-events:none;">Timeline view — read only</div>';
    host.appendChild(overlay);
    timelineOverlay = overlay;
  }

  function setTimelineOverlayVisible(on) {
    ensureTimelineOverlay();
    if (!timelineOverlay) return;
    timelineOverlay.style.display = on ? 'block' : 'none';
  }

  function setTimelineInteractionLock(on) {
    const allow = !on;
    try { manipulation?.setInteractionsEnabled?.(allow); } catch { }
    try {
      // Keep camera navigation responsive while timeline is active.
      pickingRef?.setInteractionsEnabled?.(true);
    } catch { }
    try { vr?.setInteractionsEnabled?.(allow); } catch { }
  }

  function clearTimelinePlayback() {
    timelinePlayback.stop();
    timelineState.playing = false;
  }

  function clampOffsetToBuffer(offset) {
    const offsets = frameBuffer.listOffsets();
    if (!offsets.length) return null;
    const minOffset = offsets[offsets.length - 1];
    const raw = Number.isFinite(offset) ? Math.floor(offset) : -1;
    if (raw > -1) return -1;
    if (raw < minOffset) return minOffset;
    return raw;
  }

  function getFrameCount() {
    const stats = frameBuffer.stats();
    return stats?.size || 0;
  }

  function offsetToFrameIndex(offset) {
    const size = getFrameCount();
    if (!size) return null;
    const off = Math.floor(Number(offset));
    if (!Number.isFinite(off)) return null;
    if (off >= -1) return size - 1;
    let idx = size + off;
    if (idx < 0) idx = 0;
    if (idx >= size) idx = size - 1;
    return idx;
  }

  function frameIndexToOffset(index) {
    const size = getFrameCount();
    if (!size) return null;
    const idx = Math.max(0, Math.min(size - 1, Math.floor(Number(index))));
    return idx - size;
  }

  function resolveFrameReference(ref) {
    if (!ref || typeof ref !== 'object') return null;
    if (typeof ref.frameId === 'string' && ref.frameId) {
      const off = frameBuffer.resolveOffset(ref.frameId);
      if (off != null) return clampOffsetToBuffer(off);
    }
    if (Number.isFinite(ref.frameIndex)) {
      const off = frameIndexToOffset(ref.frameIndex);
      return off != null ? clampOffsetToBuffer(off) : null;
    }
    if (Number.isFinite(ref.offset)) {
      return clampOffsetToBuffer(ref.offset);
    }
    return null;
  }

  function recomputePlaybackRuntime() {
    const size = getFrameCount();
    controlEngine.updateFrameCount(size);
    if (!size) {
      playbackRuntime.startOffset = null;
      playbackRuntime.loopStart = null;
      playbackRuntime.loopEnd = null;
      return;
    }
    const offsets = frameBuffer.listOffsets();
    const minOffset = offsets[offsets.length - 1];
    const maxOffset = -1;
    const baseConfig = timelinePlayback.getBaseConfig ? timelinePlayback.getBaseConfig() : {};
    const startOffset = resolveFrameReference(baseConfig.startFrame) ?? maxOffset;
    playbackRuntime.startOffset = clampOffsetToBuffer(startOffset ?? maxOffset);
    if (baseConfig.loop) {
      const range = baseConfig.loopRange || {};
      const startRef = range.start || (range.startFrameId ? { frameId: range.startFrameId } : null);
      const endRef = range.end || (range.endFrameId ? { frameId: range.endFrameId } : null);
      let loopStart = resolveFrameReference(startRef);
      let loopEnd = resolveFrameReference(endRef);
      if (loopStart == null) loopStart = minOffset;
      if (loopEnd == null) loopEnd = maxOffset;
      if (loopStart > loopEnd) [loopStart, loopEnd] = [loopEnd, loopStart];
      playbackRuntime.loopStart = clampOffsetToBuffer(loopStart);
      playbackRuntime.loopEnd = clampOffsetToBuffer(loopEnd);
    } else {
      playbackRuntime.loopStart = null;
      playbackRuntime.loopEnd = null;
    }
  }

  function applyTimelineFrame(offset) {
    const entry = frameBuffer.getByOffset(offset);
    if (!entry) return false;
    timelineState.offset = offset;
    const payloadResult = applyFramePayload(entry.kind || 'idle', entry, {
      isLive: false,
      allowEnergyPlot: false,
      forceAll: true,
      advanceCounters: false,
      noteFrame: false,
      updateUIC: false,
    });
    if (typeof entry.energyIndex === 'number') {
      energyPlot.setMarker(entry.energyIndex);
    } else if (typeof payloadResult?.energyIndex === 'number') {
      energyPlot.setMarker(payloadResult.energyIndex);
    } else {
      energyPlot.clearMarker();
    }
    const frameIndex = offsetToFrameIndex(offset);
    const controlResult = Number.isInteger(frameIndex)
      ? controlEngine.evaluate(frameIndex)
      : { speed: null, callout: null, callouts: [], opacity: null, activeMessages: [] };
    applyControlActions(controlResult, { frame: entry, offset, frameIndex });
    if (timelineUi) timelineUi.setActiveOffset(offset);
    if (EDIT_MODE) {
      try { refreshEditorFrame(offset); } catch { }
    }
    return true;
  }

  function applyControlActions(result, meta) {
    const speed = result?.speed || null;
    if (speed) {
      const prev = activeControlState.speed;
      if (!prev || prev.fps !== speed.fps || prev.speedMultiplier !== speed.speedMultiplier || prev.sourceId !== speed.sourceId) {
        timelinePlayback.setSpeedOverride(speed);
        activeControlState.speed = {
          sourceId: speed.sourceId || null,
          fps: speed.fps ?? null,
          speedMultiplier: speed.speedMultiplier ?? null,
        };
      }
    } else if (activeControlState.speed) {
      timelinePlayback.clearOverrides();
      activeControlState.speed = null;
    }

    const callouts = Array.isArray(result?.callouts) && result.callouts.length
      ? result.callouts
      : result?.callout
        ? [result.callout]
        : [];
    if (callouts.length) {
      const layer = ensureCalloutLayer();
      layer?.showAll?.(callouts);
      layer?.update?.();
      activeControlState.callouts = callouts.map((c) => ({ key: c.key, sourceId: c.sourceId }));
    } else {
      if (activeControlState.callouts?.length) calloutLayer?.hide();
      else calloutLayer?.update?.();
      activeControlState.callouts = [];
    }

    const opacity = result?.opacity || null;
    if (opacity) {
      try { view?.setOpacityMask?.(opacity); } catch { }
      activeControlState.opacity = opacity;
    } else if (activeControlState.opacity) {
      try { view?.setOpacityMask?.(null); } catch { }
      activeControlState.opacity = null;
    }
  }

  function enterTimelineMode(offset = -1) {
    const offsets = frameBuffer.listOffsets();
    if (!offsets.length) return false;
    const minOffset = offsets[offsets.length - 1];
    const target = Math.min(-1, Math.max(minOffset, Number.isFinite(offset) ? Math.floor(offset) : -1));
    if (!timelineState.active) {
      timelineState.resumeMode = mode;
      timelineState.active = true;
      timelineState.suppressEnergy = true;
      clearTimelinePlayback();
      setTimelineOverlayVisible(true);
      setTimelineInteractionLock(true);
      try { getWS().stopSimulation(); } catch { }
      setMode(Mode.Timeline);
    } else {
      clearTimelinePlayback();
    }
    applyTimelineFrame(target);
    if (timelineUi) {
      timelineUi.setMode('paused');
      timelineUi.refresh();
    }
    return true;
  }

  function timelineStepForward() {
    const offsets = frameBuffer.listOffsets();
    if (!offsets.length) {
      timelinePlayback.stop();
      if (timelineUi) timelineUi.setMode('paused');
      return false;
    }
    const minOffset = offsets[offsets.length - 1];
    const maxOffset = -1;
    let currentOffset = timelineState.offset;
    if (!Number.isFinite(currentOffset)) {
      currentOffset = playbackRuntime.startOffset != null ? playbackRuntime.startOffset : maxOffset;
    }
    if (currentOffset < minOffset) currentOffset = minOffset;
    const loopEnabled =
      timelinePlayback.getLoopConfig &&
      timelinePlayback.getLoopConfig().loop &&
      Number.isFinite(playbackRuntime.loopStart) &&
      Number.isFinite(playbackRuntime.loopEnd);
    let next = currentOffset >= maxOffset ? maxOffset : currentOffset + 1;
    if (loopEnabled && playbackRuntime.loopStart != null && playbackRuntime.loopEnd != null) {
      if (next > playbackRuntime.loopEnd) {
        next = playbackRuntime.loopStart;
      }
    } else if (next > maxOffset) {
      timelinePlayback.stop();
      resumeLiveFromTimeline().catch(() => { });
      return false;
    }
    next = clampOffsetToBuffer(next);
    if (!Number.isFinite(next)) return false;
    if (!applyTimelineFrame(next)) {
      timelinePlayback.stop();
      if (timelineUi) timelineUi.setMode('paused');
      return false;
    }
    if (!loopEnabled && next === maxOffset) {
      timelinePlayback.stop();
      resumeLiveFromTimeline().catch(() => { });
      return false;
    }
    return true;
  }

  function pauseTimelinePlayback() {
    clearTimelinePlayback();
    if (timelineUi) timelineUi.setMode('paused');
  }

  async function resumeLiveFromTimeline() {
    if (!timelineState.active) return;
    try {
      console.log('[Timeline][resumeLiveFromTimeline] begin', {
        resumeMode: timelineState.resumeMode,
        lastOffsets: frameBuffer?.stats?.() || null,
      });
    } catch { }
    clearTimelinePlayback();
    timelineState.active = false;
    timelineState.suppressEnergy = false;
    setTimelineOverlayVisible(false);
    setTimelineInteractionLock(false);
    energyPlot.clearMarker();
    timelineState.offset = -1;
    if (EDIT_MODE) {
      editorStatusFrame = { isLive: true };
      updateEditorStatus();
    }
    timelinePlayback.clearOverrides();
    activeControlState.speed = null;
    if (view?.setOpacityMask) {
      try { view.setOpacityMask(null); } catch { }
    }
    activeControlState.opacity = null;
    if (calloutLayer) {
      try { calloutLayer.hide(); } catch { }
    }
    activeControlState.callouts = [];
    const resumeMode = timelineState.resumeMode;
    timelineState.resumeMode = null;
    if (timelineUi) {
      timelineUi.setMode('live');
      timelineUi.setActiveOffset(-1);
    }
    setMode(Mode.Idle);
    if (timelineUi) timelineUi.refresh();
    try {
      await ensureWsInit({ allowOffline: false });
      try {
        console.log('[Timeline][resumeLiveFromTimeline] ensureWsInit done', {
          resumeMode,
          wsState: getWS()?.getState?.() || null,
        });
      } catch { }
    } catch (err) {
      if (env.apiOn()) dbg.warn('[timeline][resume] ensureWsInit failed', err?.message || err);
    }
    if (resumeMode === Mode.MD) {
      const opts = lastContinuousOpts.md ? { ...lastContinuousOpts.md } : {};
      startContinuous('md', opts).catch(() => { });
      try {
        console.log('[Timeline][resumeLiveFromTimeline] restarting MD', opts);
      } catch { }
    } else if (resumeMode === Mode.Relax) {
      const opts = lastContinuousOpts.relax ? { ...lastContinuousOpts.relax } : {};
      startContinuous('relax', opts).catch(() => { });
      try {
        console.log('[Timeline][resumeLiveFromTimeline] restarting Relax', opts);
      } catch { }
    }
  }

  function handleTimelineOffsetRequest(offset) {
    const clamped = clampOffsetToBuffer(offset);
    if (clamped == null) return;
    if (!timelineState.active) {
      const entered = enterTimelineMode(clamped);
      if (!entered) return;
      if (timelineState.offset !== clamped) {
        pauseTimelinePlayback();
        applyTimelineFrame(clamped);
      }
      return;
    }
    pauseTimelinePlayback();
    applyTimelineFrame(clamped);
    if (timelineUi) timelineUi.setMode('paused');
  }

  function handleTimelinePlayRequest(offset) {
    const clamped = clampOffsetToBuffer(Number.isFinite(offset) ? offset : timelineState.offset);
    if (clamped == null) return;
    enterTimelineMode(clamped);
    timelinePlayback.start();
  }

  function handleTimelinePauseRequest() {
    if (!timelineState.active) {
      enterTimelineMode(-1);
      return;
    }
    pauseTimelinePlayback();
  }

  function handleTimelineLiveRequest() {
    resumeLiveFromTimeline().catch(() => { });
  }

  function clearReconnectCountdown() {
    if (reconnectState.countdownTimer) {
      clearInterval(reconnectState.countdownTimer);
      reconnectState.countdownTimer = null;
    }
  }

  function hideReconnectBannerUi() {
    if (!reconnectState.bannerVisible) return;
    reconnectState.bannerVisible = false;
    clearReconnectCountdown();
    try { hideErrorBanner(); } catch { }
  }

  function updateReconnectCountdownText() {
    if (!reconnectState.bannerVisible) return;
    const label = typeof document !== 'undefined' ? document.getElementById('wsReconnectCountdown') : null;
    if (!label) return;
    const remainingMs = Math.max(0, (reconnectState.nextAttemptAt || Date.now()) - Date.now());
    const seconds = Math.max(0, Math.ceil(remainingMs / 1000));
    label.textContent = `Reconnecting in ${seconds}s (attempt ${Math.max(0, reconnectState.attempts)})`;
  }

  function showReconnectBannerUi(data) {
    reconnectState.attempts = Math.max(0, data?.attempts || reconnectState.attempts || 0);
    reconnectState.nextAttemptAt = data?.nextAttemptAt || Date.now();
    if (reconnectState.attempts < 1) {
      hideReconnectBannerUi();
      return;
    }
    reconnectState.bannerVisible = true;
    showErrorBanner('', {
      persist: true,
      className: 'ws-reconnect-banner',
      onRender: (el) => {
        try {
          el.innerHTML = `
            <span class="ws-reconnect-message">
              Connection to the simulation server was lost.
              <strong id="wsReconnectCountdown"></strong>
            </span>
            <button id="wsReconnectNow" class="ws-reconnect-button">Reconnect now</button>
          `;
          el.style.background = '#f8c9c9';
          el.style.color = '#5a0000';
          el.style.display = 'flex';
          el.style.flexWrap = 'wrap';
          el.style.justifyContent = 'center';
          el.style.alignItems = 'center';
          el.style.gap = '12px';
          el.dataset.banner = 'ws-reconnect';
          const btn = el.querySelector('#wsReconnectNow');
          if (btn) {
            btn.style.background = '#ffffff';
            btn.style.border = '1px solid #bb3d3d';
            btn.style.color = '#5a0000';
            btn.style.padding = '4px 10px';
            btn.style.cursor = 'pointer';
            btn.onclick = () => {
              try { getWS().reconnectNow(); } catch { }
            };
          }
        } catch { }
        updateReconnectCountdownText();
      },
    });
    clearReconnectCountdown();
    reconnectState.countdownTimer = setInterval(updateReconnectCountdownText, 500);
  }

  const Mode = Object.freeze({ Idle: 'idle', MD: 'md', Relax: 'relax', Timeline: 'timeline' });

  function rememberResume(kind, opts) {
    if (!kind) return;
    pendingResume.kind = kind;
    pendingResume.opts = opts ? { ...opts } : {};
  }

  function consumePendingResume() {
    const info = { kind: pendingResume.kind, opts: pendingResume.opts ? { ...pendingResume.opts } : {} };
    pendingResume.kind = null;
    pendingResume.opts = null;
    return info;
  }

  addWsStateListener((evt) => {
    const type = evt?.type;
    if (type === 'reconnect-scheduled') {
      showReconnectBannerUi(evt || {});
    } else if (type === 'open') {
      hideReconnectBannerUi();
      reconnectState.attempts = 0;
      reconnectState.nextAttemptAt = 0;
      __wsState.inited = false;
      __wsState.lastAtomCount = 0;
      __wsState.lastCellKey = 'off';
      const pending = consumePendingResume();
      Promise.resolve()
        .then(() => ensureWsInit({ allowOffline: true }))
        .then(() => {
          attachIdleWSListener().catch(() => { });
          if (pending.kind === 'md') {
            startContinuous('md', pending.opts || {}).catch(() => { });
          } else if (pending.kind === 'relax') {
            startContinuous('relax', pending.opts || {}).catch(() => { });
          }
        })
        .catch(() => { });
    } else if (type === 'close') {
      __wsState.inited = false;
      __wsState.lastAtomCount = 0;
      __wsState.lastCellKey = 'off';
      const currentMode = mode;
      if (currentMode === Mode.MD || currentMode === Mode.Relax) {
        const kind = currentMode === Mode.MD ? 'md' : 'relax';
        rememberResume(kind, lastContinuousOpts[kind]);
      }
      try { clearMdStreamListener(); } catch { }
      try { clearRelaxStreamListener(); } catch { }
      setMode(Mode.Idle);
    }
  });
  let mode = Mode.Idle;
  let running = { kind: null, abort: null, stepBudget: null }; // keep for API compatibility

  let idleUnsub = null;
  let mdStreamUnsub = null;
  let relaxStreamUnsub = null;

  const streamListenerStats = {
    md: { attach: 0, detach: 0, active: 0, maxActive: 0 },
    relax: { attach: 0, detach: 0, active: 0, maxActive: 0 },
  };

  const statsFor = (kind) => (kind === 'relax' ? streamListenerStats.relax : streamListenerStats.md);

  function resetStreamListenerStats() {
    for (const key of Object.keys(streamListenerStats)) {
      const stat = streamListenerStats[key];
      stat.attach = 0;
      stat.detach = 0;
      stat.active = 0;
      stat.maxActive = 0;
    }
  }

  function trackAttach(kind) {
    const stat = statsFor(kind);
    stat.attach += 1;
    stat.active += 1;
    if (stat.active > stat.maxActive) stat.maxActive = stat.active;
  }

  function trackDetach(kind) {
    const stat = statsFor(kind);
    stat.detach += 1;
    stat.active = Math.max(0, stat.active - 1);
  }

  function makeTrackedUnsub(kind, unsub) {
    let done = false;
    const tracked = () => {
      if (done) return;
      done = true;
      try { unsub && unsub(); } catch { }
      trackDetach(kind);
      if (kind === 'md' && mdStreamUnsub === tracked) mdStreamUnsub = null;
      if (kind === 'relax' && relaxStreamUnsub === tracked) relaxStreamUnsub = null;
    };
    return tracked;
  }

  function clearMdStreamListener() {
    if (typeof mdStreamUnsub === 'function') {
      const fn = mdStreamUnsub;
      mdStreamUnsub = null;
      fn();
    }
  }

  function clearRelaxStreamListener() {
    if (typeof relaxStreamUnsub === 'function') {
      const fn = relaxStreamUnsub;
      relaxStreamUnsub = null;
      fn();
    }
  }

  function debugStreamListenerStats({ reset = false } = {}) {
    if (reset) resetStreamListenerStats();
    return {
      md: { ...streamListenerStats.md },
      relax: { ...streamListenerStats.relax },
    };
  }
  async function attachIdleWSListener() {
    const ws = getWS();
    if (!idleUnsub) {
      idleUnsub = ws.onResult((r) => {
        if (mode !== Mode.Idle) return;
        handleStreamFrame('idle', ws, r);
      });
    }
    const ok = await ensureWsInit({ allowOffline: true });
    return ok;
  }
  function detachIdleWSListener() { try { idleUnsub && idleUnsub(); } catch { } finally { idleUnsub = null; } }

  function setMode(next) {
    if (mode === next) return;
    if (mode === Mode.MD && next !== Mode.MD) clearMdStreamListener();
    if (mode === Mode.Relax && next !== Mode.Relax) clearRelaxStreamListener();
    mode = next;
    if (mode !== Mode.MD && mode !== Mode.Timeline) {
      try { if (state?.dynamics) state.dynamics.temperature = undefined; } catch { }
      setText('instTemp', TEMP_PLACEHOLDER);
    }
    if (mode === Mode.Idle) {
      running.kind = null;
      running.stepBudget = null;
      currentStepBudget = null;
      timelineState.suppressEnergy = false;
      attachIdleWSListener().catch(() => { });
    } else if (mode === Mode.Timeline) {
      running.kind = null;
      running.stepBudget = null;
      currentStepBudget = null;
      detachIdleWSListener();
    } else {
      running.kind = (mode === Mode.MD) ? 'md' : 'relax';
      detachIdleWSListener();
    }
  }

  async function startContinuous(kind, opts = {}) {
    const isMD = kind === 'md';
    if ((isMD && !FEATURES.MD_LOOP) || (!isMD && !FEATURES.RELAX_LOOP)) { dbg.warn('[feature] loop disabled', kind); return { disabled: true }; }
    if (mode !== Mode.Idle) return { ignored: true };
    setMode(isMD ? Mode.MD : Mode.Relax);
    lastContinuousOpts[isMD ? 'md' : 'relax'] = { ...(opts || {}) };
    if (pendingResume.kind) {
      pendingResume.kind = null;
      pendingResume.opts = null;
    }

    const ws = getWS();
    try { await ws.ensureConnected({ timeoutMs: TUNABLES.WS_CONNECT_TIMEOUT_MS }); } catch { }
    let initOk = false;
    let lastErr = null;
    for (let attempt = 0; attempt < 2 && !initOk; attempt++) {
      try {
        await ensureWsInit({ allowOffline: false });
        initOk = true;
      } catch (err) {
        lastErr = err;
        if (attempt === 0) {
          try { await ws.ensureConnected({ timeoutMs: TUNABLES.WS_CONNECT_TIMEOUT_MS }); } catch { }
          await new Promise(resolve => setTimeout(resolve, 100));
          continue;
        }
      }
    }
    if (!initOk) {
      if (env.apiOn()) dbg.warn(`[${kind}WS][connect] proceeding without confirmed init`, lastErr?.message || lastErr);
    }
    try { if (posTickScheduled) posTickScheduled = false; } catch { }
    try { ws.setCounters({ userInteractionCount: userInteractionVersion, simStep: 0 }); } catch { }
    let prepSeq = null;
    try {
      const v = getVelocitiesIfFresh();
      prepSeq = ws.userInteraction({ positions: posToTriples(state), velocities: v || undefined });
      if (typeof ws.waitForClientSeq === 'function' && Number.isFinite(prepSeq) && prepSeq > 0) {
        const waitMs = Math.max(4000, Number(TUNABLES.WS_CONNECT_TIMEOUT_MS) || 0);
        await ws.waitForClientSeq(prepSeq, { timeoutMs: waitMs });
      }
    } catch (err) {
      if (env.apiOn()) dbg.warn(`[${kind}WS][prep] failed`, err?.message || err);
      if (typeof prepSeq === 'number' && prepSeq > 0) {
        try { ws.ack?.(prepSeq); } catch { }
      }
    }
    resetRPS();

    const rawLimit = isMD ? opts.steps : opts.maxSteps;
    const limitCandidate = Number.isFinite(rawLimit) ? Math.floor(rawLimit) : DEFAULTS.CONTINUOUS_STEPS;
    const stepLimit = Math.max(1, limitCandidate);
    const extraAllowance = isMD ? 5 : 0;

    currentStepBudget = {
      kind: isMD ? 'md' : 'relax',
      limit: stepLimit,
      remaining: stepLimit,
      extraAllowance,
      overshoot: 0,
    };
    running.stepBudget = currentStepBudget;

    let unsub = null, stopped = false;
    if (isMD) clearMdStreamListener();
    else clearRelaxStreamListener();
    unsub = ws.onResult((r) => {
      if (stopped || mode !== (isMD ? Mode.MD : Mode.Relax)) return;
      if (env.apiOn()) dbg.log(`[${kind}WS][frame]`, { seq: Number(r.seq) || 0, uic: Number(r.userInteractionCount) || 0, simStep: Number(r.simStep) || 0, have: { pos: !!r.positions, forces: !!r.forces, energy: typeof r.energy === 'number' } });
      handleStreamFrame(kind, ws, r);
      const budget = currentStepBudget;
      if (!budget || budget.kind !== (isMD ? 'md' : 'relax')) return;
      if (budget.remaining > 0) {
        budget.remaining -= 1;
        budget.overshoot = 0;
        return;
      }

      if (isMD) {
        const haveV = Array.isArray(state.dynamics?.velocities) && state.dynamics.velocities.length === state.elements.length;
        if (haveV || budget.overshoot >= budget.extraAllowance) {
          stopped = true; try { ws.stopSimulation(); } catch { }
          clearMdStreamListener();
          setMode(Mode.Idle); resetRPS();
        } else {
          budget.overshoot += 1;
          return;
        }
      } else {
        stopped = true; try { ws.stopSimulation(); } catch { }
        clearRelaxStreamListener();
        setMode(Mode.Idle); resetRPS();
      }
    });
    if (isMD) {
      trackAttach('md');
      mdStreamUnsub = makeTrackedUnsub('md', unsub);
    } else {
      trackAttach('relax');
      relaxStreamUnsub = makeTrackedUnsub('relax', unsub);
    }

    // live temperature from global target if present (MD)
    let temperature = isMD ? (opts.temperature ?? 1500) : undefined;
    try { if (isMD && typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null) temperature = Number(window.__MLIP_TARGET_TEMPERATURE); } catch { }
    if (isMD && Number.isFinite(temperature)) {
      const tempValue = Number(temperature);
      state.dynamics ||= {};
      state.dynamics.temperature = tempValue;
      setText('instTemp', `T: ${tempValue.toFixed(1)} K`);
    }

    ws.startSimulation({
      type: isMD ? 'md' : 'relax',
      params: isMD
        ? { calculator: opts.calculator || 'uma', temperature, timestep_fs: opts.timestep_fs ?? 1.0, friction: (Number.isFinite(opts.friction) ? opts.friction : cfg.mdFriction) }
        : { calculator: 'uma', fmax: state?.optimizer?.fmax || 0.05, max_step: MAX_STEP, optimizer: 'bfgs' },
    });
    if (env.apiOn()) dbg.log(`[${kind}WS][start]`, isMD ? { temperature, timestep_fs: opts.timestep_fs ?? 1.0 } : {});
    return { streaming: true };
  }

  const startMDContinuous = (opts) => startContinuous('md', opts);
  const startRelaxContinuous = (opts) => startContinuous('relax', opts);

  function stopSimulation() {
    try { const ws = getWS(); ws.stopSimulation(); } catch { }
    clearMdStreamListener();
    clearRelaxStreamListener();
    setMode(Mode.Idle); resetRPS();
  }

  function sendLiveMDParamsUpdate() {
    try {
      if (mode !== Mode.MD) return;
      const ws = getWS();
      let T = 1500;
      try {
        if (typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null) {
          const tRaw = Number(window.__MLIP_TARGET_TEMPERATURE);
          if (Number.isFinite(tRaw)) T = tRaw;
        }
      } catch { }
      const fr = (typeof window !== 'undefined' && Number.isFinite(window.__MLIP_CONFIG?.mdFriction))
        ? Number(window.__MLIP_CONFIG.mdFriction) : DEFAULT_MD_FRICTION;
      const dt = 1.0;
      if (Number.isFinite(T)) {
        const tempValue = Number(T);
        state.dynamics ||= {};
        state.dynamics.temperature = tempValue;
        setText('instTemp', `T: ${tempValue.toFixed(1)} K`);
      } else {
        setText('instTemp', TEMP_PLACEHOLDER);
      }
      ws.startSimulation({ type: 'md', params: { calculator: 'uma', temperature: T, timestep_fs: dt, friction: fr } });
      if (env.apiOn()) dbg.log('[mdWS][live-update]', { temperature: T, timestep_fs: dt, friction: fr });
    } catch { }
  }

  try {
    if (typeof window !== 'undefined') {
      window.addEventListener('mlip:temperature-changed', sendLiveMDParamsUpdate);
      window.addEventListener('mlip:friction-changed', sendLiveMDParamsUpdate);
    }
  } catch { }

  function setForceVectorsEnabled(on) {
    try {
      if (typeof on === 'boolean') { if (!!state.showForces !== on) state.toggleForceVectorsVisibility(); }
      else state.toggleForceVectorsVisibility();
    } catch { }
  }

  function getMetrics() { __count('index#getMetrics'); return { energy: state.dynamics?.energy, running: running.kind }; }
  function debugEnergySeriesLength() { return energyPlot.length(); }
  function debugRecordInteraction(kind) { recordInteraction(kind || 'debug'); }
  function getForceCacheVersion() { return state.forceCache?.version; }
  function getVersionInfo() { return { userInteractionVersion, totalInteractionVersion, resetEpoch }; }
  function debugGhostSnapshot() {
    const groups = view?._internals?.ghostBondGroups;
    if (!groups || typeof groups.entries !== 'function') return { ghostBondCount: 0, ghostGroups: [] };
    let ghostBondCount = 0;
    const ghostGroups = [];
    for (const [key, grp] of groups.entries()) {
      const count = grp && grp.mats ? grp.mats.length >>> 0 : 0;
      ghostBondCount += count;
      ghostGroups.push({ key, count });
    }
    return { ghostBondCount, ghostGroups };
  }
  function debugHighlightState() {
    const hi = view?._internals?.highlight;
    return {
      atomVisible: !!hi?.atom?.isVisible,
      bondVisible: !!hi?.bond?.isVisible,
    };
  }
  function debugCameraControls({ reset = false } = {}) {
    if (reset) {
      cameraControlStats.detachCalls = 0;
      cameraControlStats.attachCalls = 0;
    }
    return { detachCalls: cameraControlStats.detachCalls, attachCalls: cameraControlStats.attachCalls };
  }
  function debugBondMetrics() {
    const bonds = Array.isArray(state.bonds) ? state.bonds : [];
    const positions = Array.isArray(state.positions) ? state.positions : [];
    const elements = Array.isArray(state.elements) ? state.elements : [];
    const distance = (a, b) => {
      if (!a || !b) return NaN;
      const dx = (a.x || 0) - (b.x || 0);
      const dy = (a.y || 0) - (b.y || 0);
      const dz = (a.z || 0) - (b.z || 0);
      return Math.hypot(dx, dy, dz);
    };
    return bonds.map((bond, idx) => {
      const i = bond?.i ?? 0;
      const j = bond?.j ?? 0;
      return {
        index: idx,
        i,
        j,
        elements: [elements[i] ?? null, elements[j] ?? null],
        length: distance(positions[i], positions[j]),
        opacity: typeof bond?.opacity === 'number' ? bond.opacity : 1,
      };
    });
  }

  function handleSessionSnapshotApplied(payload = {}) {
    recomputePlaybackRuntime();
    const hasFrames = getFrameCount() > 0;
    if (hasFrames && Number.isFinite(playbackRuntime.startOffset) && timelineState.active) {
      applyTimelineFrame(playbackRuntime.startOffset);
      if (timelineUi) timelineUi.refresh();
    }
    const shouldAuto = timelinePlayback.shouldAutoPlay && timelinePlayback.shouldAutoPlay();
    if (shouldAuto && hasFrames) {
      timelinePlayback.start();
    } else {
      pauseTimelinePlayback();
    }
    if (EDIT_MODE) {
      timelineEditorPanel?.refresh?.();
      refreshEditorSelection();
      refreshEditorFrame(timelineState.offset);
    }
  }

  let libraryManifestCache = null;
  let libraryManifestPromise = null;

  async function fetchLibraryManifest() {
    if (libraryManifestCache) return libraryManifestCache;
    if (libraryManifestPromise) return libraryManifestPromise;
    libraryManifestPromise = (async () => {
      const url = '/examples/library.json';
      const res = await fetch(url, { cache: 'no-store' });
      if (!res.ok) throw new Error(`Library manifest fetch failed (${res.status})`);
      const data = await res.json();
      const normalized = Array.isArray(data)
        ? data.map((entry) => ({
            id: entry.id,
            label: entry.label || entry.id,
            description: entry.description || '',
            path: entry.path,
            tags: Array.isArray(entry.tags) ? entry.tags.slice() : [],
          }))
        : [];
      libraryManifestCache = normalized;
      libraryManifestPromise = null;
      return normalized;
    })().catch((err) => {
      libraryManifestPromise = null;
      throw err;
    });
    return libraryManifestPromise;
  }

  async function loadLibraryEntry(id) {
    if (!id) throw new Error('Missing library id');
    const manifest = await fetchLibraryManifest();
    const entry = manifest.find((item) => item.id === id);
    if (!entry) throw new Error(`Library entry not found: ${id}`);
    const res = await fetch(entry.path, { cache: 'no-store' });
    if (!res.ok) throw new Error(`Library session fetch failed (${res.status})`);
    const text = await res.text();
    if (!sessionStateManager?.loadFromFile) throw new Error('Session manager unavailable');
    await sessionStateManager.loadFromFile(text);
    return true;
  }

  if (!sessionStateManager) {
    sessionStateManager = createSessionStateManager({
      getViewerState: () => state,
      applyFullSnapshot,
      normalizeElement,
      energyPlot,
      frameBuffer,
      captureBaselineFromState: () => captureResetBaselineFromState(),
      installBaseline: installResetBaseline,
      getInteractionCounters: () => getInteractionCountersSnapshot(),
      setInteractionCounters: (counters) => setInteractionCountersFromSnapshot(counters),
      resetWsInitState,
      getWsClient: () => getWS(),
      ensureWsInit,
      posToTriples: (s) => posToTriples(s || state),
      zOf,
      stateCellToArray,
      getVelocitiesForSnapshot: () => getVelocitiesIfFresh(),
      setMode,
      Mode,
      timelineState,
      timelineUiRef,
      clearTimelinePlayback,
      setTimelineOverlayVisible,
      setTimelineInteractionLock,
      applyTimelineFrame,
      refreshTimelineUi: () => timelineUi?.refresh?.(),
      setTimelineUiMode: (modeName) => timelineUi?.setMode?.(modeName),
      lastContinuousOpts,
      rememberResume,
      stopSimulation,
      getMode: () => mode,
      controlMessageEngine: controlEngine,
      timelinePlayback,
      getTimelinePlaybackSnapshot: () => timelinePlayback.getSnapshot(),
      getTimelineControlSnapshot: () => controlEngine.getSnapshot?.() || [],
      applyTimelinePlaybackSnapshot: (cfg) => {
        timelinePlayback.applySnapshot(cfg || {});
        recomputePlaybackRuntime();
      },
      applyControlMessageSnapshot: (messages) => {
        controlEngine.applySnapshot(messages || []);
        controlEngine.refresh();
      },
      onSnapshotApplied: handleSessionSnapshotApplied,
      getMeshModes: () => ({
        atoms: view?.getAtomMeshModes?.(),
        bonds: view?.getBondMeshModes?.(),
      }),
      applyMeshModes: (modes) => view?.applyMeshModes?.(modes || {}),
      clearOpacityMask: (opts) => view?.clearOpacityMask?.(opts || { refresh: true }),
    });
    try { sessionStateManager.setBaselineFromState({ kind: 'init' }, { includeTimeline: false }); } catch { }
  }
  function setTestAutoSelectFallback(on = true) {
    if (view?._internals) {
      view._internals._debugAutoSelectFirstOnEmpty = !!on;
      return true;
    }
    return false;
  }
  function debugSelectAtom(index) {
    try {
      const idx = Number(index);
      selectionService.clickAtom(idx);
      return selectionService.get();
    } catch {
      return null;
    }
  }
  function debugSelectBond(ref) {
    try {
      if (!ref || typeof ref !== 'object') return null;
      const i = Number(ref.i);
      const j = Number(ref.j);
      if (!Number.isInteger(i) || !Number.isInteger(j)) return null;
      const payload = {
        i,
        j,
        index: ref.index != null ? ref.index : 0,
        key: ref.key != null ? ref.key : `${i}-${j}`,
      };
      selectionService.clickBond(payload);
      if (ref.orientation != null && state.selection?.kind === 'bond') {
        state.selection.data.orientation = ref.orientation;
      }
      if (ref.side && state.selection?.kind === 'bond') {
        state.selection.data.orientation = ref.side === 'i' ? 1 : 0;
      }
      return selectionService.get();
    } catch {
      return null;
    }
  }
  function debugGetSelection() {
    try {
      return selectionService.get();
    } catch {
      return null;
    }
  }

  function debugWsState() {
    return {
      reconnect: {
        attempts: reconnectState.attempts,
        nextAttemptAt: reconnectState.nextAttemptAt,
        bannerVisible: reconnectState.bannerVisible,
      },
      pendingResume: pendingResume.kind ? { kind: pendingResume.kind, opts: pendingResume.opts ? { ...pendingResume.opts } : {} } : null,
      log: wsStateLog.map((entry) => ({ ...entry })),
    };
  }

  function refreshResetBaseline(reason = 'manual') {
    seedResetBaseline(reason);
    try { sessionStateManager?.setBaselineFromState({ kind: reason }, { includeTimeline: false }); } catch { }
    return state.__resetBaseline;
  }

  function forceWsReconnect() {
    try {
      const ws = getWS();
      ws.reconnectNow?.();
      return true;
    } catch {
      return false;
    }
  }

  function simulateWsDrop({ failAttempts } = {}) {
    try {
      if (typeof failAttempts === 'number' && Number.isFinite(failAttempts)) {
        if (typeof window !== 'undefined') window.__MLIPVIEW_WS_FAIL_ATTEMPTS = Math.max(0, failAttempts | 0);
      }
    } catch { }
    try {
      const ws = getWS();
      ws.forceDisconnect?.('test-drop');
      return true;
    } catch {
      return false;
    }
  }

  async function resetToInitialPositions() {
    const t0 = env.now();
    dbg.log('[Reset] begin VR/AR-safe reset', { when: new Date().toISOString() });
    if (sessionStateManager) {
      try {
        const handled = await sessionStateManager.resetToLastSource();
        if (handled) {
          resetEpoch++; bumpUser('reset'); __suppressNextPosChange = true;
          interactions.length = 0;
          return true;
        }
      } catch (err) {
        if (env.apiOn()) dbg.warn('[Reset] session snapshot restore failed', err?.message || err);
      }
    }
    try {
      stopSimulation();
      resetEpoch++; bumpUser('reset'); __suppressNextPosChange = true;
      interactions.length = 0; energyPlot.reset();

      const baseline = state.__resetBaseline;
      let usedBaseline = false;
      if (baseline && Array.isArray(baseline.positions) && baseline.positions.length) {
        const baselinePositions = clonePositionList(baseline.positions);
        const baselineVelocities = baseline.velocities ? cloneVelocityList(baseline.velocities) : null;
        const baselineCellMatrix = baseline.cell ? cellSnapshotToMatrix(baseline.cell) : null;
        await applyFullSnapshot({
          elements: Array.isArray(baseline.elements) ? baseline.elements.map(normalizeElement) : [],
          positions: baselinePositions,
          velocities: baselineVelocities ?? undefined,
          cell: baselineCellMatrix || undefined,
        }, { updateBaseline: false });
        state.showCell = !!baseline.showCell;
        state.showGhostCells = !!baseline.showGhostCells;
        state.markCellChanged?.();
        state.__initialPositions = clonePositionList(baseline.positions);
        state.__initialCellSnapshot = cloneCellSnapshot(baseline.cell) || null;
        usedBaseline = true;
      }

      if (!usedBaseline) {
        const init = Array.isArray(state.__initialPositions) ? state.__initialPositions : null;
        if (!init || init.length !== state.positions.length) { dbg.warn('[Reset] missing or size-mismatch initial positions'); return false; }
        for (let i = 0; i < init.length; i++) { const p = init[i], tp = state.positions[i]; tp.x = p.x; tp.y = p.y; tp.z = p.z; }
        state.markPositionsChanged();
        try { if (state.dynamics) state.dynamics.velocities = []; } catch { }

        try {
          state.showCell = false;
          if (state.showGhostCells) state.showGhostCells = false;
          if (state.__initialCellSnapshot) {
            const c = state.__initialCellSnapshot;
            state.cell = {
              a: { x: c.a.x, y: c.a.y, z: c.a.z },
              b: { x: c.b.x, y: c.b.y, z: c.b.z },
              c: { x: c.c.x, y: c.c.y, z: c.c.z },
              originOffset: c.originOffset ? { x: c.originOffset.x || 0, y: c.originOffset.y || 0, z: c.originOffset.z || 0 } : { x: 0, y: 0, z: 0 },
              enabled: !!c.enabled,
            };
          } else {
            state.cell = null;
          }
          state.markCellChanged?.();
        } catch { }
      }

      try { await ff.computeForces({ sync: true }); } catch (e) { dbg.warn('[Reset] computeForces failed', e?.message || e); }
      energyPlot.reset();

      const t1 = env.now();
      dbg.log('[Reset] done', { totalMs: Math.round((t1 - t0) * 100) / 100 });
      return true;
    } catch (e) { dbg.warn('[Reset] exception', e?.message || e); return false; }
  }

  // Auto-start MD once (unless disabled or tests)
  try {
    const autoOk = (typeof window !== 'undefined') && !window.__MLIPVIEW_TEST_MODE && !window.__MLIPVIEW_NO_AUTO_MD;
    if (autoOk) {
      const MAX_AUTO_MD_ATTEMPTS = 120;
      const scheduleAutoMD = (attempt = 0, delayMs = 0) => {
        if (attempt > MAX_AUTO_MD_ATTEMPTS) return;
        setTimeout(() => {
          try {
            const api = window.viewerApi;
            if (!api || typeof api.startMDContinuous !== 'function') {
              scheduleAutoMD(attempt + 1, 75);
              return;
            }
            const metrics = api.getMetrics?.();
            if (metrics && metrics.running) return;

            try {
              const btn = document.getElementById('btnMDRun');
              if (btn) btn.textContent = 'stop';
            } catch { }

            Promise.resolve(api.startMDContinuous({}))
              .then(() => {
                try {
                  const btn = document.getElementById('btnMDRun');
                  if (btn && btn.textContent === 'stop') btn.textContent = 'run';
                } catch { }
                if (env.apiOn()) dbg.log('[autoMD] start invoked');
              })
              .catch((err) => {
                if (env.apiOn()) dbg.warn('[autoMD] start failed', err?.message || err);
                scheduleAutoMD(attempt + 1, 250);
              });
          } catch (err) {
            if (env.apiOn()) dbg.warn('[autoMD] start error', err?.message || err);
            scheduleAutoMD(attempt + 1, 150);
          }
        }, Math.max(0, delayMs | 0));
      };
      scheduleAutoMD();
    }
  } catch { }

  // Render loop
  let renderActive = true;
  const testMode = (typeof window !== 'undefined') && !!window.__MLIPVIEW_TEST_MODE;
  const rafLoop = () => {
    if (!renderActive) return;
    try { scene.render(); } catch (e) { dbg.warn('[Render][raf] render error', e); }
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(rafLoop);
  };
  try {
    if (!testMode && engine?.runRenderLoop) {
      dbg.log('[Render] starting engine.runRenderLoop (XR compatible)');
      engine.runRenderLoop(() => { if (!renderActive) return; try { scene.render(); } catch (e) { dbg.warn('[Render] loop error', e); } });
    } else {
      dbg.log('[Render] starting requestAnimationFrame loop (test mode fallback)');
      if (typeof requestAnimationFrame === 'function') requestAnimationFrame(rafLoop);
      else if (engine?.runRenderLoop) engine.runRenderLoop(() => { if (!renderActive) return; scene.render(); });
    }
  } catch (e) {
    dbg.warn('[Render] primary loop init failed; falling back to rAF', e);
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(rafLoop);
  }

  // Cleanup & dev helpers
  onWin(w => {
    w.__MLIPVIEW_CLEANUP ||= [];
    w.__MLIPVIEW_CLEANUP.push(() => {
      try {
        if (engine?.__mlipviewResizeHandler) {
          w.removeEventListener('resize', engine.__mlipviewResizeHandler);
          engine.__mlipviewResizeHandler = null;
        }
      } catch {}
      try { engine?.stopRenderLoop?.(); } catch { }
      try { scene?.dispose?.(); } catch { }
    });
    w.__MLIPVIEW_API_ENABLE = (on) => { w.__MLIPVIEW_DEBUG_API = !!on; dbg.log('[API] debug set to', w.__MLIPVIEW_DEBUG_API); };

    w.setDragThrottleMs = (ms) => {
      w.__MLIP_CONFIG.dragThrottleMs = Number(ms) || DEFAULTS.DRAG_THROTTLE_MS;
      dbg.log('[config] dragThrottleMs =', w.__MLIP_CONFIG.dragThrottleMs);
      return w.__MLIP_CONFIG.dragThrottleMs;
    };
    w.setPosEnergyDebounceMs = (ms) => {
      w.__MLIP_CONFIG.posEnergyDebounceMs = Number(ms) || DEFAULTS.POS_ENERGY_DEBOUNCE_MS;
      dbg.log('[config] posEnergyDebounceMs =', w.__MLIP_CONFIG.posEnergyDebounceMs);
      return w.__MLIP_CONFIG.posEnergyDebounceMs;
    };
    w.applyFullSnapshot = applyFullSnapshot;
    w.addAtom = addAtom;
    w.addAtomAtOrigin = addAtomAtOrigin;
    w.addAtomAtPosition = addAtomAtPosition;
    w.removeAtoms = removeAtoms;
    w.removeAtomByIndex = removeAtomByIndex;
    w.refreshResetBaseline = refreshResetBaseline;
  });

  function setForceProvider() { return 'uma'; }
  function shutdown() {
    renderActive = false;
    try { engine?.stopRenderLoop?.(); } catch { }
    try { calloutLayer?.dispose?.(); } catch { }
  }

  // Attach idle listener for idle mode
  try { attachIdleWSListener().catch(() => { }); } catch { }

  return {
    state,
    bondService,
    selection: selectionService,
    selectionService,
    ff,
    dynamics: { stepMD: () => { }, stepRelax: ({ forceFn }) => forceFn && forceFn() },
    view,
    picking,
    vr,
    applyFullSnapshot,
    addAtom,
    removeAtoms,
    recomputeBonds,
    relaxStep,
    mdStep,
    startRelaxContinuous,
    startMDContinuous,
    stopSimulation,
    setForceProvider,
    getMetrics,
    resetToInitialPositions,
    refreshResetBaseline,
    debugEnergySeriesLength,
    getEnergyMarker: () => energyPlot.getMarker(),
    debugRecordInteraction,
    manipulation: wrappedManipulation,
    scene,
    engine,
    camera,
    baselineEnergy,
    addAtomAtOrigin,
    addAtomAtPosition,
    removeAtomByIndex,
    setForceVectorsEnabled,
    getForceCacheVersion,
    getVersionInfo,
    shutdown,
    enableFeatureFlag,
    setMinStepInterval,
    requestSimpleCalculateNow,
      debugGhostSnapshot,
      debugHighlightState,
      debugStreamListenerStats,
      debugCameraControls,
    debugBondMetrics,
    setTestAutoSelectFallback,
    debugSelectAtom,
    debugSelectBond,
    debugGetSelection,
    debugWsState,
    forceWsReconnect,
    simulateWsDrop,
    session: {
      captureSnapshot: (meta) => sessionStateManager?.captureSnapshot?.(meta) ?? null,
      saveToFile: (opts) => sessionStateManager?.saveToFile?.(opts) ?? null,
      loadSnapshot: async (snapshot) => sessionStateManager ? sessionStateManager.loadSnapshot(snapshot) : null,
      loadFromFile: async (file) => sessionStateManager ? sessionStateManager.loadFromFile(file) : null,
      resetToLastLoad: () => sessionStateManager?.resetToLastSource?.() ?? false,
      noteSource: (meta) => sessionStateManager?.noteSource?.(meta) ?? null,
      getBaseline: () => sessionStateManager?.getBaseline?.() ?? null,
      getLibraryManifest: () => fetchLibraryManifest().then((list) => list.map((entry) => ({ ...entry }))),
      loadFromLibrary: (id) => loadLibraryEntry(id),
    },
    timeline: {
      select: (offset) => { handleTimelineOffsetRequest(offset); return timelineState.offset; },
      play: (offset) => { handleTimelinePlayRequest(offset ?? timelineState.offset); return timelineState.offset; },
      pause: () => { handleTimelinePauseRequest(); return timelineState.offset; },
      live: () => { handleTimelineLiveRequest(); return timelineState.offset; },
      getState: () => ({ ...(timelineUi?.getState?.() || {}), active: !!timelineState.active, playing: !!timelineState.playing, offset: timelineState.offset }),
      getSignature: (offset) => frameBuffer.getSignature(offset ?? timelineState.offset),
      bufferStats: () => frameBuffer.stats(),
      getOffsets: () => frameBuffer.listOffsets(),
      getPlaybackConfig: () => timelinePlayback.getBaseConfig(),
      getControlState: () => ({ ...activeControlState }),
      getFrameMeta: (offset) => {
        const off = Number.isFinite(offset) ? offset : timelineState.offset;
        if (!Number.isFinite(off)) return null;
        const entry = frameBuffer.getByOffset(off);
        if (!entry) return null;
        return {
          frameId: entry.id,
          offset: off,
          frameIndex: offsetToFrameIndex(off),
        };
      },
    },
    timelineEditor: EDIT_MODE ? {
      refresh: () => { timelineEditorPanel?.refresh?.(); return true; },
      getState: () => timelineEditorPanel?.getState?.() || null,
      select: (id) => timelineEditorPanel?.select?.(id) ?? false,
      getDraft: () => timelineEditorPanel?.getCurrentDraft?.() || null,
      status: () => ({
        frame: editorStatusFrame,
        selection: editorStatusSelection,
      }),
    } : null,
  };
}

/* ──────────────────────────────────────────────────────────────────────────
   Global viewerApi bootstrap
   ────────────────────────────────────────────────────────────────────────── */

onWin(w => {
  try {
    if (!w.viewerApi && typeof initNewViewer === 'function') {
      const __orig = initNewViewer;
      Object.defineProperty(w, 'initNewViewer', {
        value: async function (...args) {
          const api = await __orig(...args);
          w.viewerApi = api;
          return api;
        },
        configurable: true,
      });
    }
  } catch { }
  try {
    Object.defineProperty(w, '__dumpCurrentAtoms', {
      value: function () {
        try {
          if (!w.viewerApi) return null;
          const st = w.viewerApi.state;
          return {
            elements: st.elements.map(e => e.symbol || e.sym || e.S || e.Z || e.atomicNumber || '?'),
            atomic_numbers: st.elements.map(zOf),
            positions: st.positions.map(p => [p.x, p.y, p.z]),
            energy: st.dynamics?.energy,
          };
        } catch (e) { return { error: e?.message || String(e) }; }
      },
      writable: false,
    });
  } catch { }
});
