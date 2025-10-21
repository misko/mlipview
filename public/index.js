// index.js — streamlined viewer orchestration with WS/protobuf backend
// Behavior-compatible refactor to reduce verbosity & duplication.

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
import { SYMBOL_TO_Z } from './data/periodicTable.js';

/* ──────────────────────────────────────────────────────────────────────────
   Small utilities (logging, query flags, throttling, transforms, cache ops)
   ────────────────────────────────────────────────────────────────────────── */

const dbg = {
  apiOn: () => (typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API),
  log: (...a) => { try { console.log(...a); } catch { } },
  warn: (...a) => { try { console.warn(...a); } catch { } },
  err: (...a) => { try { console.error(...a); } catch { } },
};

const qbool = (key) => {
  try {
    const q = new URLSearchParams(window.location?.search || '');
    const v = q.get(key);
    return v === '1' || String(v).toLowerCase() === 'true';
  } catch { return false; }
};

const throttle = (fn, ms) => {
  let t = 0, id = null;
  return (...a) => {
    const now = (performance?.now?.() ?? Date.now());
    const since = now - t;
    if (id) return;
    const delay = Math.max(0, ms - since);
    id = setTimeout(() => { id = null; t = (performance?.now?.() ?? Date.now()); fn(...a); }, delay);
  };
};

const zOf = (e) =>
  (typeof e === 'number') ? e :
    (typeof e === 'string') ? (SYMBOL_TO_Z[e] || 0) :
      (e?.Z || e?.atomicNumber || e?.z ||
        SYMBOL_TO_Z[e?.symbol] || SYMBOL_TO_Z[e?.sym] || SYMBOL_TO_Z[e?.S] || 0);

const posToTriples = (state) => state.positions.map(p => [p.x, p.y, p.z]);

const applyTriples = (state, triples, { exclude } = {}) => {
  const N = Math.min(state.positions.length, triples?.length || 0);
  if (!N) return;
  if (exclude?.size) {
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

const stateCellToArray = (c) => (c && c.enabled) ? [
  [c.a.x, c.a.y, c.a.z],
  [c.b.x, c.b.y, c.b.z],
  [c.c.x, c.c.y, c.c.z],
] : null;

const setText = (id, txt) => { try { const el = document.getElementById(id); if (el) el.textContent = txt; } catch { } };

const onWin = (fn) => { try { if (typeof window !== 'undefined') fn(window); } catch { } };

/* ──────────────────────────────────────────────────────────────────────────
   Runtime config (min step interval, MD friction) — centralized & clamped
   ────────────────────────────────────────────────────────────────────────── */

if (typeof window !== 'undefined') {
  window.__MLIP_CONFIG ||= {};
  if (!Number.isFinite(window.__MLIP_CONFIG.minStepIntervalMs))
    window.__MLIP_CONFIG.minStepIntervalMs = DEFAULT_MIN_STEP_INTERVAL_MS;
  if (!Number.isFinite(window.__MLIP_CONFIG.mdFriction))
    window.__MLIP_CONFIG.mdFriction = DEFAULT_MD_FRICTION;

  // Optional toggles via query params
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
   Feature flags (mutable at runtime)
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
   Energy plot (tiny inline module)
   ────────────────────────────────────────────────────────────────────────── */

const energyPlot = (() => {
  let series = []; // {i,E,kind}
  let last = undefined;
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
    if (series.length === 0) return;
    if (series.length === 1) {
      const p = series[0]; ctx.fillStyle = '#6fc2ff'; ctx.beginPath();
      ctx.arc(W / 2, H / 2, 3, 0, Math.PI * 2); ctx.fill();
      if (label) label.textContent = 'E steps=1';
      return;
    }
    let minE = Infinity, maxE = -Infinity;
    for (const p of series) { if (p.E < minE) minE = p.E; if (p.E > maxE) maxE = p.E; }
    if (maxE - minE < 1e-12) maxE = minE + 1e-12;
    ctx.strokeStyle = '#6fc2ff'; ctx.lineWidth = 1; ctx.beginPath();
    for (let k = 0; k < series.length; k++) {
      const p = series[k];
      const x = (k / (series.length - 1)) * (W - 4) + 2;
      const y = H - 2 - ((p.E - minE) / (maxE - minE)) * (H - 4);
      if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    if (label) label.textContent = `E steps=${series.length} range=${(maxE - minE).toExponential(2)}`;
  };

  const push = (E, kind) => {
    if (typeof E !== 'number' || !isFinite(E)) return;
    if (E === last && !(typeof window !== 'undefined' && window.__ALLOW_DUPLICATE_ENERGY_TICKS)) return;
    const i = series.length; series.push({ i, E, kind }); last = E;
    if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_API) dbg.log('[energy] len=', series.length);
    setText('instEnergy', `E: ${E.toFixed(2)}`);
    draw();
  };

  const reset = () => { series = []; last = undefined; draw(); setText('instEnergy', 'E: —'); };
  const length = () => series.length;

  return { push, reset, length };
})();

/* ──────────────────────────────────────────────────────────────────────────
   Viewer initialization
   ────────────────────────────────────────────────────────────────────────── */

export async function initNewViewer(canvas, { elements, positions, bonds }) {
  __count('index#initNewViewer');

  // Debug toggles from query-string
  onWin(w => {
    if (w.__MLIPVIEW_DEBUG_API == null) w.__MLIPVIEW_DEBUG_API = qbool('debug');
    w.__MLIPVIEW_DEBUG_PICK = qbool('debugPick');
    w.__MLIPVIEW_DEBUG_SELECT = qbool('debugSelect');
    w.__MLIPVIEW_DEBUG_UI = qbool('debugUI');
  });

  const state = createMoleculeState({ elements, positions, bonds });
  // Cache initial positions for reset
  try { state.__initialPositions = (state.positions || []).map(p => ({ x: p.x, y: p.y, z: p.z })); } catch { }

  // Versioning & caches
  state.forceCache = { version: 0, energy: NaN, forces: [], stress: null, stale: true };
  let structureVersion = 0;
  let resetEpoch = 0;
  let userInteractionVersion = 0;
  let totalInteractionVersion = 0;

  const modifiedByVersion = new Map(); // v -> Set(indices)
  const draggingAtoms = new Set();
  const latchedUntil = new Map(); // idx -> ts
  const latchAtom = (i, ms = 600) => {
    const now = (performance?.now?.() ?? Date.now());
    latchedUntil.set(i, now + ms);
  };
  const currentDraggedAtom = { idx: null };

  const buildExcludeSet = () => {
    const out = new Set(draggingAtoms);
    const now = (performance?.now?.() ?? Date.now());
    for (const [i, t] of latchedUntil) { if (t > now) out.add(i); else latchedUntil.delete(i); }
    return out;
  };

  const bumpUser = (reason) => {
    userInteractionVersion++; totalInteractionVersion++; structureVersion++;
    if (state.forceCache) { state.forceCache.stale = true; state.forceCache.version = structureVersion; }
    if (dbg.apiOn()) dbg.log('[version][user]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
  };
  const bumpSim = (reason) => {
    totalInteractionVersion++;
    if (dbg.apiOn()) dbg.log('[version][sim]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
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
    if (typeof energy === 'number') {
      state.dynamics ||= {};
      state.dynamics.energy = energy;
    }
    const effectiveForces = Array.isArray(forces) ? forces : (state.forces || []);
    if (Array.isArray(forces)) updateForces(forces, { reason });
    state.forceCache = { version: structureVersion, energy: state.dynamics?.energy, forces: effectiveForces, stress, stale: false };
    if (dbg.apiOn()) dbg.log('[forces-cache]', reason, { n: effectiveForces.length, E: state.forceCache.energy });
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
  async function ensureWsInit() {
    const ws = getWS();
    await ws.ensureConnected();
    const nat = (state.elements || []).length;
    const ckey = cellKey();
    const need = !__wsState.inited || __wsState.lastAtomCount !== nat || __wsState.lastCellKey !== ckey;
    if (need) {
      const atomic_numbers = state.elements.map(zOf);
      const positionsTriples = posToTriples(state);
      const v = getVelocitiesIfFresh();
      const cell = stateCellToArray(state.cell);
      ws.userInteraction({ atomic_numbers, positions: positionsTriples, velocities: v || undefined, cell });
      __wsState.inited = true; __wsState.lastAtomCount = nat; __wsState.lastCellKey = ckey;
      if (dbg.apiOn()) dbg.log('[WS][ensureInit]', { nat, v: !!v, hasCell: !!cell });
    }
    return true;
  }

  // Simple force provider using WS idle compute
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
      await ensureWsInit();
      const positions = posToTriples(state);
      try { ws.setCounters({ userInteractionCount: userInteractionVersion }); } catch { }
      if (dbg.apiOn()) dbg.log('[WS][simple] USER_INTERACTION', { n: positions.length, uic: userInteractionVersion });
      ws.userInteraction({ positions });
      const t0 = (performance?.now?.() ?? Date.now());
      const { energy, forces, userInteractionCount: uicFromServer } = await ws.waitForEnergy({ timeoutMs: 5000 });
      const timingMs = Math.round(((performance?.now?.() ?? Date.now()) - t0) * 100) / 100;
      if (dbg.apiOn()) dbg.log('[WS][simple][recv]', { energy, forcesLen: forces?.length || 0, uicFromServer, currentUIC: userInteractionVersion, timingMs });
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
      const ws = getWS(); await ensureWsInit();
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

  // Recompute bonds when positions change
  const __origMarkPositionsChanged = state.markPositionsChanged?.bind(state) || null;
  const bondService = createBondService(state);
  const recomputeBonds = (reason = 'manual') => {
    __count('index#recomputeBonds');
    const bonds = bondService.recomputeAndStore();
    try { view.rebuildBonds(bonds); } catch { }
    if (reason !== 'markPositionsChanged') recordInteraction('rebonds');
    if (dbg.apiOn()) dbg.log(`[bonds] recomputed after ${reason} (count=${bonds ? bonds.length : 0})`);
    return bonds;
  };
  state.markPositionsChanged = (...a) => {
    structureVersion++;
    if (state.forceCache) { state.forceCache.stale = true; state.forceCache.version = structureVersion; }
    const r = __origMarkPositionsChanged ? __origMarkPositionsChanged(...a) : undefined;
    try { recomputeBonds('markPositionsChanged'); } catch { }
    return r;
  };
  if (!state.bonds?.length) { try { recomputeBonds(); } catch { } }

  // Scene & view
  const { engine, scene, camera } = await createScene(canvas);
  onWin(w => {
    w.__MLIPVIEW_DEBUG_TOUCH = !!w.__MLIPVIEW_DEBUG_TOUCH;
    w.enableTouchDebug = (on) => { w.__MLIPVIEW_DEBUG_TOUCH = !!on; dbg.log('[touchDebug] set to', w.__MLIPVIEW_DEBUG_TOUCH); };
  });
  const view = createMoleculeView(scene, state);
  const manipulation = createManipulationService(state, { bondService });

  // Picking & touch controls
  let pickingRef = null;
  const pickingProxy = { _debug: { get dragActive() { try { return !!(pickingRef?._debug?.dragActive); } catch { return false; } } } };
  const NO_TOUCH = qbool('noTouch');
  try {
    if (NO_TOUCH) { onWin(w => { w.__MLIPVIEW_NO_TOUCH = true; dbg.warn('[touchControls] skipped by ?noTouch=1'); }); }
    else { installTouchControls({ canvas, scene, camera, picking: pickingProxy }); dbg.log('[touchControls] installed'); }
  } catch (e) { dbg.warn('[touchControls] install failed', e?.message || e); }

  // Interaction bookkeeping
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
      const ws = getWS();
      ws.setCounters?.({ userInteractionCount: userInteractionVersion });
      ws.userInteraction?.({ positions: posToTriples(state) });
    } catch { }
  }, 100);

  // Wrapped manipulation to track edits
  function selectedAtomIndex() {
    try { return state.selection?.kind === 'atom' ? state.selection.data.index : null; }
    catch { return null; }
  }
  const wrappedManipulation = {
    beginDrag: (...a) => {
      const ok = manipulation.beginDrag?.(...a);
      if (ok) {
        const idx = selectedAtomIndex(); if (idx != null) { currentDraggedAtom.idx = idx; draggingAtoms.add(idx); }
      }
      return ok;
    },
    updateDrag: (...a) => {
      const r = manipulation.updateDrag?.(...a);
      if (r) {
        bumpUser('dragMove');
        const idx = selectedAtomIndex(); if (idx != null) { modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set()); modifiedByVersion.get(userInteractionVersion).add(idx); currentDraggedAtom.idx = idx; draggingAtoms.add(idx); }
        emitDuringDrag(); recordInteraction('dragMove');
      }
      return r;
    },
    endDrag: (...a) => {
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
          if (dragSource === 'vr') latchAtom(idx, 800);
        }
        try { const ws = getWS(); ws.setCounters?.({ userInteractionCount: userInteractionVersion }); ws.userInteraction?.({ positions: posToTriples(state) }); } catch { }
        currentDraggedAtom.idx = null;
      }
      recordInteraction('dragEnd');
      return r;
    },
    setDragPlane: (...a) => manipulation.setDragPlane?.(...a),
    rotateBond: (...a) => {
      const r = manipulation.rotateBond?.(...a);
      if (r) {
        let sideAtoms = null; try { sideAtoms = manipulation._debug?.getLastRotation?.()?.sideAtoms || null; } catch { }
        bumpUser('bondRotate');
        if (Array.isArray(sideAtoms) && sideAtoms.length) {
          modifiedByVersion.get(userInteractionVersion) || modifiedByVersion.set(userInteractionVersion, new Set());
          for (const i of sideAtoms) modifiedByVersion.get(userInteractionVersion).add(i);
        }
        if (!running.kind) ff.computeForces();
        recordInteraction('bondRotate');
        __suppressNextPosChange = true;
      }
      return r;
    },
  };
  try { Object.defineProperty(wrappedManipulation, '__isWrappedManipulation', { value: true }); } catch { wrappedManipulation.__isWrappedManipulation = true; }

  // Picking after wrappedManipulation available
  const picking = (pickingRef = createPickingService(
    scene,
    view,
    createSelectionService(state),
    {
      manipulation: new Proxy({}, { get: (_, k) => wrappedManipulation?.[k] }),
      camera,
      energyHook: ({ kind }) => {
        const dragging = draggingAtoms.size > 0 || currentDraggedAtom.idx != null;
        if (!running.kind) { if (dragging) emitDuringDrag(); else ff.computeForces(); }
        recordInteraction(kind || 'drag');
      },
    },
  ));

  // VR setup
  let vrPicker = null; try { vrPicker = createVRPicker({ scene, view }); } catch (e) { dbg.warn('[VR] vrPicker init failed', e?.message || e); }
  const vr = createVRSupport(scene, {
    picking: {
      ...picking,
      view,
      vrPicker,
      selectionService: picking.selectionService,
      manipulation: new Proxy({}, { get: (_, k) => wrappedManipulation?.[k] }),
      molState: state,
    },
  });
  try { vr.init().then(res => dbg.log(res.supported ? '[VR] support initialized (auto)' : '[VR] not supported')); } catch (e) { dbg.warn('[VR] auto init failed', e?.message || e); }

  // Baseline energy (seed forces/energy once)
  async function baselineEnergy() {
    __count('index#baselineEnergy');
    interactions.length = 0; energyPlot.reset();
    let res; try { res = await ff.computeForces({ sync: !!window?.__FORCE_SYNC_REMOTE__ }); state.dynamics ||= {}; state.dynamics.energy = res.energy; }
    catch { res = { energy: NaN, forces: [] }; }
    window.__RELAX_FORCES = res.forces || [];
    if (res.forces?.length) updateForces(res.forces, { reason: 'baselineEnergy' });
    // first energy point is added when API returns; keep baseline silent
  }
  try { await baselineEnergy(); } catch { }

  // Seed initial positions if missing
  try {
    const need = !state.__initialPositions?.length || state.__initialPositions.length !== (state.positions?.length || 0);
    if (need && state.positions?.length) {
      state.__initialPositions = state.positions.map(p => ({ x: p.x, y: p.y, z: p.z }));
      dbg.log('[Reset] baseline seeded after baselineEnergy', { count: state.__initialPositions.length });
    }
  } catch { }

  // Capture initial cell snapshot for reset
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

  // positionsChanged debounce (avoid overwriting step forces)
  let posEnergyTimer = null, pendingPosEnergy = false, skipFirstPosChangeForVersion = true;
  state.bus.on('positionsChanged', () => {
    const dragging = draggingAtoms.size > 0 || currentDraggedAtom.idx != null;
    if (dragging) { emitDuringDrag(); return; }
    pendingPosEnergy = true;
    if (posEnergyTimer) return;
    posEnergyTimer = setTimeout(() => {
      posEnergyTimer = null;
      if (!pendingPosEnergy) return; pendingPosEnergy = false;
      if (__suppressNextPosChange) { __suppressNextPosChange = false; return; }
      if (skipFirstPosChangeForVersion) { skipFirstPosChangeForVersion = false; if (!running.kind) ff.computeForces(); return; }
      if (!running.kind) ff.computeForces();
      bumpUser('posChange'); recordInteraction('posChange');
    }, 50);
  });

  // RPS counter (UI label optional)
  const rps = { samples: [], windowMs: 2000, value: 0 };
  const resetRPS = () => { rps.samples.length = 0; rps.value = 0; setText('rpsLabel', 'RPS: --'); };
  const noteReqDone = () => {
    const now = (performance?.now?.() ?? Date.now());
    rps.samples.push(now);
    while (rps.samples.length && now - rps.samples[0] > rps.windowMs) rps.samples.shift();
    rps.value = (rps.samples.length >= 2) ? ((rps.samples.length - 1) * 1000) / (rps.samples[rps.samples.length - 1] - rps.samples[0]) : 0;
    setText('rpsLabel', 'RPS: ' + rps.value.toFixed(1));
  };

  // Helpers for app focus (desktop)
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
    const ws = getWS(); await ensureWsInit();
    try { ws.setCounters({ userInteractionCount: userInteractionVersion, simStep: 0 }); } catch { }
    try {
      const v = getVelocitiesIfFresh();
      ws.userInteraction({ positions: posToTriples(state), velocities: v || undefined });
    } catch { }
    const r = await ws.requestSingleStep({ type: kind, params });
    const epochAtSend = resetEpoch;
    if (epochAtSend !== resetEpoch) return { stale: true, staleReason: 'staleEpoch', epochAtSend, resetEpoch };
    if (Array.isArray(r.positions) && r.positions.length === state.positions.length) {
      applyTriples(state, r.positions);
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
    try {
      const r = await doOneStepViaWS('relax', { calculator: 'uma', fmax: state?.optimizer?.fmax || 0.05, max_step: MAX_STEP, optimizer: 'bfgs' });
      if (r?.applied) { recordInteraction('relaxStep'); energyPlot.push(state.dynamics?.energy, 'relax'); }
      return r;
    } catch (e) { dbg.warn('[relaxStep] failed', e); return { error: e?.message || String(e) }; }
  }

  async function mdStep(opts = {}) {
    __count('index#mdStep');
    try {
      // live T override if present
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
     Streaming WS loops
     ──────────────────────────────────────────────────────────────────────── */

  function handleStreamFrame(kind, ws, r) {
    // discard stale UIC only if server echoes a positive lagging count
    if (Number.isFinite(r.userInteractionCount) && r.userInteractionCount > 0 && r.userInteractionCount < userInteractionVersion) {
      if (dbg.apiOn()) dbg.log(`[${kind}WS][discard-stale]`, { server: r.userInteractionCount, current: userInteractionVersion });
      return;
    }
    try { const seq = Number(r.seq) || 0; if (seq > 0) ws.ack(seq); } catch { }

    const exclude = buildExcludeSet();
    if (Array.isArray(r.positions) && r.positions.length === state.positions.length) {
      applyTriples(state, r.positions, { exclude });
      __suppressNextPosChange = true; state.markPositionsChanged(); bumpSim(kind === 'md' ? 'mdWS' : 'relaxWS');
    }
    const E = (typeof r.energy === 'number') ? r.energy : state.dynamics?.energy;
    updateEnergyForces({ energy: E, forces: r.forces?.length ? r.forces : undefined, reason: kind + 'WS' });
    if (kind === 'md') {
      if (Array.isArray(r.velocities) && r.velocities.length === state.elements.length) {
        state.dynamics ||= {}; state.dynamics.velocities = r.velocities;
      }
      if (typeof r.temperature === 'number' && isFinite(r.temperature)) {
        state.dynamics ||= {}; state.dynamics.temperature = r.temperature;
        setText('instTemp', 'T: ' + r.temperature.toFixed(1) + ' K');
      }
      // test-mode fallback: ensure velocities shape exists
      try {
        const tm = typeof window !== 'undefined' && !!window.__MLIPVIEW_TEST_MODE;
        const haveV = Array.isArray(state.dynamics?.velocities) && state.dynamics.velocities.length === state.elements.length;
        if (tm && !haveV && state.positions?.length === state.elements.length) {
          state.dynamics ||= {}; state.dynamics.velocities = state.positions.map(() => [0, 0, 0]);
        }
      } catch { }
    }
    energyPlot.push(E, kind);
    noteReqDone();
  }

  let running = { kind: null, abort: null };

  async function startRelaxContinuous({ maxSteps = 1000 } = {}) {
    if (!FEATURES.RELAX_LOOP) { dbg.warn('[feature] RELAX_LOOP disabled'); return { disabled: true }; }
    if (running.kind) return { ignored: true };
    running.kind = 'relax';

    const ws = getWS(); await ensureWsInit();
    try { if (posEnergyTimer) { clearTimeout(posEnergyTimer); posEnergyTimer = null; } pendingPosEnergy = false; } catch { }
    try { ws.setCounters({ userInteractionCount: userInteractionVersion, simStep: 0 }); } catch { }
    try {
      const v = getVelocitiesIfFresh();
      ws.userInteraction({ positions: posToTriples(state), velocities: v || undefined });
    } catch { }
    resetRPS();

    let taken = 0, unsub = null, stopped = false;
    unsub = ws.onResult((r) => {
      if (stopped || running.kind !== 'relax') return;
      if (dbg.apiOn()) dbg.log('[relaxWS][frame]', { seq: Number(r.seq) || 0, uic: Number(r.userInteractionCount) || 0, simStep: Number(r.simStep) || 0, have: { pos: !!r.positions, forces: !!r.forces, energy: typeof r.energy === 'number' } });
      handleStreamFrame('relax', ws, r);
      if (++taken >= maxSteps) { stopped = true; try { ws.stopSimulation(); } catch { } try { unsub && unsub(); } catch { } running.kind = null; resetRPS(); }
    });

    ws.startSimulation({ type: 'relax', params: { calculator: 'uma', fmax: state?.optimizer?.fmax || 0.05, max_step: MAX_STEP, optimizer: 'bfgs' } });
    if (dbg.apiOn()) dbg.log('[relaxWS][start]');
    return { streaming: true };
  }

  async function startMDContinuous({ steps = 1000, calculator = 'uma', temperature = 1500, timestep_fs = 1.0, friction } = {}) {
    if (!FEATURES.MD_LOOP) { dbg.warn('[feature] MD_LOOP disabled'); return { disabled: true }; }
    if (running.kind) return { ignored: true };
    running.kind = 'md';

    const ws = getWS(); await ensureWsInit();
    try { if (posEnergyTimer) { clearTimeout(posEnergyTimer); posEnergyTimer = null; } pendingPosEnergy = false; } catch { }
    try { ws.setCounters({ userInteractionCount: userInteractionVersion, simStep: 0 }); } catch { }
    try {
      const v = getVelocitiesIfFresh();
      ws.userInteraction({ positions: posToTriples(state), velocities: v || undefined });
    } catch { }
    // test-mode: seed zero velocities
    try {
      const tm = typeof window !== 'undefined' && !!window.__MLIPVIEW_TEST_MODE;
      if (tm && state.positions?.length === state.elements.length) {
        state.dynamics ||= {}; state.dynamics.velocities = state.positions.map(() => [0, 0, 0]);
      }
    } catch { }
    resetRPS();

    let taken = 0, unsub = null, stopped = false;
    unsub = ws.onResult((r) => {
      if (stopped || running.kind !== 'md') return;
      if (dbg.apiOn()) dbg.log('[mdWS][frame]', { seq: Number(r.seq) || 0, uic: Number(r.userInteractionCount) || 0, simStep: Number(r.simStep) || 0, have: { pos: !!r.positions, forces: !!r.forces, energy: typeof r.energy === 'number', vel: !!r.velocities } });
      handleStreamFrame('md', ws, r);
      const extraAllowance = 5;
      const haveV = Array.isArray(state.dynamics?.velocities) && state.dynamics.velocities.length === state.elements.length;
      if (++taken >= steps && (haveV || taken - steps >= extraAllowance)) {
        stopped = true; try { ws.stopSimulation(); } catch { } try { unsub && unsub(); } catch { }
        running.kind = null; resetRPS();
      }
    });

    // live temperature from global target if present
    try { if (typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null) temperature = Number(window.__MLIP_TARGET_TEMPERATURE); } catch { }
    ws.startSimulation({
      type: 'md',
      params: {
        calculator,
        temperature,
        timestep_fs,
        friction: (Number.isFinite(friction) ? friction : cfg.mdFriction),
      },
    });
    if (dbg.apiOn()) dbg.log('[mdWS][start]', { temperature, timestep_fs, friction: (Number.isFinite(friction) ? friction : cfg.mdFriction) });
    return { streaming: true };
  }

  function stopSimulation() {
    try { const ws = getWS(); ws.stopSimulation(); } catch { }
    running.kind = null; resetRPS();
  }

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

  async function resetToInitialPositions() {
    const t0 = (performance?.now?.() ?? Date.now());
    dbg.log('[Reset] begin VR/AR-safe reset', { when: new Date().toISOString() });
    try {
      stopSimulation();
      resetEpoch++; bumpUser('reset'); __suppressNextPosChange = true;
      interactions.length = 0; energyPlot.reset();

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

      const f0 = (performance?.now?.() ?? Date.now());
      try { await ff.computeForces({ sync: true }); } catch (e) { dbg.warn('[Reset] computeForces failed', e?.message || e); }
      const f1 = (performance?.now?.() ?? Date.now());
      energyPlot.reset();

      const t1 = (performance?.now?.() ?? Date.now());
      dbg.log('[Reset] done', { totalMs: Math.round((t1 - t0) * 100) / 100, recomputeMs: Math.round((f1 - f0) * 100) / 100 });
      return true;
    } catch (e) { dbg.warn('[Reset] exception', e?.message || e); return false; }
  }

  // Auto-start MD once (unless disabled or tests)
  try {
    const autoOk = (typeof window !== 'undefined') && !window.__MLIPVIEW_TEST_MODE && !window.__MLIPVIEW_NO_AUTO_MD;
    if (autoOk) {
      setTimeout(() => {
        try {
          if (!window.viewerApi) return;
          const m = window.viewerApi.getMetrics(); if (m.running) return;
          window.viewerApi.startMDContinuous({}).then(() => {
            try { const btn = document.getElementById('btnMDRun'); if (btn && btn.textContent === 'stop') btn.textContent = 'run'; } catch { }
            if (dbg.apiOn()) dbg.log('[autoMD] completed initial MD run');
          });
          try { const btn = document.getElementById('btnMDRun'); if (btn) btn.textContent = 'stop'; } catch { }
        } catch (e) { dbg.warn('[autoMD] start failed', e?.message || e); }
      }, 0);
    }
  } catch { }

  // Render loop (XR-friendly)
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

  // Cleanup for Jest/teardown
  onWin(w => {
    w.__MLIPVIEW_CLEANUP ||= [];
    w.__MLIPVIEW_CLEANUP.push(() => {
      try { engine?.stopRenderLoop?.(); } catch { }
      try { scene?.dispose?.(); } catch { }
    });
    w.__MLIPVIEW_API_ENABLE = (on) => { w.__MLIPVIEW_DEBUG_API = !!on; dbg.log('[API] debug set to', w.__MLIPVIEW_DEBUG_API); };
  });

  function setForceProvider() { return 'uma'; }
  function shutdown() { renderActive = false; try { engine?.stopRenderLoop?.(); } catch { } }

  return {
    state,
    bondService,
    selection: picking.selectionService,
    ff,
    dynamics: { stepMD: () => { }, stepRelax: ({ forceFn }) => forceFn && forceFn() },
    view,
    picking,
    vr,
    recomputeBonds,
    relaxStep,
    mdStep,
    startRelaxContinuous,
    startMDContinuous,
    stopSimulation,
    setForceProvider,
    getMetrics,
    resetToInitialPositions,
    debugEnergySeriesLength,
    debugRecordInteraction,
    manipulation: wrappedManipulation,
    scene,
    engine,
    camera,
    baselineEnergy,
    setForceVectorsEnabled,
    getForceCacheVersion,
    getVersionInfo,
    shutdown,
    enableFeatureFlag,
    setMinStepInterval,
    requestSimpleCalculateNow,
  };
}

/* ──────────────────────────────────────────────────────────────────────────
   Global viewerApi bootstrap for tests / devtools
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
  // small debug helper
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
