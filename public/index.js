import { createMoleculeState } from './domain/moleculeState.js';
import { getEndpointSync } from './api_endpoints.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { createFairChemForcefield } from './fairchem_provider.js';
import { createScene } from './render/scene.js';
import { createMoleculeView } from './render/moleculeView.js';
import { createPickingService } from './core/pickingService.js';
import installTouchControls from './ui/touchControls.js';
import { createManipulationService } from './domain/manipulationService.js';
import { createVRSupport } from './vr/setup.js';
import { createVRPicker } from './vr/vr-picker.js';
import { __count } from './util/funcCount.js';

// --- Runtime Config (pacing, etc.) ---
if (typeof window !== 'undefined') {
  // Add mdFriction default if not present
  window.__MLIP_CONFIG = window.__MLIP_CONFIG || { minStepIntervalMs: 30, mdFriction: 0.5 };
  if (window.__MLIP_CONFIG.mdFriction == null) window.__MLIP_CONFIG.mdFriction = 0.5;
}
function getConfig(){
  // Always return a valid config object with defaults; repair globals if missing
  if (typeof window === 'undefined') return { minStepIntervalMs:30, mdFriction:0.5 };
  const cfg = window.__MLIP_CONFIG || {};
  if (cfg.minStepIntervalMs == null || !Number.isFinite(cfg.minStepIntervalMs)) cfg.minStepIntervalMs = 30;
  if (cfg.mdFriction == null || !Number.isFinite(cfg.mdFriction)) cfg.mdFriction = 0.5;
  window.__MLIP_CONFIG = cfg;
  return cfg;
}
export function setMinStepInterval(ms){
  const v = Number(ms);
  if(!Number.isFinite(v) || v < 0) throw new Error('minStepIntervalMs must be non-negative number');
  if (typeof window !== 'undefined') {
    window.__MLIP_CONFIG.minStepIntervalMs = Math.max(1, Math.round(v));
    console.log('[config] set minStepIntervalMs =', window.__MLIP_CONFIG.minStepIntervalMs);
    return window.__MLIP_CONFIG.minStepIntervalMs;
  }
  return Math.max(1, Math.round(v));
}
export function setMdFriction(f){
  const v = Number(f);
  if(!Number.isFinite(v) || v < 0) throw new Error('mdFriction must be a non-negative number');
  const clamped = Math.max(0, Math.min(5, v)); // safety clamp: [0,5]
  if (typeof window !== 'undefined') {
    window.__MLIP_CONFIG.mdFriction = clamped;
    console.log('[config] set mdFriction =', window.__MLIP_CONFIG.mdFriction);
    return window.__MLIP_CONFIG.mdFriction;
  }
  return clamped;
}

// Debug no-cache mode deprecated; server-side cache removed.
function __noCacheMode(){ return true; }

// Minimal periodic table map for elements we currently load; extend as needed.
const SYMBOL_TO_Z = {
  H:1, He:2, C:6, N:7, O:8, F:9, Ne:10, S:16, Cl:17, Fe:26, Co:27, Ni:28, Cu:29, Zn:30,
  Si:14, P:15, Br:35, I:53, Au:79, Ag:47, Al:13, Ca:20, K:19, Na:11, Li:3
};
function elementToZ(e){
  if (e==null) return 0;
  if (typeof e === 'number') return e; // already Z
  if (typeof e === 'string') return SYMBOL_TO_Z[e] || 0;
  return e.Z || e.atomicNumber || e.z || (SYMBOL_TO_Z[e.symbol] || SYMBOL_TO_Z[e.sym] || SYMBOL_TO_Z[e.S] || 0) || 0;
}

// Resolve API base URL across browser + Jest (node) environments.
function __resolveApiBase(){
  // Browser environment: prefer explicit override, then same-origin, then test globals.
  if (typeof window !== 'undefined') {
    // 1. Explicit override (developer can set before viewer init): window.__MLIPVIEW_SERVER = 'https://host:port'
    if (window.__MLIPVIEW_SERVER) return String(window.__MLIPVIEW_SERVER).replace(/\/$/, '');
    // 2. If running in a real browser (or served via proxy) use same-origin so calls stay on the page host (avoids mixed-content & CORS).
    try {
      if (window.location && window.location.origin && window.location.origin !== 'null') {
        // In jsdom this is usually 'http://localhost'. We still may want a better test base; handle below.
        const origin = window.location.origin;
        // If jsdom default AND a test global is provided, defer to test global.
        if (origin === 'http://localhost' && typeof global !== 'undefined' && global.__MLIP_API_URL) {
          return global.__MLIP_API_URL.replace(/\/$/, '');
        }
        return origin.replace(/\/$/, '');
      }
    } catch { /* ignore */ }
    // 3. Fallback: if test harness placed API URL on global (jsdom env) use it.
    if (typeof global !== 'undefined' && global.__MLIP_API_URL) {
      return global.__MLIP_API_URL.replace(/\/$/, '');
    }
  }
  // Node/Jest (no window) fallback: prefer globals or env vars injected by harness.
  if (typeof global !== 'undefined') {
    const g = global.__MLIP_API_URL || process.env.MLIP_API_URL || process.env.MLIPVIEW_SERVER;
    if (g) return g.replace(/\/$/, '');
  }
  // Final hard-coded development fallback (direct Ray Serve default port)
  return 'http://127.0.0.1:8000';
}

export async function initNewViewer(canvas, { elements, positions, bonds } ) {
  __count('index#initNewViewer');
  const state = createMoleculeState({ elements, positions, bonds });
  // We'll wrap markPositionsChanged after bondService + view exist so we can recompute bonds centrally.
  const __origMarkPositionsChanged = state.markPositionsChanged ? state.markPositionsChanged.bind(state) : null;
  const bondService = createBondService(state);
  const selection = createSelectionService(state);
  // Remote UMA force provider via /simple_calculate
  let lastForceResult = { energy: NaN, forces: [] };
  let inFlight = false;
  // --- API debug instrumentation ---
  if (window.__MLIPVIEW_DEBUG_API == null) {
    // Default OFF; enable with ?debug=1 or ?debug=true
    try {
      const q = new URLSearchParams(window.location.search);
      const dbg = q.get('debug');
      window.__MLIPVIEW_DEBUG_API = (dbg === '1' || dbg === 'true');
    } catch { window.__MLIPVIEW_DEBUG_API = false; }
  }
  // Picking debug via ?debugPick=1
  try {
    const q = new URLSearchParams(window.location.search||'');
    const dp = q.get('debugPick');
    window.__MLIPVIEW_DEBUG_PICK = (dp === '1' || dp === 'true');
    if (window.__MLIPVIEW_DEBUG_PICK) console.log('[debug] pick logging enabled');
    const ds = q.get('debugSelect');
    window.__MLIPVIEW_DEBUG_SELECT = (ds === '1' || ds === 'true');
    if (window.__MLIPVIEW_DEBUG_SELECT) console.log('[debug] selection logging enabled');
    const du = q.get('debugUI');
    window.__MLIPVIEW_DEBUG_UI = (du === '1' || du === 'true');
    if (window.__MLIPVIEW_DEBUG_UI) console.log('[debug] UI logging enabled');
  } catch {}
  // Latency debugging gate comes from index.html ?debugLatency=1 wrapper
  const __latencyOn = (typeof window!=='undefined') && !!window.__MLIP_DEBUG_LATENCY;
  function noteLatency(kind, phase, meta){
    if(!__latencyOn) return;
    try {
      const s = new Date().toISOString();
      console.log(`[latency:${kind}] ${phase} ${s}`, meta||{});
    } catch {}
  }
  let __apiSeq = window.__MLIPVIEW_API_SEQ || 0;
  function debugApi(kind, phase, meta){
    __count('index#debugApi');
    if(!window.__MLIPVIEW_DEBUG_API) return;
    try {
      const stamp = new Date().toISOString();
      // We copy a shallow subset to avoid huge spam (positions trimmed)
      if(meta && meta.body && meta.body.coordinates && meta.body.coordinates.length>30){
        meta = { ...meta, body: { ...meta.body, coordinates:`<${meta.body.coordinates.length} coords>` } };
      }
      if(meta && meta.response && meta.response.positions && meta.response.positions.length>30){
        meta = { ...meta, response: { ...meta.response, positions:`<${meta.response.positions.length} positions>` } };
      }
      // Energy summary convenience
      if(meta && meta.response && typeof meta.response.final_energy === 'number'){
        meta.energy = meta.response.final_energy;
      } else if(meta && meta.response && meta.response.results && typeof meta.response.results.energy === 'number'){
        meta.energy = meta.response.results.energy;
      }
      console.log(`[API][${kind}][${phase}]#${meta.seq} t=${meta.timingMs!=null?meta.timingMs+'ms':''} ${stamp}`, meta);
    } catch(e){ /* ignore logging errors */ }
  }
  window.__MLIPVIEW_API_ENABLE = function(on){ window.__MLIPVIEW_DEBUG_API = !!on; console.log('[API] debug set to', window.__MLIPVIEW_DEBUG_API); };
  window.__MLIPVIEW_API_SEQ = __apiSeq;

  // Versioned force cache: state.forceCache = { version, energy, forces, stress, stale }
  state.forceCache = { version:0, energy:NaN, forces:[], stress:null, stale:true };
  let structureVersion = 0; // increments when positions (or elements) change
  // Server-side atoms cache removed: no cache key tracking
  // Version counters:
  // userInteractionVersion: increments ONLY on user geometry edits (drag, bond rotate, debounced posChange)
  // totalInteractionVersion: increments on user edits AND accepted simulation (relax/md) steps.
  let userInteractionVersion = 0;
  let totalInteractionVersion = 0;
  // Track which atoms were modified by each user interaction version increment.
  // Map: version -> Set(atomIndices)
  const __modifiedAtomsByVersion = new Map();
  // Track atoms currently being dragged (by index). Only atom-drags are considered here.
  const __draggingAtoms = new Set();
  // Recently interacted atoms latch (atomIndex -> expiresAt ms). Prevents in-flight MD/relax from overriding
  // the just-dragged atom if a response lands immediately after stick release.
  const __latchedAtomsUntil = new Map();
  function __latchAtom(i, ms=600){
    try {
      const now = (typeof performance!=='undefined' && performance.now) ? performance.now() : Date.now();
      __latchedAtomsUntil.set(i, now + ms);
    } catch {}
  }
  let __currentDraggedAtomIndex = null;
  function __addModifiedAtoms(indices){
    if(!indices || !indices.length) return;
    const v = userInteractionVersion;
    let s = __modifiedAtomsByVersion.get(v);
    if(!s){ s = new Set(); __modifiedAtomsByVersion.set(v, s); }
    for(const i of indices){ if(Number.isInteger(i) && i>=0 && i<state.positions.length) s.add(i); }
  }
  function __modifiedSetBetween(vFromExclusive, vToInclusive){
    const out = new Set();
    for(let v=vFromExclusive+1; v<=vToInclusive; v++){
      const s = __modifiedAtomsByVersion.get(v);
      if(s) for(const i of s) out.add(i);
    }
    return out;
  }
  function bumpUserInteractionVersion(reason){
    userInteractionVersion++; totalInteractionVersion++; structureVersion++;
    if(state.forceCache) state.forceCache.stale = true;
    if(window.__MLIPVIEW_DEBUG_API) console.debug('[version][user]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
  }
  function bumpSimulationVersion(reason){
    // Simulation application affects totalInteractionVersion but structureVersion already bumped via markPositionsChanged.
    totalInteractionVersion++;
    if(window.__MLIPVIEW_DEBUG_API) console.debug('[version][sim]', { reason, userInteractionVersion, totalInteractionVersion, structureVersion });
  }
  function __updateForces(forces, { reason }={}) {
    const DBG = (typeof window !== 'undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search||''));
    try {
      if (Array.isArray(forces) && forces.length) {
        state.forces = forces;
        if (DBG) console.log('[Forces][update]', reason||'?', 'len=', forces.length);
        state.bus?.emit && state.bus.emit('forcesChanged');
      } else if (DBG) {
        console.warn('[Forces][update] empty set for reason', reason);
      }
    } catch(e){ if (DBG) console.warn('[Forces][update] error', e); }
  }
  async function fetchRemoteForces({ awaitResult=false }={}){
    __count('index#fetchRemoteForces');
    if(inFlight && !awaitResult) return;
    while(inFlight && awaitResult) { await new Promise(r=>setTimeout(r,10)); }
    inFlight = true;
    try {
  // Server-side cache removed; do not short-circuit on local cache here.
  // Skip remote call until we actually have atoms; prevents 500 errors on empty initial state
  if(!state.elements || state.elements.length === 0){
        return lastForceResult;
      }
  const atomic_numbers = state.elements.map(e=> elementToZ(e));
      const coordinates = state.positions.map(p=>[p.x,p.y,p.z]);
      const body = { atomic_numbers, coordinates, properties:['energy','forces'], calculator:'uma' };
    const base = __resolveApiBase();
    const url = base + getEndpointSync('simple');
      const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
      const build0 = performance.now();
      // payload already prepared above
      const t0 = performance.now();
      debugApi('simple_calculate','request',{ seq, url, body });
      let resp, json;
      try {
  resp = await fetch(url, { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body) });
      } catch(netErr){
        debugApi('simple_calculate','network-error',{ seq, url, error: netErr?.message||String(netErr) });
        throw netErr;
      }
      const t1 = performance.now();
      const timingMs = Math.round((t1 - t0)*100)/100;
      noteLatency('simple', 'network', { seq, ms: timingMs });
      if(resp.ok){
        try {
          const p0=performance.now();
          json = await resp.json();
          const p1=performance.now(); noteLatency('simple','parse',{ seq, ms: Math.round((p1-p0)*100)/100 });
        } catch(parseErr){ debugApi('simple_calculate','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); }
        debugApi('simple_calculate','response',{ seq, url, timingMs, status: resp.status, response: json });
      } else {
        const txt = await resp.text();
        debugApi('simple_calculate','http-error',{ seq, url, timingMs, status: resp.status, statusText: resp.statusText, bodyText: txt });
      }
      if(resp.ok){
        const res = json?.results||{};
        if(typeof res.energy === 'number'){
          lastForceResult = { energy: res.energy, forces: res.forces||[] };
          state.dynamics = state.dynamics || {}; state.dynamics.energy = res.energy;
          window.__RELAX_FORCES = res.forces||[];
          state.forceCache = { version: structureVersion, energy: res.energy, forces: res.forces||[], stress: res.stress||null, stale:false };
          // Server cache key removed
          // Commit forces to state & notify renderer each successful fetch so visualization updates.
          if (Array.isArray(res.forces) && res.forces.length) {
            __updateForces(res.forces, { reason:'fetchRemoteForces' });
          }
          maybePlotEnergy('forces');
        }
      }
    } catch(_e) { /* silent for now */ } finally { inFlight=false; }
    return lastForceResult;
  }
  const ff = { computeForces: ({ sync }={})=>{ if(sync) return fetchRemoteForces({ awaitResult:true }); fetchRemoteForces(); return lastForceResult; } };

  // Helper: construct precomputed object if force cache is fresh and matches structure
  function buildPrecomputedIfFresh(){
    const fc = state.forceCache;
    if(!fc || fc.stale) return null;
    if(fc.version !== structureVersion) return null;
    const n = state.elements.length;
    if(!Array.isArray(fc.forces) || fc.forces.length !== n) return null;
    if(!Number.isFinite(fc.energy)) return null;
    const pre = { energy: fc.energy, forces: fc.forces };
    if (fc.stress && Array.isArray(fc.stress) && fc.stress.length === 6) pre.stress = fc.stress;
    return pre;
  }

  function getVelocitiesIfFresh(){
    const dyn = state.dynamics || {};
    const v = dyn.velocities;
    if(!v) return null;
    if(!Array.isArray(v) || v.length !== state.elements.length) return null;
    // Basic finiteness check without allocating large arrays
    for(let i=0;i<v.length;i++){ const row=v[i]; if(!row || row.length!==3) return null; for(let k=0;k<3;k++){ if(!Number.isFinite(row[k])) return null; } }
    return v;
  }

  const dynamics = { stepMD: ()=>{}, stepRelax: ({ forceFn })=>{ forceFn(); } };
  // --- Requests-per-second (RPS) counter for MD/Relax runs ---
  const __rps = { samples: [], windowMs: 2000, value: 0 };
  function __resetRPS(){
    __rps.samples.length = 0; __rps.value = 0;
    try { const el = document.getElementById('rpsLabel'); if (el) el.textContent = 'RPS: --'; } catch {}
  }
  function __noteRequestCompleted(){
    const now = (typeof performance!=='undefined' && performance.now) ? performance.now() : Date.now();
    __rps.samples.push(now);
    // Drop samples outside the window
    while(__rps.samples.length && (now - __rps.samples[0]) > __rps.windowMs){ __rps.samples.shift(); }
    if(__rps.samples.length >= 2){
      const dt = (__rps.samples[__rps.samples.length-1] - __rps.samples[0]);
      __rps.value = dt > 0 ? (__rps.samples.length - 1) * 1000 / dt : 0;
    } else {
      __rps.value = 0;
    }
    try { const el = document.getElementById('rpsLabel'); if (el) el.textContent = 'RPS: ' + __rps.value.toFixed(1); } catch {}
  }
  // Remote relaxation: call backend /relax
  async function callRelaxEndpoint(steps=1){
    __count('index#callRelaxEndpoint');
    const pos = state.positions.map(p=>[p.x,p.y,p.z]);
  const atomic_numbers = state.elements.map(e=> elementToZ(e));
  const body = { atomic_numbers, coordinates: pos, steps, calculator:'uma' };
  const maybePre = buildPrecomputedIfFresh();
  if(maybePre){ body.precomputed = maybePre; debugApi('relax','precomputed-attach',{ seq: __apiSeq+1, keys:Object.keys(maybePre) }); }
  const uivAtSend = userInteractionVersion; const tivAtSend = totalInteractionVersion;
    const base = __resolveApiBase();
    const url = base + getEndpointSync('relax');
    const payload = body;
    const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
    const build0 = performance.now();
    const t0 = performance.now();
  debugApi('relax','request',{ seq, url, body: payload });
    let resp, json;
    try {
  resp = await fetch(url || '/relax', { method:'POST', headers:{ 'Content-Type':'application/json' }, body: JSON.stringify(payload) });
    } catch(netErr){
      debugApi('relax','network-error',{ seq, url, error: netErr?.message||String(netErr) });
      throw netErr;
    }
    const t1 = performance.now();
    const timingMs = Math.round((t1 - t0)*100)/100; // network time only
    noteLatency('relax','network',{ seq, ms: timingMs });
    if(!resp.ok){
      const txt = await resp.text();
      debugApi('relax','http-error',{ seq, url, timingMs, status: resp.status, statusText: resp.statusText, bodyText: txt });
      throw new Error('Relax request failed '+resp.status+': '+txt);
    }
    let parseMs = 0;
    try { const p0=performance.now(); json = await resp.json(); const p1=performance.now(); parseMs = Math.round((p1-p0)*100)/100; noteLatency('relax','parse',{ seq, ms: parseMs }); }
    catch(parseErr){ debugApi('relax','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); throw parseErr; }
    debugApi('relax','response',{ seq, url, timingMs, status: resp.status, response: json });
  // no server cache key
    return { data: json, uivAtSend, tivAtSend, netMs: timingMs, parseMs };
  }
  const { engine, scene, camera } = await createScene(canvas);
  // Enable verbose touch debug if requested
  try {
    if (typeof window !== 'undefined') {
      window.__MLIPVIEW_DEBUG_TOUCH = !!(window.__MLIPVIEW_DEBUG_TOUCH);
      window.enableTouchDebug = (on)=>{ window.__MLIPVIEW_DEBUG_TOUCH = !!on; console.log('[touchDebug] set to', window.__MLIPVIEW_DEBUG_TOUCH); };
    }
  } catch {}
  const view = createMoleculeView(scene, state);
  const manipulation = createManipulationService(state, { bondService });
  // We'll define wrappedManipulation below; temporarily pass placeholder then rebind after definition
  let wrappedManipulationRef = null;
  // Install custom touch controls BEFORE picking service so it can skip binding its own touch handlers.
  // Support ?noTouch=1 flag to disable touch installation for debugging bond selection issues.
  // Pass a proxy that will read dragActive from the picking service once created.
  let __pickingRef = null;
  const __pickingProxy = { _debug: { get dragActive(){ try { return !!(__pickingRef && __pickingRef._debug && __pickingRef._debug.dragActive); } catch { return false; } } } };
  let __noTouch = false;
  try {
    const q = new URLSearchParams(window.location.search||'');
    __noTouch = (q.get('noTouch') === '1' || q.get('noTouch') === 'true');
  } catch {}
  try {
    if (__noTouch) {
      window.__MLIPVIEW_NO_TOUCH = true;
      console.warn('[touchControls] skipped by ?noTouch=1');
    } else {
      installTouchControls({ canvas, scene, camera, picking: __pickingProxy });
      console.log('[touchControls] installed');
    }
  } catch(e) { console.warn('[touchControls] install failed', e?.message||e); }

  const picking = (__pickingRef = createPickingService(scene, view, selection, {
    manipulation: new Proxy({}, { get: (_, k) => wrappedManipulationRef?.[k] }),
    camera,
    // Do not trigger a standalone simple_calculate while a simulation is running;
    // wait for the next MD/relax response to update forces and energy.
    energyHook: ({ kind }) => {
      try {
        if (!running.kind) ff.computeForces();
      } catch {}
      recordInteraction(kind || 'drag');
    },
  }));
  // Attach VR semantic picker (bond-first) so VR layer can use it without legacy imports
  let vrPicker = null;
  try { vrPicker = createVRPicker({ scene, view }); } catch (e) { console.warn('[VR] vrPicker init failed', e?.message||e); }
  function recomputeBonds(reason='manual') {
    __count('index#recomputeBonds');
    const bonds = bondService.recomputeAndStore();
    try { view.rebuildBonds(bonds); } catch {}
    // Avoid inflating energy step count for internal position change cascades. When
    // markPositionsChanged fires (e.g., during a bond rotation or drag), we already
    // log an interaction specific to the user action (bondRotate, dragMove/dragEnd, etc.).
    // Logging an additional 'rebonds' event creates duplicate energy steps for a single
    // conceptual action. We still log manual or external recomputes so they appear in the
    // trace. If future workflows need explicit bond recompute steps, this guard can be
    // refined (e.g., per-flag) but keeps default UX clean.
    if (reason !== 'markPositionsChanged') {
      recordInteraction('rebonds');
    }
  if(window.__MLIPVIEW_DEBUG_API) console.log(`[bonds] recomputed after ${reason} (count=${bonds?bonds.length:0})`);
    return bonds;
  }
  // Wrap markPositionsChanged now that bondService/view exist
  state.markPositionsChanged = (...a)=>{
    structureVersion++;
    if(state.forceCache) state.forceCache.stale = true;
    const r = __origMarkPositionsChanged? __origMarkPositionsChanged(...a): undefined;
    // Always recompute bonds when positions change
    try { recomputeBonds('markPositionsChanged'); } catch {}
    return r;
  };
  // If no bonds yet, ensure at least a recompute so initial energy uses bond term if applicable
  if (!state.bonds || state.bonds.length===0) {
    try { recomputeBonds(); } catch {}
  }
  async function relaxStep() {
    __count('index#relaxStep');
    try {
  const step0 = performance.now();
  const { data, uivAtSend, tivAtSend, netMs, parseMs } = await callRelaxEndpoint(1); // single step
  const netDone = performance.now();
      // If user interactions occurred during in-flight relax, apply partial update:
      // update only atoms NOT edited by user between send and receive.
      let partialApplied = false;
      let modifiedDuring = null;
      const userChangedDuringFlight = (uivAtSend !== userInteractionVersion);
      if(userChangedDuringFlight){
        modifiedDuring = __modifiedSetBetween(uivAtSend, userInteractionVersion);
      }
      // Build exclusion set: union of modifiedDuring and currently dragged atoms
      const exclude = new Set();
      if(modifiedDuring) for(const i of modifiedDuring) exclude.add(i);
      if(__draggingAtoms.size) for(const i of __draggingAtoms) exclude.add(i);
      // Also exclude any latched atoms whose latch window has not expired
      try {
        const now = (typeof performance!=='undefined' && performance.now) ? performance.now() : Date.now();
        for(const [i, t] of __latchedAtomsUntil){ if(t > now) exclude.add(i); else __latchedAtomsUntil.delete(i); }
      } catch {}
      // If a newer simulation result has already applied AND there were no user edits during flight,
      // mark this response as superseded. If user edits happened, we'll attempt a partial apply below.
      if((tivAtSend !== totalInteractionVersion) && !userChangedDuringFlight){
        if(window.__MLIPVIEW_DEBUG_API) console.debug('[staleStep][relax] supersededSimulation', { uivAtSend, userInteractionVersion, tivAtSend, totalInteractionVersion });
        return { stale:true, staleReason:'supersededSimulation', userInteractionVersionAtSend:uivAtSend, totalInteractionVersionAtSend:tivAtSend, currentUserInteractionVersion:userInteractionVersion, currentTotalInteractionVersion: totalInteractionVersion };
      }
      const { positions, forces, final_energy } = data;
      if(Array.isArray(positions) && positions.length === state.positions.length){
        // If modifiedDuring is set, apply only for indices not in the set
        const N = positions.length;
        const apply0 = performance.now();
        if(exclude && exclude.size>0){
          for(let i=0;i<N;i++){
            if(exclude.has(i)) continue; // keep user's current edit or drag
            const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
          }
          partialApplied = true;
        } else {
          for(let i=0;i<N;i++){
            const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
          }
        }
        state.markPositionsChanged();
        // Suppress the debounced generic posChange tick from counting as a user interaction
        __suppressNextPosChange = true;
        bumpSimulationVersion(partialApplied ? 'relaxStepApplyPartial' : 'relaxStepApply');
        const apply1 = performance.now();
        noteLatency('relax','apply',{ ms: Math.round((apply1-apply0)*100)/100 });
      }
      window.__RELAX_FORCES = forces;
      state.dynamics = state.dynamics || {}; state.dynamics.energy = final_energy;
      if(forces && forces.length){
        let mergedForces = forces;
        if(partialApplied && exclude && exclude.size>0){
          // Merge: keep current forces for user-edited atoms if available
          const cur = Array.isArray(state.forces) ? state.forces : null;
          mergedForces = forces.map((f,i)=> exclude.has(i) ? (cur && cur[i] ? cur[i] : [0,0,0]) : f);
        }
        lastForceResult = { energy: final_energy, forces: mergedForces };
        state.forceCache = { version: structureVersion, energy: final_energy, forces: mergedForces, stress: data.stress||null, stale:false };
        try { __updateForces(mergedForces, { reason: partialApplied ? 'relaxStepPartial' : 'relaxStep' }); } catch {}
      }
      if (partialApplied) {
        const step1 = performance.now();
        noteLatency('relax','step-total',{ ms: Math.round((step1-step0)*100)/100 });
        return { applied:true, partial:true, energy: final_energy, userInteractionVersion, totalInteractionVersion, stepType:'relax', netMs, parseMs };
      }
      { const step1 = performance.now(); noteLatency('relax','step-total',{ ms: Math.round((step1-step0)*100)/100 }); }
      return { applied:true, energy: final_energy, userInteractionVersion, totalInteractionVersion, stepType:'relax', netMs, parseMs };
    } catch(e){
      console.warn('[relaxStep] failed', e);
      return { error: e?.message||String(e) };
    }
  }

  // --- MD step (backend-only logic, but callable from UI) ---
  // Track MD request ordering so concurrent responses can be resolved as latest-wins for application
  let __mdReqCounter = 0;
  let __mdLastApplied = 0;
  async function callMDEndpoint({ steps=1, calculator='uma', temperature=1500, timestep_fs=1.0, friction }={}){
    __count('index#callMDEndpoint');
    const pos = state.positions.map(p=>[p.x,p.y,p.z]);
    const atomic_numbers = state.elements.map(e=> elementToZ(e));
  const frictionFinal = (typeof friction === 'number' && Number.isFinite(friction)) ? friction : (getConfig().mdFriction ?? 0.5);
  const body = { atomic_numbers, coordinates: pos, steps, temperature, timestep_fs, friction: frictionFinal, calculator };
  const maybePre = buildPrecomputedIfFresh();
  if(maybePre){ body.precomputed = maybePre; debugApi('md','precomputed-attach',{ seq: __apiSeq+1, keys:Object.keys(maybePre) }); }
  const maybeV = getVelocitiesIfFresh();
  if(maybeV){ body.velocities = maybeV; }
  const uivAtSend = userInteractionVersion; const tivAtSend = totalInteractionVersion;
    const base = __resolveApiBase();
    const url = base + getEndpointSync('md');
    const payload = body;
    const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
    const build0 = performance.now();
    const t0 = performance.now();
  debugApi('md','request',{ seq, url, body: payload });
  let resp, json; try { resp = await fetch(url,{ method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(payload) }); } catch(netErr){ debugApi('md','network-error',{ seq, url, error: netErr?.message||String(netErr) }); throw netErr; }
    const t1 = performance.now(); const timingMs = Math.round((t1 - t0)*100)/100; // network time only
    noteLatency('md','network',{ seq, ms: timingMs });
    if(!resp.ok){ const txt = await resp.text(); debugApi('md','http-error',{ seq, url, timingMs, status: resp.status, statusText: resp.statusText, bodyText: txt }); throw new Error('MD request failed '+resp.status+': '+txt); }
    let parseMs = 0; try { const p0=performance.now(); json = await resp.json(); const p1=performance.now(); parseMs = Math.round((p1-p0)*100)/100; noteLatency('md','parse',{ seq, ms: parseMs }); } catch(parseErr){ debugApi('md','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); throw parseErr; }
    debugApi('md','response',{ seq, url, timingMs, status: resp.status, response: json });
  // no server cache key
    return { data: json, uivAtSend, tivAtSend, netMs: timingMs, parseMs };
  }
  async function mdStep(opts={}){
    __count('index#mdStep');
    try {
      const t0_step = performance.now();
      const reqId = (++__mdReqCounter);
  const { data, uivAtSend, tivAtSend, netMs, parseMs } = await callMDEndpoint({ steps:1, ...opts });
      // If a newer MD result has already been applied, this one is stale regardless of other checks
      if (reqId < __mdLastApplied) {
        if(window.__MLIPVIEW_DEBUG_API) console.debug('[staleStep][md] superseded by newer applied req', { reqId, lastApplied: __mdLastApplied });
        return { stale:true, staleReason:'supersededSimulation', userInteractionVersionAtSend:uivAtSend, totalInteractionVersionAtSend:tivAtSend, currentUserInteractionVersion:userInteractionVersion, currentTotalInteractionVersion: totalInteractionVersion };
      }
      // If user interactions occurred during in-flight MD, apply partial update for non-edited atoms.
      let partialApplied = false; let modifiedDuring = null;
      if(uivAtSend !== userInteractionVersion){ modifiedDuring = __modifiedSetBetween(uivAtSend, userInteractionVersion); }
      const exclude = new Set();
      if(modifiedDuring) for(const i of modifiedDuring) exclude.add(i);
      if(__draggingAtoms.size) for(const i of __draggingAtoms) exclude.add(i);
      try {
        const now = (typeof performance!=='undefined' && performance.now) ? performance.now() : Date.now();
        for(const [i, t] of __latchedAtomsUntil){ if(t > now) exclude.add(i); else __latchedAtomsUntil.delete(i); }
      } catch {}
  const { positions, forces, final_energy, velocities, temperature: instT } = data;
      if(Array.isArray(positions) && positions.length === state.positions.length){
        const N = positions.length;
        const apply0 = performance.now();
        if(exclude && exclude.size>0){
          for(let i=0;i<N;i++){
            if(exclude.has(i)) continue;
            const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
          }
          partialApplied = true;
        } else {
          for(let i=0;i<N;i++){
            const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
          }
        }
        state.markPositionsChanged();
        __suppressNextPosChange = true;
        bumpSimulationVersion(partialApplied ? 'mdStepApplyPartial' : 'mdStepApply');
        const apply1 = performance.now();
        noteLatency('md','apply',{ ms: Math.round((apply1-apply0)*100)/100 });
      }
      window.__RELAX_FORCES = forces;
      state.dynamics = state.dynamics || {}; state.dynamics.energy = final_energy;
      if(typeof instT === 'number' && isFinite(instT)){
        state.dynamics.temperature = instT;
        try { const el=document.getElementById('instTemp'); if(el) el.textContent = 'T: '+instT.toFixed(1)+' K'; } catch {}
      }
      if(Array.isArray(velocities) && velocities.length === state.elements.length){
        // Shallow validate
        let ok = true; for(let i=0;i<velocities.length && ok;i++){ const r=velocities[i]; if(!r || r.length!==3) ok=false; }
        if(ok){
          if(partialApplied && exclude && exclude.size>0){
            const curV = Array.isArray(state.dynamics?.velocities) ? state.dynamics.velocities : null;
            const mergedV = velocities.map((v,i)=> exclude.has(i) ? (curV && curV[i] ? curV[i] : [0,0,0]) : v);
            state.dynamics.velocities = mergedV;
          } else {
            state.dynamics.velocities = velocities;
          }
        }
      }
      if(forces && forces.length){
        let mergedForces = forces;
        if(partialApplied && exclude && exclude.size>0){
          const cur = Array.isArray(state.forces) ? state.forces : null;
          mergedForces = forces.map((f,i)=> exclude.has(i) ? (cur && cur[i] ? cur[i] : [0,0,0]) : f);
        }
        lastForceResult = { energy: final_energy, forces: mergedForces };
        state.forceCache = { version: structureVersion, energy: final_energy, forces: mergedForces, stress: null, stale:false };
        // Commit forces so visualization updates during MD sequences
        try { __updateForces(mergedForces, { reason: partialApplied ? 'mdStepPartial' : 'mdStep' }); } catch {}
      }
      // Mark this response as the latest applied
      __mdLastApplied = Math.max(__mdLastApplied, reqId);
      if(partialApplied) { noteLatency('md','step-total',{ ms: Math.round((performance.now()-t0_step)*100)/100 }); return { applied:true, partial:true, energy: final_energy, userInteractionVersion, totalInteractionVersion, stepType:'md', netMs, parseMs }; }
      noteLatency('md','step-total',{ ms: Math.round((performance.now()-t0_step)*100)/100 });
      return { applied:true, energy: final_energy, userInteractionVersion, totalInteractionVersion, stepType:'md', netMs, parseMs };
    } catch(e){
      console.warn('[mdStep] failed', e);
      return { error: e?.message||String(e) };
    }
  }

  // Interaction & energy time series
  const interactions = [];
  const energySeries = []; // {i,E,kind}
  let lastPlottedEnergy = undefined;
  function maybePlotEnergy(kind){
    try {
      const E = state.dynamics?.energy;
      if (typeof E !== 'number' || !isFinite(E)) return;
      if (E === lastPlottedEnergy) return; // skip duplicate energy values
      const idx = energySeries.length;
      energySeries.push({ i: idx, E, kind });
      lastPlottedEnergy = E;
      window.__RELAX_TRACE = energySeries.map(e=>e.E);
      drawEnergy();
    } catch(_e) { /* ignore plot errors */ }
  }
  // Suppress the very next debounced generic 'posChange' interaction if a higher-level
  // action (e.g. bondRotate) already recorded an energy step in the same logical edit.
  let __suppressNextPosChange = false;
  // recordInteraction: central funnel for adding user/system events that should reflect
  // a (potential) energy change. Any geometry-changing action (drag move/end, bond rotate,
  // relax step, MD step, provider recompute) should either call recordInteraction directly
  // or be wrapped (see wrappers below). This keeps the energy sparkline in sync without
  // scattering plotting logic across services. We seed an initial 'init' interaction after
  // the first force compute so a single point is drawn (dot) before any steps.
  function recordInteraction(kind){
    __count('index#recordInteraction');
    const E = state.dynamics.energy ?? 0;
    const idx = interactions.length;
    interactions.push({ i:idx, kind, E });
    // Energy plot no longer advances here; only API energy returns call maybePlotEnergy.
    // Notify force visualization (now handled by moleculeView thin instances)
    if (kind && kind !== 'dragMove') {
      try {
        if (state.bus.emit) {
          const DBG = (typeof window !== 'undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search||''));
          if (!state.forces && window.__RELAX_FORCES && window.__RELAX_FORCES.length) {
            __updateForces(window.__RELAX_FORCES, { reason:'lateAttachInteraction' });
          }
          if (DBG) console.log('[Forces][emit] forcesChanged due to interaction kind=', kind, 'haveForces=', !!state.forces);
          state.bus.emit('forcesChanged');
        }
      } catch {}
    }
  }
  async function baselineEnergy(){
    __count('index#baselineEnergy');
    interactions.length = 0;
    energySeries.length = 0; // cleared; first API energy response will add first point
    let res;
    try {
      res = await ff.computeForces({ sync: !!window.__FORCE_SYNC_REMOTE__ });
      state.dynamics = state.dynamics || {}; state.dynamics.energy = res.energy;
    } catch(e){ res = { energy: NaN, forces: [] }; }
    window.__RELAX_FORCES = res.forces||[];
    if (res.forces && res.forces.length) {
      __updateForces(res.forces, { reason:'baselineEnergy' });
    } else {
      const DBG=(typeof window!=='undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search||''));
      if (DBG) {
        const synth = state.positions.map(p=>{ const r=Math.hypot(p.x,p.y,p.z)||1; return [p.x/r*0.3, p.y/r*0.3, p.z/r*0.3]; });
        window.__RELAX_FORCES = synth;
        __updateForces(synth, { reason:'baselineEnergySynth' });
      }
    }
    // Do not create an 'init' energy tick. First real API response will create the first point.
  }
  // Wrap atomic-changing operations (relax/MD steps already wrapped below)
  const origRelaxStep = relaxStep; relaxStep = async ()=>{ const r = await origRelaxStep(); if(r && r.applied){ recordInteraction('relaxStep'); maybePlotEnergy('relax'); } return r; };
  const origMdStep = mdStep; mdStep = async (o)=>{ const r = await origMdStep(o); if(r && r.applied){ recordInteraction('mdStep'); maybePlotEnergy('md'); } return r; };
  // Wrap manipulation (drag & bond rotation) so any geometry change recomputes energy and updates plot.
  function __selectedAtomIndex(){ try { return (state.selection && state.selection.kind==='atom') ? state.selection.data.index : null; } catch { return null; } }
  const wrappedManipulation = {
    beginDrag: (...a)=> { const ok = manipulation.beginDrag(...a); if(ok){ const idx = __selectedAtomIndex(); if(idx!=null){ __currentDraggedAtomIndex = idx; __draggingAtoms.add(idx); } } return ok; },
    updateDrag: (...a)=> { const r = manipulation.updateDrag(...a); if (r) { bumpUserInteractionVersion('dragMove'); const idx = __selectedAtomIndex(); if(idx!=null){ __addModifiedAtoms([idx]); __currentDraggedAtomIndex = idx; __draggingAtoms.add(idx); } if(!running.kind) ff.computeForces(); recordInteraction('dragMove'); } return r; },
    endDrag: (...a)=> {
      const hadDrag = (__currentDraggedAtomIndex!=null) || (__draggingAtoms.size>0);
      const r = manipulation.endDrag(...a);
      if(hadDrag){
        bumpUserInteractionVersion('dragEnd');
        const idx = (__currentDraggedAtomIndex!=null)? __currentDraggedAtomIndex : __selectedAtomIndex();
        if(idx!=null){
          __addModifiedAtoms([idx]);
          __draggingAtoms.delete(idx);
          // Latch this atom briefly so any in-flight MD/relax response doesn’t snap it
          __latchAtom(idx, 800);
        }
        // Clear current pointer
        __currentDraggedAtomIndex = null;
      }
      if(!running.kind) ff.computeForces();
      recordInteraction('dragEnd');
      return r;
    },
    setDragPlane: (...a)=> manipulation.setDragPlane(...a),
    rotateBond: (...a)=> { const r = manipulation.rotateBond(...a); if (r) { // capture rotated side atoms via manipulation debug hook
        let sideAtoms = null; try { sideAtoms = manipulation._debug?.getLastRotation?.()?.sideAtoms || null; } catch{}
        bumpUserInteractionVersion('bondRotate'); if(Array.isArray(sideAtoms) && sideAtoms.length) __addModifiedAtoms(sideAtoms); if(!running.kind) ff.computeForces(); recordInteraction('bondRotate'); __suppressNextPosChange = true; }
      return r; }
  };
  // Mark for debug/diagnostics so other subsystems (VR) can assert they received the wrapped instance
  try { Object.defineProperty(wrappedManipulation, '__isWrappedManipulation', { value: true, enumerable: false }); } catch { wrappedManipulation.__isWrappedManipulation = true; }
  wrappedManipulationRef = wrappedManipulation;

  // Provide explicit cleanup to help Jest teardown
  if (typeof window !== 'undefined') {
    window.__MLIPVIEW_CLEANUP = window.__MLIPVIEW_CLEANUP || [];
    window.__MLIPVIEW_CLEANUP.push(()=>{
      try { engine && engine.stopRenderLoop && engine.stopRenderLoop(); } catch {}
      try { if(scene && scene.dispose) scene.dispose(); } catch {}
    });
  }
  // Selection or manipulation events could be hooked similarly via bus later

  let energyCtx=null; let energyCanvas=null; let energyLabel=null;
  function initEnergyCanvas(){
    __count('index#initEnergyCanvas');
    if (energyCanvas) return;
    energyCanvas = document.getElementById('energyCanvas');
    if (!energyCanvas) return;
    energyCtx = energyCanvas.getContext('2d');
    energyLabel = document.getElementById('energyLabel');
  }
  function drawEnergy(){
    __count('index#drawEnergy');
    initEnergyCanvas(); if(!energyCtx) return;
    const W=energyCanvas.width, H=energyCanvas.height;
    energyCtx.clearRect(0,0,W,H);
    if (energySeries.length===0) return;
    if (energySeries.length===1) {
      const p=energySeries[0];
      energyCtx.fillStyle='#6fc2ff';
      energyCtx.beginPath();
      energyCtx.arc(W/2, H/2, 3, 0, Math.PI*2);
      energyCtx.fill();
      if (energyLabel) energyLabel.textContent = 'E steps=1';
      return;
    }
    let minE=Infinity,maxE=-Infinity; for(const p of energySeries){ if(p.E<minE)minE=p.E; if(p.E>maxE)maxE=p.E; }
    if (maxE-minE < 1e-12) { maxE=minE+1e-12; }
    energyCtx.strokeStyle='#6fc2ff'; energyCtx.lineWidth=1; energyCtx.beginPath();
    for (let k=0;k<energySeries.length;k++){
      const p=energySeries[k];
      const x = (k/(energySeries.length-1))* (W-4) + 2;
      const y = H - 2 - ((p.E - minE)/(maxE-minE))*(H-4);
      if(k===0) energyCtx.moveTo(x,y); else energyCtx.lineTo(x,y);
    }
    energyCtx.stroke();
    if (energyLabel) energyLabel.textContent = `E steps=${energySeries.length} range=${(maxE-minE).toExponential(2)}`;
  }

  // Debounced auto energy update when any positionsChanged event fires (catch-all for drag paths)
  let posEnergyTimer=null; let pendingPosEnergy=false;
  state.bus.on('positionsChanged', () => {
    pendingPosEnergy = true;
    if (posEnergyTimer) return;
    posEnergyTimer = setTimeout(()=>{
      posEnergyTimer=null;
      if (!pendingPosEnergy) return;
      pendingPosEnergy=false;
      // If this positionsChanged was triggered by applying a simulation step (relax/md),
      // skip the immediate force recompute so we don't overwrite the step's returned forces.
      if (__suppressNextPosChange) {
        __suppressNextPosChange = false; // consume the suppression flag
        // Intentionally do NOT call ff.computeForces or record an interaction here.
        return;
      }
      try { if(!running.kind) ff.computeForces(); } catch{}
      bumpUserInteractionVersion('posChange');
      recordInteraction('posChange');
    }, 50); // slight debounce to batch rapid pointer move events
  });

  // --- Continuous simulation orchestration (relax / MD) ---
  // Feature flags (can be toggled at build time via define plugin):
  // Feature flags are now read dynamically each access so they can be enabled after init.
  // Previous implementation captured a snapshot causing UI buttons to remain disabled if flags
  // were set after viewer initialization. Use getter to always reflect current window.__MLIP_FEATURES.
  function featureFlags(){
    if (typeof window === 'undefined') return { RELAX_LOOP:true, MD_LOOP:true, ENERGY_TRACE:true, FORCE_VECTORS:true };
    if(!window.__MLIP_FEATURES) window.__MLIP_FEATURES = { RELAX_LOOP:true, MD_LOOP:true, ENERGY_TRACE:true, FORCE_VECTORS:true };
    return window.__MLIP_FEATURES;
  }
  function enableFeatureFlag(name, value=true){
    if (typeof window === 'undefined') return false;
    if(!window.__MLIP_FEATURES) window.__MLIP_FEATURES={};
    window.__MLIP_FEATURES[name]=value; return true;
  }
  let running = { kind: null, abort: null };
  function setForceProvider(){ return 'uma'; }
  // Continuous loops with request pacing and exponential backoff on server overload.
  // Contract:
  //  - Make at most one network step request per configured pacing interval (minStepIntervalMs, default 30ms) in normal operation.
  //  - Adjustable at runtime via setMinStepInterval(ms).
  //  - If a step returns an error (http or parse caught inside relaxStep/mdStep), apply
  //    exponential backoff starting at 200ms (200,400,800,... up to 5s) before retrying next step.
  //  - Abort if another simulation is already running or stopSimulation() called.
  //  - For relax: attempt up to maxSteps (default 200) or until an error streak exceeds threshold.
  //  - For md: run fixed number of steps unless aborted.
  // Desktop focus helpers: only send relax/MD requests when window is focused (desktop mode)
  function isAppFocused(){
    try {
      if (typeof window === 'undefined' || typeof document === 'undefined') return true; // non-browser/test
      if (window.__MLIPVIEW_TEST_MODE) return true; // do not gate tests
      if (document.hidden) return false;
      if (typeof document.hasFocus === 'function') return !!document.hasFocus();
      return true;
    } catch { return true; }
  }
  function waitForFocus(){
    return new Promise(resolve => {
      if (isAppFocused()) return resolve();
      function cleanup(){ try { window.removeEventListener('focus', onFocus); document.removeEventListener('visibilitychange', onVis); } catch {} }
      function onFocus(){ cleanup(); resolve(); }
      function onVis(){ try { if (!document.hidden && isAppFocused()) { cleanup(); resolve(); } } catch { /* ignore */ } }
      try { window.addEventListener('focus', onFocus); document.addEventListener('visibilitychange', onVis); } catch { resolve(); }
    });
  }

  async function startRelaxContinuous({ maxSteps=1000 }={}) {
  if(!featureFlags().RELAX_LOOP) { console.warn('[feature] RELAX_LOOP disabled'); return { disabled:true }; }
    __count('index#startRelaxContinuous');
    if (running.kind) return { ignored:true };
    running.kind='relax';
    __resetRPS();
  const minInterval = getConfig().minStepIntervalMs; // ms pacing per API request (configurable)
    let lastTime=0; // retained for focus-wait reset only
    let backoffMs=0; // 0 means no backoff active
    const baseBackoff=200; const maxBackoff=5000;
    let errorStreak=0; const maxErrorStreak=10;
    let stepsDone=0;
  while(running.kind==='relax' && stepsDone<maxSteps){
      // Pause if app not focused (desktop). Resume on focus.
      if (!isAppFocused()) {
        const f0 = performance.now();
        await waitForFocus();
        if (running.kind!=='relax') break;
        // reset pacing reference to avoid burst immediately after focus
        lastTime = performance.now();
        noteLatency('relax','focus-wait',{ ms: Math.round((lastTime - f0)*100)/100 });
        continue;
      }
      if(backoffMs>0){
        const w0 = performance.now();
        await new Promise(r=>setTimeout(r, backoffMs));
        noteLatency('relax','wait-backoff',{ ms: Math.round((performance.now()-w0)*100)/100, backoffMs });
      }
      if(running.kind!=='relax') break;
  const res = await relaxStep();
  if(res && res.applied){ stepsDone++; __noteRequestCompleted(); } // count only applied (non-stale) steps
      lastTime = performance.now();
      // Post-request pacing: only wait if API network time was less than minInterval
      if(!res?.error && backoffMs===0){
        const net = Number(res?.netMs)||0;
        const need = Math.max(0, minInterval - net);
        if(need>0){ const w0=performance.now(); await new Promise(r=>setTimeout(r, need)); noteLatency('relax','wait-pacing',{ ms: Math.round((performance.now()-w0)*100)/100, needed: Math.round(need), netMs: net }); }
      }
      if(res && res.error){
        errorStreak++;
        if(errorStreak>=maxErrorStreak){ console.warn('[relaxRun] aborting due to error streak'); break; }
        backoffMs = backoffMs? Math.min(backoffMs*2, maxBackoff) : baseBackoff;
      } else {
        errorStreak=0; backoffMs=0;
      }
    }
    const converged = errorStreak===0 && stepsDone>=maxSteps; // simplistic criterion
    running.kind=null;
    __resetRPS();
    return { converged, steps: stepsDone };
  }
  async function startMDContinuous({ steps=1000, calculator='uma', temperature=1500, timestep_fs=1.0, friction }={}){
  if(!featureFlags().MD_LOOP){ console.warn('[feature] MD_LOOP disabled'); return { disabled:true }; }
    if(running.kind) return { ignored:true };
    running.kind='md';
    __resetRPS();
  const minInterval = getConfig().minStepIntervalMs; let lastTime=0; // retained for focus-wait reset only
    let backoffMs=0; const baseBackoff=200; const maxBackoff=5000; let errorStreak=0; const maxErrorStreak=10;
    let i=0;
    while(i<steps && running.kind==='md'){
      if (!isAppFocused()) {
        const f0 = performance.now();
        await waitForFocus();
        if (running.kind!=='md') break;
        lastTime = performance.now();
        noteLatency('md','focus-wait',{ ms: Math.round((lastTime - f0)*100)/100 });
        continue;
      }
  if(backoffMs>0){ const w0=performance.now(); await new Promise(r=>setTimeout(r, backoffMs)); noteLatency('md','wait-backoff',{ ms: Math.round((performance.now()-w0)*100)/100, backoffMs }); }
  if(running.kind!=='md') break;
      // Re-read target temperature dynamically to allow live slider adjustments mid-run.
      let dynT = temperature; // fallback to initial argument
      try {
        if (typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null) {
          dynT = window.__MLIP_TARGET_TEMPERATURE;
        }
      } catch {}
  // Clamp to configured slider bounds (0-3000 K)
      if(!(Number.isFinite(dynT))) dynT = temperature;
  if(dynT < 0) dynT = 0; else if(dynT > 3000) dynT = 3000;
  const res = await mdStep({ calculator, temperature: dynT, timestep_fs, friction });
  try { const el=document.getElementById('instTemp'); if(el && state.dynamics && typeof state.dynamics.temperature==='number') el.textContent='T: '+state.dynamics.temperature.toFixed(1)+' K'; } catch {}
      lastTime=performance.now();
      // Post-request pacing: only wait if API network time was less than minInterval
      if(!res?.error && backoffMs===0){
        const net = Number(res?.netMs)||0;
        const need = Math.max(0, minInterval - net);
        if(need>0){ const w0=performance.now(); await new Promise(r=>setTimeout(r, need)); noteLatency('md','wait-pacing',{ ms: Math.round((performance.now()-w0)*100)/100, needed: Math.round(need), netMs: net }); }
      }
      if(res && res.error){
        errorStreak++; if(errorStreak>=maxErrorStreak){ console.warn('[mdRun] aborting due to error streak'); break; }
        backoffMs = backoffMs? Math.min(backoffMs*2, maxBackoff) : baseBackoff;
      } else {
        // Only advance step counter for applied (non-stale) MD results
        if(res && res.applied){ i++; __noteRequestCompleted(); }
        errorStreak=0; backoffMs=0;
      }
    }
    const completed = (i>=steps && errorStreak===0);
    running.kind=null;
    __resetRPS();
    return { completed, steps: i };
  }
  function stopSimulation(){ running.kind=null; __resetRPS(); }

  const lastMetrics = { energy:null, maxForce:null, maxStress:null };
  function getMetrics(){ __count('index#getMetrics'); return { energy: state.dynamics?.energy, running: running.kind }; }

  // Initial energy baseline (compute once to seed plot so first interaction can draw a segment)
  try { await baselineEnergy(); } catch(e) { /* ignore */ }

  // --- Auto-start MD (optional) ---
  // Requirement: Start continuous MD automatically after the very first successful energy/force
  // acquisition (baselineEnergy above). We only run this in normal browser mode (not test mode)
  // and allow users/tests to opt out via window.__MLIPVIEW_NO_AUTO_MD = true before init.
  // We trigger via the public API so UI state (run/stop button text) can reflect the run.
  try {
    const autoOk = (typeof window !== 'undefined') && !window.__MLIPVIEW_TEST_MODE && !window.__MLIPVIEW_NO_AUTO_MD;
    if (autoOk) {
      // Defer a tick to allow index.html script (buttons & handlers) to finish wiring viewerApi
      setTimeout(()=>{
        try {
          // Avoid starting if another loop already active
            if(!window.viewerApi) return; // safety
            const m = window.viewerApi.getMetrics();
            if(m.running) return; // already running something
            // Start MD with default parameters; UI interval/metrics loop will update status label.
            window.viewerApi.startMDContinuous({}).then(()=>{
              try {
                const btn = document.getElementById && document.getElementById('btnMDRun');
                if(btn && btn.textContent === 'stop') btn.textContent='run';
                const statusEl = document.getElementById && document.getElementById('status');
                if(statusEl && /MD running/.test(statusEl.textContent||'')) statusEl.textContent='MD stopped';
              } catch{}
              if (window.__MLIPVIEW_DEBUG_API) console.log('[autoMD] completed initial MD run');
            });
            // Update button text immediately if present to mimic user click path.
            try {
              const btn = document.getElementById && document.getElementById('btnMDRun');
              if(btn){ btn.textContent = 'stop'; }
              const statusEl = document.getElementById && document.getElementById('status');
              if(statusEl){ statusEl.textContent = 'MD running'; }
            } catch {}
        } catch(e){ console.warn('[autoMD] start failed', e?.message||e); }
      }, 0);
    }
  } catch(_e){ /* ignore auto start errors */ }

  // --- Render loop (critical for Babylon WebXR) ---
  // NOTE: WebXR integration in Babylon expects engine.runRenderLoop to own the frame pump so it can
  // swap to XR's requestAnimationFrame internally. Replacing it with a raw requestAnimationFrame can
  // cause enterXRAsync to hang or never resolve. Re-introduce runRenderLoop as default, keeping a
  // test-mode fallback for environments (like unit tests / jsdom) that lack a proper RAF or canvas.
  let __renderActive = true;
  const __testMode = (typeof window !== 'undefined') && !!window.__MLIPVIEW_TEST_MODE;
  function __rafLoop(){
    if(!__renderActive) return;
    try { scene.render(); } catch(e){ console.warn('[Render][raf] render error', e); }
    if(typeof requestAnimationFrame === 'function') requestAnimationFrame(__rafLoop);
  }
  try {
    if(!__testMode && engine?.runRenderLoop){
      console.log('[Render] starting engine.runRenderLoop (XR compatible)');
      engine.runRenderLoop(()=>{
        if(!__renderActive){ return; }
        try { scene.render(); } catch(e){ console.warn('[Render] loop error', e); }
      });
    } else {
      console.log('[Render] starting requestAnimationFrame loop (test mode fallback)');
      if(typeof requestAnimationFrame === 'function') requestAnimationFrame(__rafLoop);
      else if(engine?.runRenderLoop){ engine.runRenderLoop(()=>{ if(!__renderActive) return; scene.render(); }); }
    }
  } catch(e){
    console.warn('[Render] primary loop init failed; falling back to rAF', e);
    if(typeof requestAnimationFrame === 'function') requestAnimationFrame(__rafLoop);
  }

  // Force vectors now rendered via moleculeView thin instances (see moleculeView.rebuildForces).
  // Keep a small compatibility layer: setting true will ensure showForces=true, false hides; undefined toggles.
  function setForceVectorsEnabled(on){
    try {
      if (typeof on === 'boolean') { if (!!state.showForces !== on) state.toggleForceVectorsVisibility(); }
      else { state.toggleForceVectorsVisibility(); }
    } catch {}
  }

  // VR support is lazy; user can call vr.init() explicitly later.
  // IMPORTANT: Route VR interactions through the wrapped manipulation so dragged atoms are tracked
  // and excluded from MD/relax position application while held.
  const vr = createVRSupport(scene, { picking: { ...picking, view, vrPicker, selectionService: selection, manipulation: new Proxy({}, { get: (_, k) => wrappedManipulationRef?.[k] }), molState: state } });
  // Auto-init VR support so controllers & debug logging are ready; actual immersive session still requires user gesture.
  try {
    vr.init().then(res => { if (res.supported) { console.log('[VR] support initialized (auto)'); } else { console.log('[VR] not supported'); } });
  } catch (e) { console.warn('[VR] auto init failed', e); }
  function debugEnergySeriesLength(){ return energySeries.length; }
  function debugRecordInteraction(kind){ recordInteraction(kind||'debug'); }
  function getForceCacheVersion(){ return state.forceCache?.version; }
  function getVersionInfo(){ return { userInteractionVersion, totalInteractionVersion }; }
  function shutdown(){ __renderActive=false; try{ engine.stopRenderLoop && engine.stopRenderLoop(); }catch{} }
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, startRelaxContinuous, startMDContinuous, stopSimulation, setForceProvider, getMetrics, debugEnergySeriesLength, debugRecordInteraction, manipulation: wrappedManipulation, scene, engine, camera, baselineEnergy, setForceVectorsEnabled, getForceCacheVersion, getVersionInfo, shutdown, enableFeatureFlag, setMinStepInterval };
}

// Ensure global viewerApi for tests if not already set
if (typeof window !== 'undefined') {
  try {
    if (!window.viewerApi && typeof initNewViewer === 'function') {
      // Will be assigned after user calls initNewViewer; patch the function to set it.
      const __origInit = initNewViewer;
      // Redefine only once
      Object.defineProperty(window, 'initNewViewer', { value: async function(...args){
        const api = await __origInit(...args);
        window.viewerApi = api; // expose globally for tests & debug inspector
        return api;
      }, configurable: true });
    }
  } catch(_) { /* ignore */ }
}

// Debug / test helpers injected after definition (non-enumerable minimal surface impact)
Object.defineProperty(window, '__dumpCurrentAtoms', { value: function(){
  try {
    if(!window.viewerApi) return null;
    const st = window.viewerApi.state;
    return {
      elements: st.elements.map(e=> e.symbol||e.sym||e.S||e.Z||e.atomicNumber||'?'),
      atomic_numbers: st.elements.map(e=> (typeof e==='string'? (SYMBOL_TO_Z[e]||0) : (e && (e.Z||e.atomicNumber||e.z|| (SYMBOL_TO_Z[e.symbol]||SYMBOL_TO_Z[e.sym]||SYMBOL_TO_Z[e.S]||0))) || 0)),
      positions: st.positions.map(p=> [p.x,p.y,p.z]),
      energy: st.dynamics?.energy
    };
  } catch(e){ return { error: e?.message||String(e) }; }
}, writable:false });
