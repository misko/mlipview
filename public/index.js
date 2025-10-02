import { createMoleculeState } from './domain/moleculeState.js';
import { getEndpointSync } from './api_endpoints.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { createFairChemForcefield } from './fairchem_provider.js';
import { createScene } from './render/scene.js';
import { createMoleculeView } from './render/moleculeView.js';
import { createPickingService } from './core/pickingService.js';
import { createManipulationService } from './domain/manipulationService.js';
import { createVRSupport } from './vr/setup.js';
import { createVRPicker } from './vr/vr-picker.js';
import { __count } from './util/funcCount.js';

// --- Runtime Config (pacing, etc.) ---
if (typeof window !== 'undefined') {
  window.__MLIP_CONFIG = window.__MLIP_CONFIG || { minStepIntervalMs: 30 };
}
function getConfig(){ if (typeof window === 'undefined') return { minStepIntervalMs:30 }; return window.__MLIP_CONFIG; }
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
    // default ON until user disables explicitly
    window.__MLIPVIEW_DEBUG_API = true;
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
  async function fetchRemoteForces({ awaitResult=false }={}){
    __count('index#fetchRemoteForces');
    if(inFlight && !awaitResult) return;
    while(inFlight && awaitResult) { await new Promise(r=>setTimeout(r,10)); }
    inFlight = true;
    try {
  // If cache is fresh (not stale and matches current structureVersion) return it.
  if(!state.forceCache.stale && state.forceCache.version === structureVersion){
        lastForceResult = { energy: state.forceCache.energy, forces: state.forceCache.forces };
        return lastForceResult;
      }
  // Skip remote call until we actually have atoms; prevents 500 errors on empty initial state
  if(!state.elements || state.elements.length === 0){
        return lastForceResult;
      }
  const atomic_numbers = state.elements.map(e=> elementToZ(e));
      const coordinates = state.positions.map(p=>[p.x,p.y,p.z]);
      const body = { atomic_numbers, coordinates, properties:['energy','forces'], calculator:'uma' };
    const base = __resolveApiBase();
    const ep = getEndpointSync('simple');
    const url = base + ep;
      const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
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
      if(resp.ok){
        try { json = await resp.json(); } catch(parseErr){ debugApi('simple_calculate','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); }
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
        }
      }
    } catch(_e) { /* silent for now */ } finally { inFlight=false; }
    return lastForceResult;
  }
  const ff = { computeForces: ({ sync }={})=>{ if(sync) return fetchRemoteForces({ awaitResult:true }); if(!inFlight) fetchRemoteForces(); return lastForceResult; } };
  const dynamics = { stepMD: ()=>{}, stepRelax: ({ forceFn })=>{ forceFn(); } };
  // Remote relaxation: call backend /relax
  async function callRelaxEndpoint(steps=1){
    __count('index#callRelaxEndpoint');
    const pos = state.positions.map(p=>[p.x,p.y,p.z]);
  const atomic_numbers = state.elements.map(e=> elementToZ(e));
    const body = { atomic_numbers, coordinates: pos, steps, calculator:'uma' };
    const base = __resolveApiBase();
    const ep = getEndpointSync('relax');
    const url = base + ep;
    const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
    const t0 = performance.now();
    debugApi('relax','request',{ seq, url, body });
    let resp, json;
    try {
      resp = await fetch(url || '/relax', { method:'POST', headers:{ 'Content-Type':'application/json' }, body: JSON.stringify(body) });
    } catch(netErr){
      debugApi('relax','network-error',{ seq, url, error: netErr?.message||String(netErr) });
      throw netErr;
    }
    const t1 = performance.now();
    const timingMs = Math.round((t1 - t0)*100)/100;
    if(!resp.ok){
      const txt = await resp.text();
      debugApi('relax','http-error',{ seq, url, timingMs, status: resp.status, statusText: resp.statusText, bodyText: txt });
      throw new Error('Relax request failed '+resp.status+': '+txt);
    }
    try { json = await resp.json(); } catch(parseErr){ debugApi('relax','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); throw parseErr; }
    debugApi('relax','response',{ seq, url, timingMs, status: resp.status, response: json });
    return json;
  }
  const { engine, scene, camera } = await createScene(canvas);
  const view = createMoleculeView(scene, state);
  const manipulation = createManipulationService(state, { bondService });
  // We'll define wrappedManipulation below; temporarily pass placeholder then rebind after definition
  let wrappedManipulationRef = null;
  const picking = createPickingService(scene, view, selection, { manipulation: new Proxy({}, { get:(_,k)=> wrappedManipulationRef?.[k] }), camera, energyHook: ({ kind }) => { try { ff.computeForces(); } catch{} recordInteraction(kind||'drag'); } });
  // Attach VR semantic picker (bond-first) so VR layer can use it without legacy imports
  let vrPicker = null;
  try { vrPicker = createVRPicker({ scene, view }); } catch (e) { console.warn('[VR] vrPicker init failed', e?.message||e); }
  function recomputeBonds(reason='manual') {
    __count('index#recomputeBonds');
    const bonds = bondService.recomputeAndStore();
    try { view.rebuildBonds(bonds); } catch {}
    recordInteraction('rebonds');
    console.log(`[bonds] recomputed after ${reason} (count=${bonds?bonds.length:0})`);
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
      const data = await callRelaxEndpoint(1); // single step
  const { positions, forces, final_energy } = data;
      // Apply positions
      if(Array.isArray(positions) && positions.length === state.positions.length){
        for(let i=0;i<positions.length;i++){
          const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2]; }
        state.markPositionsChanged(); // central wrapper handles bonds & logging
      }
      window.__RELAX_FORCES = forces;
      state.dynamics = state.dynamics || {}; state.dynamics.energy = final_energy;
      if(forces && forces.length){
        lastForceResult = { energy: final_energy, forces };
        state.forceCache = { version: structureVersion, energy: final_energy, forces, stress: data.stress||null, stale:false };
      }
      return { energy: final_energy };
    } catch(e){
      console.warn('[relaxStep] failed', e);
      return { error: e?.message||String(e) };
    }
  }

  // --- MD step (backend-only logic, but callable from UI) ---
  async function callMDEndpoint({ steps=1, calculator='uma', temperature=298, timestep_fs=1.0, friction=0.02 }={}){
    __count('index#callMDEndpoint');
    const pos = state.positions.map(p=>[p.x,p.y,p.z]);
    const atomic_numbers = state.elements.map(e=> elementToZ(e));
    const body = { atomic_numbers, coordinates: pos, steps, temperature, timestep_fs, friction, calculator };
    const base = __resolveApiBase();
    const ep = getEndpointSync('md');
    const url = base + ep;
    const seq = ++__apiSeq; window.__MLIPVIEW_API_SEQ = __apiSeq;
    const t0 = performance.now();
    debugApi('md','request',{ seq, url, body });
    let resp, json; try { resp = await fetch(url,{ method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body) }); } catch(netErr){ debugApi('md','network-error',{ seq, url, error: netErr?.message||String(netErr) }); throw netErr; }
    const t1 = performance.now(); const timingMs = Math.round((t1 - t0)*100)/100;
    if(!resp.ok){ const txt = await resp.text(); debugApi('md','http-error',{ seq, url, timingMs, status: resp.status, statusText: resp.statusText, bodyText: txt }); throw new Error('MD request failed '+resp.status+': '+txt); }
    try { json = await resp.json(); } catch(parseErr){ debugApi('md','parse-error',{ seq, url, timingMs, error: parseErr?.message||String(parseErr) }); throw parseErr; }
    debugApi('md','response',{ seq, url, timingMs, status: resp.status, response: json });
    return json;
  }
  async function mdStep(opts={}){
    __count('index#mdStep');
    try {
      const data = await callMDEndpoint({ steps:1, ...opts });
      const { positions, forces, final_energy } = data;
      if(Array.isArray(positions) && positions.length === state.positions.length){
        for(let i=0;i<positions.length;i++){
          const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2]; }
        state.markPositionsChanged(); // central wrapper handles bonds & logging
      }
      window.__RELAX_FORCES = forces;
      state.dynamics = state.dynamics || {}; state.dynamics.energy = final_energy;
      if(forces && forces.length){
        lastForceResult = { energy: final_energy, forces };
        state.forceCache = { version: structureVersion, energy: final_energy, forces, stress: null, stale:false };
      }
      return { energy: final_energy };
    } catch(e){
      console.warn('[mdStep] failed', e);
      return { error: e?.message||String(e) };
    }
  }

  // Interaction & energy time series
  const interactions = [];
  const energySeries = []; // {i, E}
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
    if (typeof E === 'number') energySeries.push({ i:idx, E });
    // Expose for test harness
    window.__RELAX_TRACE = energySeries.map(e=>e.E);
    drawEnergy();
    // Update force vectors on meaningful geometry/energy events (skip high-frequency dragMove handled separately)
    if (forceVis.enabled && kind && kind !== 'dragMove') {
      try { updateForceVectors(); } catch {}
    }
  }
  async function baselineEnergy(){
    __count('index#baselineEnergy');
    interactions.length = 0;
    energySeries.length = 0;
    let res;
    try {
      res = await ff.computeForces({ sync: !!window.__FORCE_SYNC_REMOTE__ });
      state.dynamics = state.dynamics || {}; state.dynamics.energy = res.energy;
    } catch(e){ res = { energy: NaN, forces: [] }; }
    window.__RELAX_FORCES = res.forces||[];
    recordInteraction('init');
    if (forceVis.enabled) { try { updateForceVectors(); } catch {} }
  }
  // Wrap atomic-changing operations (relax/MD steps already wrapped below)
  const origRelaxStep = relaxStep; relaxStep = async ()=>{ recordInteraction('relaxStep:pending'); const r = await origRelaxStep(); recordInteraction('relaxStep'); return r; };
  const origMdStep = mdStep; mdStep = async (o)=>{ recordInteraction('mdStep:pending'); const r = await origMdStep(o); recordInteraction('mdStep'); return r; };
  // Wrap manipulation (drag & bond rotation) so any geometry change recomputes energy and updates plot.
  const wrappedManipulation = {
    beginDrag: (...a)=> manipulation.beginDrag(...a),
    updateDrag: (...a)=> { const r = manipulation.updateDrag(...a); if (r) { ff.computeForces(); recordInteraction('dragMove'); } return r; },
    endDrag: (...a)=> { const r = manipulation.endDrag(...a); structureVersion++; if(state.forceCache) state.forceCache.stale=true; ff.computeForces(); recordInteraction('dragEnd'); return r; },
    setDragPlane: (...a)=> manipulation.setDragPlane(...a),
    rotateBond: (...a)=> { const r = manipulation.rotateBond(...a); if (r) { structureVersion++; if(state.forceCache) state.forceCache.stale=true; ff.computeForces(); recordInteraction('bondRotate'); } return r; }
  };
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
      try { ff.computeForces(); } catch{}
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
  async function startRelaxContinuous({ maxSteps=200 }={}) {
  if(!featureFlags().RELAX_LOOP) { console.warn('[feature] RELAX_LOOP disabled'); return { disabled:true }; }
    __count('index#startRelaxContinuous');
    if (running.kind) return { ignored:true };
    running.kind='relax';
  const minInterval = getConfig().minStepIntervalMs; // ms pacing between successful requests (configurable)
    let lastTime=0;
    let backoffMs=0; // 0 means no backoff active
    const baseBackoff=200; const maxBackoff=5000;
    let errorStreak=0; const maxErrorStreak=10;
    let stepsDone=0;
    while(running.kind==='relax' && stepsDone<maxSteps){
      const now=performance.now();
      const since= now - lastTime;
      if(backoffMs>0){
        await new Promise(r=>setTimeout(r, backoffMs));
      } else if (since < minInterval){
        await new Promise(r=>setTimeout(r, minInterval - since));
      }
      if(running.kind!=='relax') break;
      const res = await relaxStep();
      stepsDone++;
      lastTime = performance.now();
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
    return { converged, steps: stepsDone };
  }
  async function startMDContinuous({ steps=200, calculator='uma', temperature=298, timestep_fs=1.0, friction=0.02 }={}){
  if(!featureFlags().MD_LOOP){ console.warn('[feature] MD_LOOP disabled'); return { disabled:true }; }
    if(running.kind) return { ignored:true };
    running.kind='md';
  const minInterval = getConfig().minStepIntervalMs; let lastTime=0;
    let backoffMs=0; const baseBackoff=200; const maxBackoff=5000; let errorStreak=0; const maxErrorStreak=10;
    let i=0;
    while(i<steps && running.kind==='md'){
      const now=performance.now(); const since=now-lastTime;
      if(backoffMs>0){ await new Promise(r=>setTimeout(r, backoffMs)); }
      else if(since<minInterval){ await new Promise(r=>setTimeout(r, minInterval - since)); }
      if(running.kind!=='md') break;
      const res = await mdStep({ calculator, temperature, timestep_fs, friction });
      lastTime=performance.now();
      if(res && res.error){
        errorStreak++; if(errorStreak>=maxErrorStreak){ console.warn('[mdRun] aborting due to error streak'); break; }
        backoffMs = backoffMs? Math.min(backoffMs*2, maxBackoff) : baseBackoff;
      } else { errorStreak=0; backoffMs=0; i++; }
    }
    const completed = (i>=steps && errorStreak===0);
    running.kind=null;
    return { completed, steps: i };
  }
  function stopSimulation(){ running.kind=null; }

  const lastMetrics = { energy:null, maxForce:null, maxStress:null };
  function getMetrics(){ __count('index#getMetrics'); return { energy: state.dynamics?.energy, running: running.kind }; }

  // Initial energy baseline (compute once to seed plot so first interaction can draw a segment)
  try { await baselineEnergy(); } catch(e) { /* ignore */ }

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

  // --- Force vector visualization (legacy per-arrow implementation) ---
  // Force vectors: always enabled (requirement). We'll auto-refresh on positionsChanged and each frame.
  const forceVis = { enabled: true, root:null, arrows:[], maxLen:1.0 };
  function ensureForceRoot(){
    __count('index#ensureForceRoot');
    if(!forceVis.root){ forceVis.root = new BABYLON.TransformNode('forceArrowsRoot', scene); }
    // If molecule top-level transform exists (e.g., an atom master mesh parent) attach to it so arrows inherit VR rotations.
    if(forceVis.root && !forceVis.root.parent){
      const anyAtom = scene.meshes.find(m=>m && /^atom_/i.test(m.name));
      if(anyAtom && anyAtom.parent){ forceVis.root.parent = anyAtom.parent; }
    }
    return forceVis.root;
  }
  const forceMaterial = new BABYLON.StandardMaterial('forceMat', scene);
  forceMaterial.diffuseColor = new BABYLON.Color3(0.9,0.1,0.1);
  forceMaterial.emissiveColor = new BABYLON.Color3(0.6,0.05,0.05);
  forceMaterial.specularColor = new BABYLON.Color3(0.2,0.2,0.2);
  forceMaterial.disableLighting = false;
  function disposeArrows(){ __count('index#disposeArrows'); if(forceVis.arrows){ for(const a of forceVis.arrows){ try{ a.dispose(); }catch{} } } forceVis.arrows=[]; }
  function createArrow(name){
    __count('index#createArrow');
    const body = BABYLON.MeshBuilder.CreateCylinder(name+'_body',{height:1,diameter:0.04,tessellation:12},scene);
    const tip  = BABYLON.MeshBuilder.CreateCylinder(name+'_tip',{height:0.25,diameterTop:0,diameterBottom:0.12,tessellation:12},scene);
    tip.position.y = 0.5 + 0.125; tip.parent = body;
    body.material = forceMaterial; tip.material = forceMaterial;
    body.isPickable=false; tip.isPickable=false; body.visibility=1; return body;
  }
  function updateForceVectors(){
    __count('index#updateForceVectors');
    if(!forceVis.enabled) return;
    const forces = window.__RELAX_FORCES || [];
    if(!forces.length || !state.positions.length) return;
    const root = ensureForceRoot();
    if(!root.parent){
      const sample = scene.meshes.find(m=>m && m.name && /^atom_/i.test(m.name));
      if(sample && sample.parent){ root.parent = sample.parent; }
    }
    if(forceVis.arrows.length !== forces.length){
      disposeArrows();
      for(let i=0;i<forces.length;i++){ const a=createArrow('forceArrow'+i); a.parent = root; forceVis.arrows.push(a); }
    }
    // Compute max magnitude for scaling
    let maxMag=0; const mags=new Array(forces.length);
    for(let i=0;i<forces.length;i++){ const f=forces[i]; const mag=Math.hypot(f[0],f[1],f[2]); mags[i]=mag; if(mag>maxMag) maxMag=mag; }
    const baseScale = maxMag>0 ? 0.6/maxMag : 0;
    for(let i=0;i<forces.length;i++){
      const f=forces[i]; const mag=mags[i]; const p=state.positions[i]; const arrow=forceVis.arrows[i];
      if(!p || mag<1e-8){ arrow.setEnabled(false); continue; } else { arrow.setEnabled(true); }
      const dir = new BABYLON.Vector3(f[0],f[1],f[2]); dir.normalize();
      const up = BABYLON.Vector3.Up(); let q; const dot=BABYLON.Vector3.Dot(up,dir);
      if(Math.abs(dot-1)<1e-6) q = BABYLON.Quaternion.Identity();
      else if(Math.abs(dot+1)<1e-6) q = BABYLON.Quaternion.RotationAxis(BABYLON.Vector3.Right(), Math.PI);
      else { const axis = BABYLON.Vector3.Cross(up,dir).normalize(); const angle=Math.acos(dot); q = BABYLON.Quaternion.RotationAxis(axis, angle); }
      arrow.rotationQuaternion = q;
      const length = mag*baseScale; const clampedLen = Math.min(Math.max(length,0.05),1.2);
      arrow.scaling = new BABYLON.Vector3(1, clampedLen, 1);
      const offset = new BABYLON.Vector3(0, clampedLen/2, 0);
      const rotatedOffset = offset.rotateByQuaternionToRef(q, new BABYLON.Vector3());
      arrow.position.set(p.x + rotatedOffset.x, p.y + rotatedOffset.y, p.z + rotatedOffset.z);
    }
    root.setEnabled(forceVis.enabled);
  }
  function setForceVectorsEnabled(on){
    __count('index#setForceVectorsEnabled');
    // Always-on policy now; ignore external off requests but keep API for backward compatibility
    forceVis.enabled = true;
    updateForceVectors();
  }
  // schedule initial draw & continuous updates
  setTimeout(()=>{ try{ updateForceVectors(); }catch{} },0);
  // Refresh arrows when atom positions mutate
  state.bus.on('positionsChanged', () => { try { updateForceVectors(); } catch {} });
  // Per-frame (lightweight) orientation sync in case root parenting occurs after dynamic scene changes
  if(scene.onBeforeRenderObservable){
    scene.onBeforeRenderObservable.add(()=>{ if(forceVis.enabled) { try { updateForceVectors(); } catch {} } });
  }

  // VR support is lazy; user can call vr.init() explicitly later.
  const vr = createVRSupport(scene, { picking: { ...picking, view, vrPicker, selectionService: selection, manipulation, molState: state } });
  // Auto-init VR support so controllers & debug logging are ready; actual immersive session still requires user gesture.
  try {
    vr.init().then(res => { if (res.supported) { console.log('[VR] support initialized (auto)'); } else { console.log('[VR] not supported'); } });
  } catch (e) { console.warn('[VR] auto init failed', e); }
  function debugEnergySeriesLength(){ return energySeries.length; }
  function debugRecordInteraction(kind){ recordInteraction(kind||'debug'); }
  function getForceCacheVersion(){ return state.forceCache?.version; }
  function shutdown(){ __renderActive=false; try{ engine.stopRenderLoop && engine.stopRenderLoop(); }catch{} }
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, startRelaxContinuous, startMDContinuous, stopSimulation, setForceProvider, getMetrics, debugEnergySeriesLength, debugRecordInteraction, manipulation: wrappedManipulation, scene, engine, camera, baselineEnergy, setForceVectorsEnabled, getForceCacheVersion, shutdown, enableFeatureFlag, setMinStepInterval };
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
