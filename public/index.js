import { createMoleculeState } from './domain/moleculeState.js';
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

export async function initNewViewer(canvas, { elements, positions, bonds } ) {
  __count('index#initNewViewer');
  const state = createMoleculeState({ elements, positions, bonds });
  // Wrap markPositionsChanged to invalidate force cache deterministically
  const __origMarkPositionsChanged = state.markPositionsChanged ? state.markPositionsChanged.bind(state) : null;
  state.markPositionsChanged = (...a)=>{ structureVersion++; state.forceCache && (state.forceCache.stale = true); if(__origMarkPositionsChanged) return __origMarkPositionsChanged(...a); };
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
      const base = (window.__MLIPVIEW_SERVER||'').replace(/\/$/,'');
      const url = (base? base : '') + '/simple_calculate'; // relative when base empty so Vite proxy handles it
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
    const base = (window.__MLIPVIEW_SERVER || '').replace(/\/$/, '');
    const url = (base? base : '') + '/relax';
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
  function recomputeBonds() {
    __count('index#recomputeBonds');
    const bonds = bondService.recomputeAndStore();
    view.rebuildBonds(bonds);
    recordInteraction('rebonds');
    // Bond topology change invalidates force cache
    structureVersion++; if(state.forceCache) state.forceCache.stale = true;
  }
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
          const p=positions[i]; const tp=state.positions[i];
          tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
        }
        state.markPositionsChanged();
        try { view.rebuildBonds(); } catch {}
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
    const base = (window.__MLIPVIEW_SERVER || '').replace(/\/$/,'');
    const url = (base? base : '') + '/md';
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
          const p=positions[i]; const tp=state.positions[i]; tp.x=p[0]; tp.y=p[1]; tp.z=p[2];
        }
        state.markPositionsChanged();
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
  const origRelaxStep = relaxStep; relaxStep = async ()=>{ const r = await origRelaxStep(); recordInteraction('relaxStep'); return r; };
  const origMdStep = mdStep; mdStep = async (o)=>{ const r = await origMdStep(o); recordInteraction('mdStep'); return r; };
  // Wrap manipulation (drag & bond rotation) so any geometry change recomputes energy and updates plot.
  const wrappedManipulation = {
    beginDrag: (...a)=> manipulation.beginDrag(...a),
    updateDrag: (...a)=> { const r = manipulation.updateDrag(...a); if (r) { ff.computeForces(); recordInteraction('dragMove'); } return r; },
    endDrag: (...a)=> { const r = manipulation.endDrag(...a); structureVersion++; if(state.forceCache) state.forceCache.stale=true; ff.computeForces(); recordInteraction('dragEnd'); return r; },
    setDragPlane: (...a)=> manipulation.setDragPlane(...a),
    rotateBond: (...a)=> { const r = manipulation.rotateBond(...a); if (r) { structureVersion++; if(state.forceCache) state.forceCache.stale=true; ff.computeForces(); recordInteraction('bondRotate'); } return r; }
  };
  wrappedManipulationRef = wrappedManipulation;
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
  const FEATURES = (typeof window !== 'undefined' && window.__MLIP_FEATURES) || { RELAX_LOOP:false, MD_LOOP:false, ENERGY_TRACE:true, FORCE_VECTORS:true };
  let running = { kind: null, abort: null };
  function setForceProvider(){ return 'uma'; }
  async function startRelaxContinuous() {
    if(!FEATURES.RELAX_LOOP) { console.warn('[feature] RELAX_LOOP disabled'); return { disabled:true }; }
    __count('index#startRelaxContinuous');
    if (running.kind) return;
    running.kind='relax';
    for (let i=0;i<200;i++){ await relaxStep(); recordInteraction('relaxStep'); }
    running.kind=null;
    return { converged:true };
  }
  async function startMDContinuous({ steps=200, delayMs=0, calculator='uma', temperature=298, timestep_fs=1.0, friction=0.02 }={}){
    if(!FEATURES.MD_LOOP){ console.warn('[feature] MD_LOOP disabled'); return { disabled:true }; }
    if(running.kind) return;
    running.kind='md';
    for(let i=0;i<steps && running.kind==='md'; i++){
      await mdStep({ calculator, temperature, timestep_fs, friction });
      if(delayMs>0 && running.kind==='md') await new Promise(r=>setTimeout(r,delayMs));
    }
    const finished = running.kind==='md';
    running.kind=null;
    return { completed: finished };
  }
  function stopSimulation(){ running.kind=null; }

  const lastMetrics = { energy:null, maxForce:null, maxStress:null };
  function getMetrics(){ __count('index#getMetrics'); return { energy: state.dynamics?.energy, running: running.kind }; }

  // Initial energy baseline (compute once to seed plot so first interaction can draw a segment)
  try { await baselineEnergy(); } catch(e) { /* ignore */ }

  let __renderActive = true;
  function __renderLoop(){ if(!__renderActive) return; scene.render(); requestAnimationFrame(__renderLoop); }
  // Use requestAnimationFrame to allow graceful shutdown in tests
  if (typeof requestAnimationFrame === 'function') { requestAnimationFrame(__renderLoop); } else { engine.runRenderLoop(()=>{ scene.render(); }); }

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
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, startRelaxContinuous, startMDContinuous, stopSimulation, setForceProvider, getMetrics, debugEnergySeriesLength, debugRecordInteraction, manipulation: wrappedManipulation, scene, engine, camera, baselineEnergy, setForceVectorsEnabled, getForceCacheVersion, shutdown };
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
