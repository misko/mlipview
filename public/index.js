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
  const state = createMoleculeState({ elements, positions, bonds });
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

  async function fetchRemoteForces({ awaitResult=false }={}){
    if(inFlight && !awaitResult) return;
    while(inFlight && awaitResult) { await new Promise(r=>setTimeout(r,10)); }
    inFlight = true;
    try {
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
        }
      }
    } catch(_e) { /* silent for now */ } finally { inFlight=false; }
    return lastForceResult;
  }
  const ff = { computeForces: ({ sync }={})=>{ if(sync) return fetchRemoteForces({ awaitResult:true }); if(!inFlight) fetchRemoteForces(); return lastForceResult; } };
  const dynamics = { stepMD: ()=>{}, stepRelax: ({ forceFn })=>{ forceFn(); } };
  // Remote relaxation: call backend /relax
  async function callRelaxEndpoint(steps=1){
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
    const bonds = bondService.recomputeAndStore();
    view.rebuildBonds(bonds);
    recordInteraction('rebonds');
  }
  // If no bonds yet, ensure at least a recompute so initial energy uses bond term if applicable
  if (!state.bonds || state.bonds.length===0) {
    try { recomputeBonds(); } catch {}
  }
  async function relaxStep() {
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
      return { energy: final_energy };
    } catch(e){
      console.warn('[relaxStep] failed', e);
      return { error: e?.message||String(e) };
    }
  }
  function mdStep() { /* no-op in minimal build */ }

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
  const origRelaxStep = relaxStep; relaxStep = async ()=>{ await origRelaxStep(); recordInteraction('relaxStep'); };
  const origMdStep = mdStep; mdStep = ()=>{ origMdStep(); recordInteraction('mdStep'); };
  // Wrap manipulation (drag & bond rotation) so any geometry change recomputes energy and updates plot.
  const wrappedManipulation = {
    beginDrag: (...a)=> manipulation.beginDrag(...a),
    updateDrag: (...a)=> { const r = manipulation.updateDrag(...a); if (r) { ff.computeForces(); recordInteraction('dragMove'); } return r; },
    endDrag: (...a)=> { const r = manipulation.endDrag(...a); ff.computeForces(); recordInteraction('dragEnd'); return r; },
    setDragPlane: (...a)=> manipulation.setDragPlane(...a),
    rotateBond: (...a)=> { const r = manipulation.rotateBond(...a); if (r) { ff.computeForces(); recordInteraction('bondRotate'); } return r; }
  };
  wrappedManipulationRef = wrappedManipulation;
  // Selection or manipulation events could be hooked similarly via bus later

  let energyCtx=null; let energyCanvas=null; let energyLabel=null;
  function initEnergyCanvas(){
    if (energyCanvas) return;
    energyCanvas = document.getElementById('energyCanvas');
    if (!energyCanvas) return;
    energyCtx = energyCanvas.getContext('2d');
    energyLabel = document.getElementById('energyLabel');
  }
  function drawEnergy(){
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
  let running = { kind: null, abort: null };
  function setForceProvider(){ return 'uma'; }
  // Simple async loop wrappers
  async function startRelaxContinuous() {
    if (running.kind) return;
    running.kind='relax';
  for (let i=0;i<200;i++){ await relaxStep(); recordInteraction('relaxStep'); }
    running.kind=null;
    return { converged:true };
  }
  async function startMDContinuous(){ /* no-op */ }
  function stopSimulation(){ running.kind=null; }

  const lastMetrics = { energy:null, maxForce:null, maxStress:null };
  function getMetrics(){ return { energy: state.dynamics?.energy, running: running.kind }; }

  // Initial energy baseline (compute once to seed plot so first interaction can draw a segment)
  try { await baselineEnergy(); } catch(e) { /* ignore */ }

  engine.runRenderLoop(()=>{ scene.render(); });

  // --- Force vector visualization (legacy per-arrow implementation) ---
  const forceVis = { enabled:true, root:null, arrows:[], maxLen:1.0 };
  function ensureForceRoot(){
    if(!forceVis.root){ forceVis.root = new BABYLON.TransformNode('forceArrowsRoot', scene); }
    return forceVis.root;
  }
  const forceMaterial = new BABYLON.StandardMaterial('forceMat', scene);
  forceMaterial.diffuseColor = new BABYLON.Color3(0.9,0.1,0.1);
  forceMaterial.emissiveColor = new BABYLON.Color3(0.6,0.05,0.05);
  forceMaterial.specularColor = new BABYLON.Color3(0.2,0.2,0.2);
  forceMaterial.disableLighting = false;
  function disposeArrows(){ if(forceVis.arrows){ for(const a of forceVis.arrows){ try{ a.dispose(); }catch{} } } forceVis.arrows=[]; }
  function createArrow(name){
    const body = BABYLON.MeshBuilder.CreateCylinder(name+'_body',{height:1,diameter:0.04,tessellation:12},scene);
    const tip  = BABYLON.MeshBuilder.CreateCylinder(name+'_tip',{height:0.25,diameterTop:0,diameterBottom:0.12,tessellation:12},scene);
    tip.position.y = 0.5 + 0.125; tip.parent = body;
    body.material = forceMaterial; tip.material = forceMaterial;
    body.isPickable=false; tip.isPickable=false; body.visibility=1; return body;
  }
  function updateForceVectors(){
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
    forceVis.enabled = !!on;
    if(!forceVis.enabled){ if(forceVis.root) forceVis.root.setEnabled(false); }
    else { updateForceVectors(); }
  }
  // schedule initial draw
  setTimeout(()=>{ if(forceVis.enabled) try{ updateForceVectors(); }catch{} },0);

  // VR support is lazy; user can call vr.init() explicitly later.
  const vr = createVRSupport(scene, { picking: { ...picking, view, vrPicker, selectionService: selection, manipulation, molState: state } });
  // Auto-init VR support so controllers & debug logging are ready; actual immersive session still requires user gesture.
  try {
    vr.init().then(res => { if (res.supported) { console.log('[VR] support initialized (auto)'); } else { console.log('[VR] not supported'); } });
  } catch (e) { console.warn('[VR] auto init failed', e); }
  function debugEnergySeriesLength(){ return energySeries.length; }
  function debugRecordInteraction(kind){ recordInteraction(kind||'debug'); }
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, startRelaxContinuous, startMDContinuous, stopSimulation, setForceProvider, getMetrics, debugEnergySeriesLength, debugRecordInteraction, manipulation: wrappedManipulation, scene, engine, camera, baselineEnergy, setForceVectorsEnabled };
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
