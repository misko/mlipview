import { createMoleculeState } from './domain/moleculeState.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { createForceField } from './physics/forcefield.js';
import { createForceProvider } from './physics/force-provider.js';
import { createRelaxationRunner } from './physics/relaxation-runner.js';
import { createDynamics } from './physics/integrators.js';
import { createScene } from './render/scene.js';
import { createMoleculeView } from './render/moleculeView.js';
import { createPickingService } from './core/pickingService.js';
import { createManipulationService } from './domain/manipulationService.js';
import { createVRSupport } from './vr/setup.js';
import { createVRPicker } from './vr/vr-picker.js';

export async function initNewViewer(canvas, { elements, positions, bonds } ) {
  const state = createMoleculeState({ elements, positions, bonds });
  const bondService = createBondService(state);
  const selection = createSelectionService(state);
  const ff = createForceField(state, {});
  const dynamics = createDynamics(state, {});
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
  function relaxStep() { dynamics.stepRelax({ forceFn: ff.computeForces }); }
  function mdStep() { dynamics.stepMD({ forceFn: ff.computeForces }); }

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
    drawEnergy();
  }
  // Wrap atomic-changing operations (relax/MD steps already wrapped below)
  const origRelaxStep = relaxStep; relaxStep = ()=>{ origRelaxStep(); recordInteraction('relaxStep'); };
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
  let activeProviderKind = 'local';
  const providers = { local: createForceProvider('local', { molState: state }) };
  function ensureFairChemProvider(baseUrl){
    if (!providers.fairchem) providers.fairchem = createForceProvider('fairchem', { baseUrl });
    return providers.fairchem;
  }
  function setForceProvider(kind, options={}) {
    if (kind === 'fairchem') ensureFairChemProvider(options.baseUrl||'http://127.0.0.1:8000');
    if (!providers[kind]) throw new Error('Unknown provider '+kind);
    activeProviderKind = kind;
    return activeProviderKind;
  }
  // Simple async loop wrappers
  async function startRelaxContinuous({ forceTol=1e-3, stressTol=5e-3 }={}) {
    if (running.kind) return; // already running
    running.kind = 'relax';
    const provider = providers[activeProviderKind];
    // Create lightweight simulation proxy referencing real state arrays
    const proxy = {
      Z: Int32Array.from(state.elements.map(e=> (typeof e==='number'? e : (e==='C'?6:e==='H'?1:6)))),
      get pos(){ // expose live Float64Array view
        if (!this._pos || this._pos.length !== state.positions.length*3) {
          this._pos = new Float64Array(state.positions.length*3);
          for (let i=0;i<state.positions.length;i++){ const p=state.positions[i]; this._pos[3*i]=p.x; this._pos[3*i+1]=p.y; this._pos[3*i+2]=p.z; }
        }
        return this._pos;
      },
      set pos(v){ this._pos = v; },
      get box(){ return state.cell && state.cell.enabled ? { a:state.cell.a, b:state.cell.b, c:state.cell.c } : null; }
    };
    const runner = createRelaxationRunner(proxy, { forceProvider: provider, settings:{ forceTol, stressTol, onUpdate: (u)=>{
      // sync proxy positions back into real state
      for (let i=0;i<state.positions.length;i++){ const p=state.positions[i]; p.x = proxy.pos[3*i]; p.y = proxy.pos[3*i+1]; p.z = proxy.pos[3*i+2]; }
      state.markPositionsChanged();
      if (typeof u.maxForce === 'number') lastMetrics.maxForce = u.maxForce;
      if (typeof u.energy === 'number') lastMetrics.energy = u.energy;
      if (typeof u.maxStress === 'number' && u.maxStress != null) lastMetrics.maxStress = u.maxStress;
      if (settings?.onUpdate) settings.onUpdate(u);
    } } });
    running.abort = runner.abort;
    const res = await runner.run();
    recordInteraction('relaxRunEnd');
    running.kind = null; running.abort = null;
    return res;
  }
  async function startMDContinuous({ dt=0.5, targetTemp=300 }={}) {
    if (running.kind) return;
    running.kind = 'md';
    let stopped = false;
    running.abort = ()=>{ stopped = true; running.kind=null; running.abort=null; };
    while(!stopped) {
      dynamics.stepMD({ forceFn: ff.computeForces, targetTemp });
      state.markPositionsChanged();
      recordInteraction('mdLoop');
      await new Promise(r=>setTimeout(r, 16)); // ~60Hz
    }
  }
  function stopSimulation(){ if (running.abort) running.abort(); }

  const lastMetrics = { energy:null, maxForce:null, maxStress:null };
  function getMetrics(){ return { ...lastMetrics, running: running.kind }; }

  // Initial energy baseline (compute once to seed plot so first interaction can draw a segment)
  try { ff.computeForces(); recordInteraction('init'); } catch(e) { /* ignore */ }

  engine.runRenderLoop(()=>{ scene.render(); });
  // VR support is lazy; user can call vr.init() explicitly later.
  const vr = createVRSupport(scene, { picking: { ...picking, view, vrPicker, selectionService: selection, manipulation, molState: state } });
  // Auto-init VR support so controllers & debug logging are ready; actual immersive session still requires user gesture.
  try {
    vr.init().then(res => { if (res.supported) { console.log('[VR] support initialized (auto)'); } else { console.log('[VR] not supported'); } });
  } catch (e) { console.warn('[VR] auto init failed', e); }
  function debugEnergySeriesLength(){ return energySeries.length; }
  function debugRecordInteraction(kind){ recordInteraction(kind||'debug'); }
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, startRelaxContinuous, startMDContinuous, stopSimulation, setForceProvider, getMetrics, debugEnergySeriesLength, debugRecordInteraction, manipulation: wrappedManipulation, scene, engine, camera };
}
