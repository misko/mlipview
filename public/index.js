import { createMoleculeState } from './domain/moleculeState.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { ljEnergyForces, createBFGSStepper } from './lj_bfgs.js';
import { createFairChemForcefield } from './fairchem_provider.js';
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
  // Minimal dynamics placeholder
  let providerKind = 'lj';
  let ff = { computeForces: ()=>{ const pos=state.positions.map(p=>[p.x,p.y,p.z]); const { energy, forces } = ljEnergyForces(pos); state.dynamics = state.dynamics||{}; state.dynamics.energy = energy; return { energy, forces }; } };
  const dynamics = { stepMD: ()=>{}, stepRelax: ({ forceFn })=>{ forceFn(); } };
  // BFGS stepper (lazy)
  let bfgsStepper = null;
  function ensureStepper(){
    if(!bfgsStepper){
      const pos = state.positions.map(p=>[p.x,p.y,p.z]);
      const compute = async p => {
        if (providerKind === 'fairchem') {
          // Evaluate FairChem energy/forces at trial positions p without mutating state until accepted.
          return await ff.computeForces(p);
        }
        return ljEnergyForces(p);
      };
      bfgsStepper = createBFGSStepper({ positions: pos, fmax:0.05, maxStep:0.2, compute });
      for(let i=0;i<state.positions.length;i++){ state.positions[i].x=pos[i][0]; state.positions[i].y=pos[i][1]; state.positions[i].z=pos[i][2]; }
    }
    return bfgsStepper;
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
    const stepper = ensureStepper();
  const res = await stepper.step();
  const cur = await stepper.getCurrent();
  // Commit energy for accepted geometry (FairChem state energy committed here rather than in compute)
  state.dynamics = state.dynamics || {}; state.dynamics.energy = cur.energy;
    window.__RELAX_FORCES = cur.forces;
    // Sync updated numeric positions from stepper back into state objects
    if (stepper.positions && stepper.positions.length === state.positions.length) {
      for (let i=0;i<state.positions.length;i++) {
        const sp = stepper.positions[i];
        const tp = state.positions[i];
        if (tp.x !== sp[0] || tp.y !== sp[1] || tp.z !== sp[2]) {
          tp.x = sp[0]; tp.y = sp[1]; tp.z = sp[2];
        }
      }
      // Emit position change so view updates atom thin instances and highlights
      state.markPositionsChanged();
      // Bonds: if bond cylinders depend on positions only (length/orientation) we can trigger a lightweight rebuild
      // Instead of full recompute of connectivity, just refresh transform buffers
      try { view.rebuildBonds(); } catch {}
    }
    return { energy: cur.energy, fmax: cur.fmax, converged: res?.converged };
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
      if (providerKind === 'fairchem') {
        const pos = state.positions.map(p=>[p.x,p.y,p.z]);
        // Use provider with explicit positions; commit energy to state
        res = await ff.computeForces(pos);
        state.dynamics = state.dynamics || {}; state.dynamics.energy = res.energy;
      } else {
        res = ff.computeForces();
        state.dynamics = state.dynamics || {}; state.dynamics.energy = res.energy;
      }
    } catch(e){ res = { energy: NaN, forces: [] }; }
    window.__RELAX_FORCES = res.forces;
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
  function setForceProvider(kind){
    if (kind === providerKind) return providerKind;
    providerKind = kind === 'fairchem' ? 'fairchem' : 'lj';
    // Replace forcefield
    if (providerKind === 'fairchem') {
      ff = createFairChemForcefield(state);
    } else {
      ff = { computeForces: ()=>{ const pos=state.positions.map(p=>[p.x,p.y,p.z]); const { energy, forces } = ljEnergyForces(pos); state.dynamics = state.dynamics||{}; state.dynamics.energy = energy; return { energy, forces }; } };
    }
    // Reset stepper so it uses new compute
    bfgsStepper = null;
  try { baselineEnergy(); } catch {}
    return providerKind;
  }
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

  // --- Force vector visualization ---
  const forceVis = { enabled: true, arrows: [], root: null, maxLen: 1.0 };
  function ensureForceRoot(){
    if (!forceVis.root){
      forceVis.root = new BABYLON.TransformNode('forceArrowsRoot', scene);
    }
    return forceVis.root;
  }
  function disposeArrows(){
    if (forceVis.arrows){ for (const a of forceVis.arrows){ try { a.parent = null; a.dispose(); } catch{} } }
    forceVis.arrows = [];
  }
  function createArrow(baseName){
    // Body
    const body = BABYLON.MeshBuilder.CreateCylinder(baseName+'_body', { height:1, diameter:0.04, tessellation:12 }, scene);
    // Tip
    const tip = BABYLON.MeshBuilder.CreateCylinder(baseName+'_tip', { height:0.25, diameterTop:0, diameterBottom:0.12, tessellation:12 }, scene);
    tip.position.y = 0.5 + 0.125; // body half + tip half
    tip.parent = body;
    body.material = forceMaterial;
    tip.material = forceMaterial;
    body.isPickable = false; tip.isPickable = false;
    return body;
  }
  const forceMaterial = new BABYLON.StandardMaterial('forceMat', scene);
  forceMaterial.diffuseColor = new BABYLON.Color3(0.9,0.1,0.1);
  forceMaterial.emissiveColor = new BABYLON.Color3(0.6,0.05,0.05);
  forceMaterial.specularColor = new BABYLON.Color3(0.2,0.2,0.2);
  forceMaterial.disableLighting = false;
  function updateForceVectors(){
    if (!forceVis.enabled) return;
    const forces = window.__RELAX_FORCES || [];
    if (!forces.length || !state.positions.length) return;
    const root = ensureForceRoot();
    // Ensure arrow count matches atom count
    if (forceVis.arrows.length !== forces.length){
      disposeArrows();
      for (let i=0;i<forces.length;i++){
        const arrow = createArrow('fArrow'+i);
        arrow.parent = root;
        forceVis.arrows.push(arrow);
      }
    }
    // Determine scale (magnitude encoded directly, clamp overly large values)
    let maxMag = 0;
    const mags = new Array(forces.length);
    for (let i=0;i<forces.length;i++){
      const f = forces[i];
      const mag = Math.hypot(f[0],f[1],f[2]);
      mags[i] = mag;
      if (mag>maxMag) maxMag = mag;
    }
    // Avoid zero division; choose a base length so moderate forces are visible
    const baseScale = maxMag > 0 ? 0.6 / maxMag : 0.0; // largest arrow body ~0.6 units
    for (let i=0;i<forces.length;i++){
      const f = forces[i];
      const mag = mags[i];
      const arrow = forceVis.arrows[i];
      const p = state.positions[i];
      // Force direction: arrows should point along the physical force vector (f[0],f[1],f[2])
      // If displayed opposite, invert here (current forces appear to point atom->atom, likely correct sign already?)
      const dir = new BABYLON.Vector3(f[0], f[1], f[2]);
      // If near zero magnitude, hide
      if (mag < 1e-8){ arrow.setEnabled(false); continue; } else { arrow.setEnabled(true); }
      // Normalize direction
      dir.normalize();
      // Position arrow so its base sits at atom center; cylinder created along +Y, so orient via quaternion
      // Compute quaternion from up (0,1,0) to dir
      const up = BABYLON.Vector3.Up();
      let q;
      if (Math.abs(BABYLON.Vector3.Dot(up, dir) - 1) < 1e-6) {
        q = BABYLON.Quaternion.Identity();
      } else if (Math.abs(BABYLON.Vector3.Dot(up, dir) + 1) < 1e-6) {
        // Opposite direction
        q = new BABYLON.Quaternion();
        BABYLON.Quaternion.RotationAxisToRef(BABYLON.Vector3.Right(), Math.PI, q);
      } else {
        const axis = BABYLON.Vector3.Cross(up, dir);
        axis.normalize();
        const angle = Math.acos(BABYLON.Vector3.Dot(up, dir));
        q = BABYLON.Quaternion.RotationAxis(axis, angle);
      }
      arrow.rotationQuaternion = q;
      // Scale: body height originally 1; tip extends additional ~0.25; scale y to reflect magnitude
      const length = mag * baseScale; // encoded magnitude
      // Prevent absurdly tiny or huge visuals
      const clampedLen = Math.min(Math.max(length, 0.05), 1.2);
      arrow.scaling = new BABYLON.Vector3(1, clampedLen, 1);
      // Move arrow so base at atom: body height is scaled; origin of cylinder is center, so shift by scaledHeight/2 along +Y then rotate
      // We'll compute offset before rotation: (0, clampedLen/2, 0)
      const offset = new BABYLON.Vector3(0, clampedLen/2, 0);
      const rotatedOffset = offset.rotateByQuaternionToRef(q, new BABYLON.Vector3());
      arrow.position.set(p.x + rotatedOffset.x, p.y + rotatedOffset.y, p.z + rotatedOffset.z);
    }
    root.setEnabled(forceVis.enabled);
  }
  function setForceVectorsEnabled(on){
    forceVis.enabled = !!on;
    if (!forceVis.enabled){ if (forceVis.root) forceVis.root.setEnabled(false); }
    else { updateForceVectors(); }
  }

  // Initial attempt to render after baseline
  setTimeout(()=>{ if(forceVis.enabled) try { updateForceVectors(); } catch{} }, 0);

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
