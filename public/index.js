import { createMoleculeState } from './domain/moleculeState.js';
import { createBondService } from './domain/bondService.js';
import { createSelectionService } from './domain/selectionService.js';
import { createForceField } from './physics/forcefield.js';
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
  const picking = createPickingService(scene, view, selection, { manipulation, camera });
  // Attach VR semantic picker (bond-first) so VR layer can use it without legacy imports
  let vrPicker = null;
  try { vrPicker = createVRPicker({ scene, view }); } catch (e) { console.warn('[VR] vrPicker init failed', e?.message||e); }
  function recomputeBonds() {
    const bonds = bondService.recomputeAndStore();
    view.rebuildBonds(bonds);
  }
  function relaxStep() { dynamics.stepRelax({ forceFn: ff.computeForces }); }
  function mdStep() { dynamics.stepMD({ forceFn: ff.computeForces }); }
  engine.runRenderLoop(()=>{ scene.render(); });
  // VR support is lazy; user can call vr.init() explicitly later.
  const vr = createVRSupport(scene, { picking: { ...picking, view, vrPicker, selectionService: selection, manipulation, molState: state } });
  // Auto-init VR support so controllers & debug logging are ready; actual immersive session still requires user gesture.
  try {
    vr.init().then(res => { if (res.supported) { console.log('[VR] support initialized (auto)'); } else { console.log('[VR] not supported'); } });
  } catch (e) { console.warn('[VR] auto init failed', e); }
  return { state, bondService, selection, ff, dynamics, view, picking, vr, recomputeBonds, relaxStep, mdStep, manipulation, scene, engine, camera };
}
