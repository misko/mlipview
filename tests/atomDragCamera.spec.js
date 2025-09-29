import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.js';
import { createPickingService } from '../public/core/pickingService.js';

// BABYLON stubs
if (!global.BABYLON) global.BABYLON = {}; // test runner likely already sets
if (!BABYLON.PointerEventTypes) BABYLON.PointerEventTypes = { POINTERDOWN:1 };

function makeSceneStub(pickAtom){
  return {
    pointerX:50,
    pointerY:50,
    pick(){ return pickAtom ? { hit:true, pickedMesh:{}, thinInstanceIndex:0 } : { hit:false }; },
    onPointerObservable:{ add(fn){ this._cb=fn; } },
    getEngine(){ return { getRenderingCanvas(){ return { addEventListener(){} }; } }; }
  };
}

describe('atom drag does not move camera', () => {
  test('camera params unchanged while atom moves', () => {
    const st = createMoleculeState({ elements:['C','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0}], bonds:[{i:0,j:1}] });
    // selection service stub
    const selection = { get:()=>st.selection, clickAtom:(i)=>{ st.selection={ kind:'atom', data:{ index:i } }; st.markSelectionChanged(); }, clickBond:()=>{}, clear:()=>{ st.selection={ kind:null, data:null }; st.markSelectionChanged(); } };
    const bondService = { recomputeAndStore: ()=>{} };
    const manip = createManipulationService(st, { bondService });
    const cameraState = { alpha:1, beta:2, radius:10, target:{ x:0, y:0, z:0 } };
    const camera = {
      detachControl(){ /* simulate detaching inputs only */ },
      attachControl(){ /* simulate reattaching */ },
      get alpha(){ return cameraState.alpha; },
      get beta(){ return cameraState.beta; },
      get radius(){ return cameraState.radius; },
      get target(){ return cameraState.target; }
    };

    const scene = makeSceneStub(true);
    const view = { resolveAtomPick:()=>({ kind:'atom', index:0 }), resolveBondPick:()=>null };
    createPickingService(scene, view, selection, { manipulation: manip, camera });
    // Trigger pointer down (starts drag)
    scene.onPointerObservable._cb({ type: BABYLON.PointerEventTypes.POINTERDOWN });

    // Simulate drag updates
  // Provide hit point that shifts atom in negative X (current simplistic intersector subtracts pointer move)
  // The manipulation service computes grabOffset = startPos - firstHit.
  // We simulated beginDrag with hit (0,0,0) implicitly (scene pointer 50 -> intersector returns ~1?).
  // For this simplified test, just assert it changed (direction may depend on grab offset math)
  manip.updateDrag(()=>({ x:-0.5, y:0, z:0 }));
  expect(st.positions[0].x).not.toBeCloseTo(0,5);
    expect(camera.alpha).toBe(1);
    expect(camera.beta).toBe(2);
    expect(camera.radius).toBe(10);
    expect(camera.target.x).toBe(0);
  });
});
