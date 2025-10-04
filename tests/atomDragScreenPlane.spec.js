import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

// Babylon test stubs
if (!global.BABYLON) global.BABYLON = {};
if (!BABYLON.PointerEventTypes) BABYLON.PointerEventTypes = { POINTERDOWN:1 };
if (!BABYLON.Matrix) BABYLON.Matrix = { Identity: () => ({}) };

// Simple perspective-ish ray generator for tests
function makeRay(camPos, pointerX, pointerY, viewport){
  // Map pointer into normalized screen coords centered at 0 (orthographic approximation sufficient for plane tests)
  const nx = (pointerX - viewport.cx) / viewport.cx;
  const ny = (pointerY - viewport.cy) / viewport.cy;
  // Ray Dir: (nx, -ny, 1) from camera position (camera looks toward +Z)
  const dir = { x: nx, y: -ny, z: 1 };
  const len = Math.hypot(dir.x, dir.y, dir.z);
  return { origin: { ...camPos }, direction: { x:dir.x/len, y:dir.y/len, z:dir.z/len } };
}

describe('camera-aligned drag plane', () => {
  function setup(){
    const st = createMoleculeState({ elements:['C'], positions:[{x:0,y:0,z:0}], bonds:[] });
    const selection = { get:()=>st.selection, clickAtom:(i)=>{ st.selection={ kind:'atom', data:{ index:i } }; st.markSelectionChanged(); }, clickBond:()=>{}, clear:()=>{ st.selection={ kind:null, data:null }; st.markSelectionChanged(); } };
    const manip = createManipulationService(st, {});
    const cameraState = { position:{ x:0,y:0,z:-10 }, alpha:0, beta:0, radius:10, target:{x:0,y:0,z:0} };
    const camera = {
      get position(){ return cameraState.position; },
      detachControl(){}, attachControl(){},
      get alpha(){ return cameraState.alpha; }, get beta(){ return cameraState.beta; }, get radius(){ return cameraState.radius; }, get target(){ return cameraState.target; }
    };
    const viewport = { cx:100, cy:100 };
    const scene = {
      pointerX:100, pointerY:100,
      pick(){ return { hit:true, pickedMesh:{}, thinInstanceIndex:0 }; },
      onPointerObservable:{ add(fn){ this._cb=fn; } },
      getEngine(){ return { getRenderingCanvas(){ return { addEventListener(){}, removeEventListener(){} }; } }; },
      createPickingRay(x,y){ const r = makeRay(camera.position, x, y, viewport); return { origin: r.origin, direction: r.direction }; }
    };
    const view = { resolveAtomPick:()=>({ kind:'atom', index:0 }), resolveBondPick:()=>null };
    // Manual plane computation identical to production logic
    const atomPos = st.positions[0];
    const camPos = camera.position;
    const nx = atomPos.x - camPos.x, ny = atomPos.y - camPos.y, nz = atomPos.z - camPos.z;
    const len = Math.hypot(nx,ny,nz)||1; const planeNormal = { x:nx/len, y:ny/len, z:nz/len }; const planePoint = { ...atomPos };
    // Intersector using scene ray creator
    const intersector = (pp, pn)=>{
      const ray = scene.createPickingRay(scene.pointerX, scene.pointerY);
      const denom = pn.x*ray.direction.x + pn.y*ray.direction.y + pn.z*ray.direction.z;
      if (Math.abs(denom) < 1e-8) return null;
      const vx = pp.x - ray.origin.x, vy = pp.y - ray.origin.y, vz = pp.z - ray.origin.z;
      const t = (pn.x*vx + pn.y*vy + pn.z*vz)/denom;
      if (!isFinite(t) || t < 0) return null;
      return { x: ray.origin.x + ray.direction.x * t, y: ray.origin.y + ray.direction.y * t, z: ray.origin.z + ray.direction.z * t };
    };
    // Simulate selection then begin drag
    selection.clickAtom(0);
    manip.beginDrag(intersector, { planePoint, planeNormal });
    return { st, manip, scene, camera };
  }

  test('plane normal aligns with camera->atom and atom stays under pointer for X move', () => {
    const { st, manip, scene } = setup();
    // TEMP debug
    // eslint-disable-next-line no-console
    console.log('manip keys', Object.keys(manip));
    const dragState = manip._debug.getDragState();
    expect(dragState).toBeTruthy();
    // Camera->atom vector (0,0,0) - (0,0,-10) = (0,0,10) normalized -> (0,0,1)
    expect(dragState.planeNormal.z).toBeCloseTo(1, 5);
    expect(Math.abs(dragState.planeNormal.x)).toBeLessThan(1e-6);
    expect(Math.abs(dragState.planeNormal.y)).toBeLessThan(1e-6);
    // Initial grab offset ~ 0
    expect(Math.hypot(dragState.grabOffset.x, dragState.grabOffset.y, dragState.grabOffset.z)).toBeLessThan(1e-6);

    // Move pointer right
    scene.pointerX = 140; // +40 px
    manip.updateDrag((pp, pn) => {
      // Recompute intersection manually (same math as pickingService intersector simplified for test): intersect ray with plane
      const ray = { origin:{ x:0,y:0,z:-10 }, direction:{ x:(scene.pointerX-100)/100, y:-(scene.pointerY-100)/100, z:1 } };
      const dlen = Math.hypot(ray.direction.x, ray.direction.y, ray.direction.z); ray.direction.x/=dlen; ray.direction.y/=dlen; ray.direction.z/=dlen;
      const denom = pn.x*ray.direction.x + pn.y*ray.direction.y + pn.z*ray.direction.z;
      const vx = pp.x - ray.origin.x, vy = pp.y - ray.origin.y, vz = pp.z - ray.origin.z;
      const t = (pn.x*vx + pn.y*vy + pn.z*vz)/denom;
      return { x: ray.origin.x + ray.direction.x * t, y: ray.origin.y + ray.direction.y * t, z: ray.origin.z + ray.direction.z * t };
    });
    // Atom should have moved in +X; Z ~ 0 (remains on plane z=0)
    expect(st.positions[0].x).toBeGreaterThan(0);
    expect(st.positions[0].z).toBeCloseTo(0, 5);
  });

  test('pointer Y move changes atom Y with plane constraint', () => {
    const { st, manip, scene } = setup();
    scene.pointerY = 60; // move up
    manip.updateDrag((pp, pn) => {
      const ray = { origin:{ x:0,y:0,z:-10 }, direction:{ x:(scene.pointerX-100)/100, y:-(scene.pointerY-100)/100, z:1 } };
      const dlen = Math.hypot(ray.direction.x, ray.direction.y, ray.direction.z); ray.direction.x/=dlen; ray.direction.y/=dlen; ray.direction.z/=dlen;
      const denom = pn.x*ray.direction.x + pn.y*ray.direction.y + pn.z*ray.direction.z;
      const vx = pp.x - ray.origin.x, vy = pp.y - ray.origin.y, vz = pp.z - ray.origin.z;
      const t = (pn.x*vx + pn.y*vy + pn.z*vz)/denom;
      return { x: ray.origin.x + ray.direction.x * t, y: ray.origin.y + ray.direction.y * t, z: ray.origin.z + ray.direction.z * t };
    });
    expect(st.positions[0].y).toBeGreaterThan(0);
    expect(st.positions[0].z).toBeCloseTo(0,5);
  });
});
