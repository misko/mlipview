import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Minimal BABYLON stubs required
global.BABYLON = global.BABYLON || {};
if (!BABYLON.StandardMaterial) BABYLON.StandardMaterial = function(){ this.diffuseColor={}; this.emissiveColor={}; this.alpha=1; };
if (!BABYLON.Color3) BABYLON.Color3 = function(r,g,b){ this.r=r; this.g=g; this.b=b; this.clone=()=>new BABYLON.Color3(r,g,b); this.scale=f=>new BABYLON.Color3(r*f,g*f,b*f); };
if (!BABYLON.MeshBuilder) BABYLON.MeshBuilder = { CreateSphere:()=>({ thinInstanceSetBuffer(){}, }), CreateCylinder:()=>({ thinInstanceSetBuffer(){} }) };
if (!BABYLON.Matrix) BABYLON.Matrix = class { constructor(){ this.m=new Float32Array(16);} static Compose(){ return new BABYLON.Matrix(); } };
if (!BABYLON.Vector3) BABYLON.Vector3 = class { constructor(x,y,z){ this.x=x; this.y=y; this.z=z;} length(){return Math.hypot(this.x,this.y,this.z);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new BABYLON.Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} normalizeToNew(){ const L=this.length()||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L);} };
if (!BABYLON.Quaternion) BABYLON.Quaternion = { Identity:()=>({}), RotationAxis:()=>({}) };

function makeSceneStub(){
  return { onPointerObservable:{ add(){} }, }; // highlight tests don't need picking
}

describe('selection highlight', () => {
  test('atom highlight appears and bond highlight hidden then switched', () => {
    const state = createMoleculeState({ elements:['C','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0}], bonds:[{i:0,j:1}] });
    const scene = makeSceneStub();
    const view = createMoleculeView(scene, state);
    const sel = createSelectionService(state);
  sel.clickAtom(0); // triggers selectionChanged event -> highlight update
    const { highlight } = view._internals;
    expect(highlight.atom.isVisible).toBe(true);
    expect(highlight.bond.isVisible).toBe(false);
  // Check material alpha roughly in expected range (0.3 - 0.5)
  expect(highlight.atom.material.alpha).toBeGreaterThan(0.3);
  expect(highlight.atom.material.alpha).toBeLessThan(0.6);
    sel.clickBond({ i:0,j:1,key:'C-H', index:0 });
    expect(highlight.atom.isVisible).toBe(false);
    expect(highlight.bond.isVisible).toBe(true);
    // Bond highlight alpha also within expected visual range
    expect(highlight.bond.material.alpha).toBeGreaterThan(0.3);
    expect(highlight.bond.material.alpha).toBeLessThan(0.6);
    sel.clear();
    expect(highlight.atom.isVisible).toBe(false);
    expect(highlight.bond.isVisible).toBe(false);
  });
});
