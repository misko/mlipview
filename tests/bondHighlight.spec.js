import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Reuse highlight spec style BABYLON stubs but include minimal picking support for bond
// Only what we need for view building & selectionChanged events.
if (!global.BABYLON) global.BABYLON = {};
const B = global.BABYLON;
if (!B.StandardMaterial) B.StandardMaterial = function(){ this.diffuseColor={}; this.emissiveColor={}; this.alpha=1; };
if (!B.Color3) B.Color3 = function(r,g,b){ this.r=r; this.g=g; this.b=b; this.clone=()=>new B.Color3(r,g,b); this.scale=f=>new B.Color3(r*f,g*f,b*f); };
if (!B.MeshBuilder) B.MeshBuilder = { CreateSphere:(n,o,s)=>({ name:n, thinInstanceSetBuffer(){}, }), CreateCylinder:(n,o,s)=>({ name:n, thinInstanceSetBuffer(){}, }) };
if (!B.Matrix) B.Matrix = class { constructor(){ this.m=new Float32Array(16);} static Compose(){ return new B.Matrix(); } };
if (!B.Vector3) B.Vector3 = class { constructor(x,y,z){ this.x=x; this.y=y; this.z=z;} length(){return Math.hypot(this.x,this.y,this.z);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new B.Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} normalizeToNew(){ const L=this.length()||1; return new B.Vector3(this.x/L,this.y/L,this.z/L);} };
if (!B.Quaternion) B.Quaternion = { Identity:()=>({}), RotationAxis:()=>({}) };

function makeSceneStub(){
  return { onPointerObservable:{ add(){} } };
}

describe('bond highlight selection', () => {
  test('selecting a bond makes bond highlight visible and atom highlight hidden', () => {
    const state = createMoleculeState({
      elements:['C','H'],
      positions:[{x:0,y:0,z:0},{x:1,y:0,z:0}],
      bonds:[{ i:0, j:1 }]
    });
    const scene = makeSceneStub();
    const view = createMoleculeView(scene, state);
    const selection = createSelectionService(state);
    // First select atom to ensure atom highlight path works
    selection.clickAtom(0);
    expect(view._internals.highlight.atom.isVisible).toBe(true);
    expect(view._internals.highlight.bond.isVisible).toBe(false);
  // Now select bond
    selection.clickBond({ i:0, j:1, key:'C-H', index:0 });
    // EXPECTATION: atom highlight hidden, bond highlight visible
    expect(view._internals.highlight.atom.isVisible).toBe(false);
    expect(view._internals.highlight.bond.isVisible).toBe(true);
  });

  test('bond master is pickable enabling future pick-based selection', () => {
    const state = createMoleculeState({
      elements:['C','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0}], bonds:[{ i:0,j:1 }]
    });
    const scene = makeSceneStub();
    const view = createMoleculeView(scene, state);
    // Find bond group master
    const bondGroups = view._internals.bondGroups;
    const first = bondGroups.values().next().value;
    expect(first.master.isPickable).toBe(true);
    // Simulate a Babylon pick structure hitting thin instance 0
    const fakePick = { hit:true, pickedMesh:first.master, thinInstanceIndex:0 };
    const resolved = view.resolveBondPick(fakePick);
    expect(resolved).toBeTruthy();
    expect(resolved.kind).toBe('bond');
  });
});
