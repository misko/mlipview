// Force rendering test: verifies that force vectors are instantiated as thin instances
// with fixed-length red cylinders when forces are present and 'forcesChanged' is emitted.

// Minimal Babylon stubs sufficient for moleculeView logic
if (!global.BABYLON) global.BABYLON = {};
const BABYLON = global.BABYLON;

class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } length(){ return Math.hypot(this.x,this.y,this.z); } addInPlace(v){ this.x+=v.x; this.y+=v.y; this.z+=v.z; return this; } clone(){ return new Vector3(this.x,this.y,this.z); } static Dot(a,b){ return a.x*b.x + a.y*b.y + a.z*b.z; } static Cross(a,b){ return new Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); } static Up(){ return new Vector3(0,1,0); } static Right(){ return new Vector3(1,0,0); } normalize(){ const L=this.length()||1; this.x/=L; this.y/=L; this.z/=L; return this; } normalizeToNew(){ return this.clone().normalize(); } }
class Quaternion { constructor(x=0,y=0,z=0,w=1){ this.x=x; this.y=y; this.z=z; this.w=w; } static Identity(){ return new Quaternion(0,0,0,1); } static RotationAxis(axis, angle){ const s=Math.sin(angle/2); return new Quaternion(axis.x*s,axis.y*s,axis.z*s,Math.cos(angle/2)); } }
class Matrix { constructor(){ this.m=new Float32Array(16); } static Compose(scale, rot, pos){ // Very simplified compose: scale -> rotation (ignored here) -> translation
  // We only need distinct matrix objects; actual numeric correctness not required for test assertions currently.
  const mat = new Matrix(); mat.m[12]=pos.x; mat.m[13]=pos.y; mat.m[14]=pos.z; return mat; } }
class Color3 { constructor(r=0,g=0,b=0){ this.r=r; this.g=g; this.b=b; } clone(){ return new Color3(this.r,this.g,this.b); } scale(f){ return new Color3(this.r*f,this.g*f,this.b*f); } }
class StandardMaterial { constructor(name){ this.name=name; this.diffuseColor=new Color3(1,1,1); this.emissiveColor=new Color3(0,0,0); this.specularColor=new Color3(0,0,0); this.alpha=1; } }
class Mesh { constructor(name){ this.name=name; this.isPickable=true; this.thinInstanceEnablePicking=false; this.material=null; this._buffers={}; this.isVisible=true; }
  thinInstanceSetBuffer(kind, arr){ this._buffers[kind]=arr; }
  setEnabled(on){ this.isVisible=on; }
}
const MeshBuilder = { CreateSphere:(n)=>new Mesh(n), CreateCylinder:(n)=>new Mesh(n), CreateLines:(n)=>new Mesh(n) };
BABYLON.Vector3=Vector3; BABYLON.Quaternion=Quaternion; BABYLON.Matrix=Matrix; BABYLON.Color3=Color3; BABYLON.StandardMaterial=StandardMaterial; BABYLON.MeshBuilder=MeshBuilder; BABYLON.Material={ MATERIAL_ALPHABLEND:2 };

// Simple event bus mock
function createBus(){ const ls={}; return { on(ev,fn){ (ls[ev]=ls[ev]||[]).push(fn); }, emit(ev){ (ls[ev]||[]).forEach(f=>f()); } }; }

import { createMoleculeView } from '../public/render/moleculeView.js';

describe('force rendering', () => {
  test('creates force thin instances after forcesChanged', () => {
    const scene = { meshes:[], onBeforeRenderObservable:{ add(){} } };
    const molState = {
      elements:['C','H','H'],
      positions:[ {x:0,y:0,z:0},{x:1,y:0,z:0},{x:-1,y:0,z:0} ],
      bonds:[],
      bus: createBus(),
      selection:null
    };
    // Seed forces (arbitrary directions)
    molState.forces = [ [0,1,0], [1,1,0], [ -1, 0.5, 0 ] ];

    const view = createMoleculeView(scene, molState);
    // Emit forcesChanged to trigger force rebuild
    molState.bus.emit('forcesChanged');

    const fg = view._internals.forceGroups.get('force');
    expect(fg).toBeTruthy();
    // matrix buffer should exist and have entries (each matrix is 16 floats)
    const buf = fg.master._buffers['matrix'];
    expect(buf).toBeTruthy();
    expect(buf.length).toBeGreaterThan(0);
    // Color buffer should be solid red
    const colors = fg.master._buffers['color'];
    expect(colors).toBeTruthy();
    // First force instance RGBA approx (0.95,0.05,0.05,1)
    expect(colors[0]).toBeCloseTo(0.95, 2);
    expect(colors[1]).toBeCloseTo(0.05, 2);
    expect(colors[2]).toBeCloseTo(0.05, 2);
    expect(colors[3]).toBeCloseTo(1.0, 2);
  });
});
