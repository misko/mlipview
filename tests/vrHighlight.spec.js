// Basic non-VR highlight sanity to ensure file not empty and highlight sphere creation works.
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

global.BABYLON = global.BABYLON || {};
if (!BABYLON.StandardMaterial) BABYLON.StandardMaterial = function(){};
if (!BABYLON.Color3) BABYLON.Color3 = function(r,g,b){ this.r=r; this.g=g; this.b=b; this.clone=()=>new BABYLON.Color3(r,g,b); this.scale=f=>new BABYLON.Color3(r*f,g*f,b*f); };
if (!BABYLON.MeshBuilder) BABYLON.MeshBuilder = { CreateSphere:()=>({ thinInstanceSetBuffer(){}, }), CreateCylinder:()=>({ thinInstanceSetBuffer(){} }) };
if (!BABYLON.Matrix) BABYLON.Matrix = class { constructor(){ this.m=new Float32Array(16);} static Compose(){ return new BABYLON.Matrix(); } };
if (!BABYLON.Vector3) BABYLON.Vector3 = class { constructor(x,y,z){ this.x=x; this.y=y; this.z=z;} length(){return Math.hypot(this.x,this.y,this.z);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new BABYLON.Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} normalizeToNew(){ const L=this.length()||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L);} };
if (!BABYLON.Quaternion) BABYLON.Quaternion = { Identity:()=>({}), RotationAxis:()=>({}) };

describe('vrHighlight placeholder', () => {
	test('atom selection shows highlight sphere', () => {
		const state = createMoleculeState({ elements:['H','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0}], bonds:[] });
		const scene = { onPointerObservable:{ add(){} } };
		const view = createMoleculeView(scene, state);
		const sel = createSelectionService(state);
		sel.clickAtom(1);
		expect(view._internals.highlight.atom.isVisible).toBe(true);
	});
});
