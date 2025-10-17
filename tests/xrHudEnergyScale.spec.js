/**
 * xrHudEnergyScale.spec
 * Ensures topEnergyScaleMult multiplies plane width/height proportionally.
 */
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';
if(!global.window) global.window = { location:{ search:'' } };
if(!global.document) global.document = { createElement: (tag)=> tag==='canvas'?{ width:0,height:0,getContext:()=>({}), toDataURL:()=>'' } : {} };
if(!global.BABYLON) global.BABYLON = {};
const B=global.BABYLON; function v(x=0,y=0,z=0){ return { x,y,z, add(o){return v(x+o.x,y+o.y,z+o.z);}, subtract(o){return v(x-o.x,y-o.y,z-o.z);}, scale(s){return v(x*s,y*s,z*s);}, normalize(){ const L=Math.hypot(x,y,z)||1; return v(x/L,y/L,z/L); }, copyFrom(o){ this.x=o.x; this.y=o.y; this.z=o.z; } }; }
if(!B.Axis) B.Axis={ Z:v(0,0,1), Y:v(0,1,0) };
if(!B.MeshBuilder) B.MeshBuilder={}; if(!B.MeshBuilder.CreatePlane) B.MeshBuilder.CreatePlane=(n,o,s)=>({ name:n, opts:o, scene:s, position:v(), rotation:{y:0}, isPickable:true });
if(!B.GUI) B.GUI={}; if(!B.GUI.Control) B.GUI.Control={ HORIZONTAL_ALIGNMENT_CENTER:0 };
if(!B.GUI.AdvancedDynamicTexture) B.GUI.AdvancedDynamicTexture=class { constructor(){ this._rootContainer={ children:[] }; } static CreateForMesh(){ return new B.GUI.AdvancedDynamicTexture(); } addControl(c){ this._rootContainer.children.push(c); } };
if(!B.GUI.StackPanel) B.GUI.StackPanel=class { constructor(){ this.children=[]; } addControl(c){ this.children.push(c); } };
if(!B.GUI.Rectangle) B.GUI.Rectangle=class { constructor(){ this.children=[]; } addControl(c){ this.children.push(c); } };
if(!B.GUI.Button) B.GUI.Button=class { static CreateSimpleButton(id,text){ return { id,text,background:'', onPointerDownObservable:{add:()=>{}}, onPointerUpObservable:{add:()=>{}}, onPointerEnterObservable:{add:()=>{}}, onPointerOutObservable:{add:()=>{}} }; } };
if(!B.GUI.Image){ B.GUI.Image=class { constructor(id,src){ this.id=id; this.source=src; this.width=''; this.height=''; this.stretch=0; } markAsDirty(){} }; B.GUI.Image.STRETCH_UNIFORM=0; }
if(!B.GUI.TextBlock) B.GUI.TextBlock=class { constructor(id,text){ this.id=id; this.text=text; } };

function makeScene(){ return { activeCamera:{ fov:Math.PI/2, position:v(), getDirection(ax){ if(ax===B.Axis.Z) return v(0,0,1); if(ax===B.Axis.Y) return v(0,1,0); return v(); } }, onBeforeRenderObservable:{ _l:[], add(fn){ this._l.push(fn); return fn; } } }; }

describe('XR HUD top energy scale multiplier', () => {
	test('scale multiplier doubles width/height roughly', () => {
		const s1 = makeScene(); if(window.__XR_HUD_FALLBACK){ delete window.__XR_HUD_FALLBACK; delete window.__XR_HUD_FALLBACK_TOP; }
		const r1 = ensureWorldHUD({ scene:s1, topEnergyScaleMult:1 }); s1.onBeforeRenderObservable._l.forEach(fn=>fn());
		if(window.__XR_HUD_FALLBACK){ delete window.__XR_HUD_FALLBACK; delete window.__XR_HUD_FALLBACK_TOP; }
		const s2 = makeScene(); const r2 = ensureWorldHUD({ scene:s2, topEnergyScaleMult:2 }); s2.onBeforeRenderObservable._l.forEach(fn=>fn());
		expect(r2.planeTop.opts.width / r1.planeTop.opts.width).toBeGreaterThan(1.9);
		expect(r2.planeTop.opts.width / r1.planeTop.opts.width).toBeLessThan(2.1);
		expect(r2.planeTop.opts.height / r1.planeTop.opts.height).toBeGreaterThan(1.9);
		expect(r2.planeTop.opts.height / r1.planeTop.opts.height).toBeLessThan(2.1);
	});
});
