// Verifies force vectors update over multiple MD steps.
if(!global.window) global.window = {}; var window = global.window;
if(!global.document) global.document = { createElement:(t)=>({ tagName:t.toUpperCase(), style:{}, getContext:()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){} }) }), getElementById:()=>null, body:{ appendChild(){} } };

if (!global.BABYLON) global.BABYLON = {};
const BABYLON = global.BABYLON;
class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } length(){ return Math.hypot(this.x,this.y,this.z); } addInPlace(v){ this.x+=v.x; this.y+=v.y; this.z+=v.z; return this; } scale(f){ return new Vector3(this.x*f,this.y*f,this.z*f); } normalize(){ const L=this.length()||1; this.x/=L; this.y/=L; this.z/=L; return this; } normalizeToNew(){ return new Vector3(this.x,this.y,this.z).normalize(); } static Dot(a,b){ return a.x*b.x + a.y*b.y + a.z*b.z; } static Cross(a,b){ return new Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); } }
class Quaternion { constructor(x=0,y=0,z=0,w=1){ this.x=x; this.y=y; this.z=z; this.w=w; } static Identity(){ return new Quaternion(); } static RotationAxis(axis, angle){ const s=Math.sin(angle/2); return new Quaternion(axis.x*s,axis.y*s,axis.z*s,Math.cos(angle/2)); } }
class Matrix { constructor(){ this.m=new Float32Array(16);} static Compose(scale,_r,pos){ const m=new Matrix(); m.m[12]=pos.x; m.m[13]=pos.y; m.m[14]=pos.z; return m; } }
class Color3 { constructor(r=0,g=0,b=0){ this.r=r; this.g=g; this.b=b; } clone(){ return new Color3(this.r,this.g,this.b);} scale(f){ return new Color3(this.r*f,this.g*f,this.b*f);} }
class StandardMaterial { constructor(){ this.diffuseColor=new Color3(1,1,1); this.emissiveColor=new Color3(); this.specularColor=new Color3(); this.alpha=1; } }
class Mesh { constructor(){ this._buffers={}; this.isVisible=true; this.material=null; this.thinInstanceEnablePicking=false; this.parent=null; } thinInstanceSetBuffer(k,a){ this._buffers[k]=a; } setEnabled(on){ this.isVisible=on; } dispose(){} }
const MeshBuilder={ CreateSphere:()=>new Mesh(), CreateCylinder:()=>new Mesh(), CreateLines:()=>new Mesh() };
class TransformNode { constructor(){ this.position=new Vector3(); this.scaling=new Vector3(1,1,1); this.parent=null; } }
BABYLON.Vector3=Vector3; BABYLON.Quaternion=Quaternion; BABYLON.Matrix=Matrix; BABYLON.Color3=Color3; BABYLON.StandardMaterial=StandardMaterial; BABYLON.MeshBuilder=MeshBuilder; BABYLON.TransformNode=TransformNode; BABYLON.Material={ MATERIAL_ALPHABLEND:2 }; BABYLON.PointerEventTypes={ POINTERDOWN:1 };

window.__MLIPVIEW_TEST_MODE = true;
// Lower minimum force arrow length so incremental force magnitude changes alter matrices
window.FORCE_MIN = 0.01;

// Mock scene creation
jest.mock('../public/render/scene.js', ()=>{
  class mockVec3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }
  class mockEngine { runRenderLoop(cb){ for(let i=0;i<2;i++) cb(); } stopRenderLoop(){} }
  class mockScene { constructor(){ this.onBeforeRenderObservable={ add(){} }; this.onPointerObservable={ add(){} }; this.pointerX=0; this.pointerY=0; } render(){} getEngine(){ return { getRenderingCanvas:()=>({ addEventListener(){} }) }; } createPickingRay(){ return { origin:new mockVec3(0,0,0), direction:new mockVec3(0,1,0) }; } pick(){ return { hit:false }; } }
  return { createScene: async ()=> ({ engine:new mockEngine(), scene:new mockScene(), camera:{ detachControl(){}, attachControl(){}, alpha:0,beta:0,radius:5,target:{x:0,y:0,z:0} } }) };
});

// WebSocket-only: inject baseline idle forces and then changing forces for each mdStep()
import { getWS } from '../public/fairchem_ws_client.js';
// Deterministic force pattern
const baseForces=[ [0.05,0.0,0.0],[0.0,0.06,0.0],[0.0,0.0,0.07] ];
function mdForces(step){ return baseForces.map(f=>[ f[0]*(1+0.25*step), f[1]*(1+0.18*step), f[2]*(1+0.12*step) ]); }

let initNewViewer;
function forceInstanceCount(api){ const fg=api.view._internals.forceGroups.get('force'); if(!fg) return 0; const buf=fg.master._buffers['matrix']; return buf? buf.length/16 : 0; }

describe('md forces visualization', () => {
  test('forces update over 10 md steps', async () => {
    // Stub WebSocket and inject baseline + md frames
    const origWS = global.WebSocket;
    class FakeWS { constructor(){ this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0);} send(){} close(){} onopen(){} onmessage(){} onerror(){} }
    global.WebSocket = FakeWS;
    const ws = getWS();
    let mdStepCount = 0;
    ws.setTestHook((msg)=>{
      try {
        if (msg && msg.simulationParams) {
          mdStepCount += 1;
          setTimeout(()=>{
            ws.injectTestResult({ positions:[[0,0,0],[0.95,0,0],[-0.24,0.93,0]], forces: mdForces(mdStepCount), energy: -3.0 - 0.05*mdStepCount });
          }, 0);
        } else if (msg && msg.type) {
          // Idle compute baseline after init
          setTimeout(()=>{ ws.injectTestResult({ energy:-3.0, forces: mdForces(0) }); }, 0);
        }
      } catch {}
    });

    const mod = await import('../public/index.js'); initNewViewer = mod.initNewViewer;
    const api = await initNewViewer({ width:800,height:600,getContext:()=>({}) }, { elements:['O','H','H'], positions:[{x:0,y:0,z:0},{x:0.95,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] });
    await new Promise(r=>setTimeout(r,20));
    api.state.bus.emit('forcesChanged');
    for(let t=0;t<20 && forceInstanceCount(api)===0;t++){ try{ api.view.rebuildForces(); }catch{} await new Promise(r=>setTimeout(r,10)); }
    const initialMatrix = api.view._internals.forceGroups.get('force').master._buffers['matrix'].slice();
  for(let i=0;i<10;i++){ await api.mdStep(); api.state.bus.emit('forcesChanged'); try{ api.view.rebuildForces(); }catch{} await new Promise(r=>setTimeout(r,10)); }
  let afterMatrix; for(let p=0;p<20;p++){ afterMatrix = api.view._internals.forceGroups.get('force').master._buffers['matrix']; if(afterMatrix && afterMatrix.length) break; await new Promise(r=>setTimeout(r,10)); }
    let changed=false; for(let i=0;i<afterMatrix.length && i<initialMatrix.length;i++){ if(afterMatrix[i]!==initialMatrix[i]) { changed=true; break; } }
    expect(forceInstanceCount(api)).toBeGreaterThan(0);
    expect(changed).toBe(true);
    ws.setTestHook(null); global.WebSocket = origWS;
  });
});
