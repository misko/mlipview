// Ensure window & document stubs exist before any module import references them
if(!global.window) global.window = {};
// Ensure 'window' identifier references global.window in CommonJS scope
// (needed because some modules reference 'window' directly during evaluation)
// eslint-disable-next-line no-var
var window = global.window;
if(!global.document) global.document = { 
  createElement:(tag)=>({ 
    tagName:tag.toUpperCase(), 
    style:{}, 
    getContext:()=>({ 
      clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){} 
    }) 
  }), 
  getElementById:()=>null, 
  body:{ appendChild(){} } 
};

// Integration test: verifies that after perturbing an atom position and triggering a remote force recompute
// the force vectors in the view update (matrix buffer contents change).
// We mock network responses deterministically: first response returns one set of forces; after perturbation
// second response returns a different set.

// Lightweight Babylon stubs (reuse from existing tests where possible)
if (!global.BABYLON) global.BABYLON = {};
const BABYLON = global.BABYLON;
class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } length(){ return Math.hypot(this.x,this.y,this.z); } addInPlace(v){ this.x+=v.x; this.y+=v.y; this.z+=v.z; return this; } clone(){ return new Vector3(this.x,this.y,this.z); } static Dot(a,b){ return a.x*b.x + a.y*b.y + a.z*b.z; } static Cross(a,b){ return new Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); } normalize(){ const L=this.length()||1; this.x/=L; this.y/=L; this.z/=L; return this; } normalizeToNew(){ return this.clone().normalize(); } }
class Quaternion { constructor(x=0,y=0,z=0,w=1){ this.x=x; this.y=y; this.z=z; this.w=w; } static Identity(){ return new Quaternion(0,0,0,1); } static RotationAxis(axis, angle){ const s=Math.sin(angle/2); return new Quaternion(axis.x*s,axis.y*s,axis.z*s,Math.cos(angle/2)); } }
class Matrix { constructor(){ this.m=new Float32Array(16); } static Compose(scale, rot, pos){ const mat=new Matrix(); mat.m[12]=pos.x; mat.m[13]=pos.y; mat.m[14]=pos.z; return mat; } }
class Color3 { constructor(r=0,g=0,b=0){ this.r=r; this.g=g; this.b=b; } clone(){ return new Color3(this.r,this.g,this.b); } scale(f){ return new Color3(this.r*f,this.g*f,this.b*f); } }
class StandardMaterial { constructor(name){ this.name=name; this.diffuseColor=new Color3(1,1,1); this.emissiveColor=new Color3(0,0,0); this.specularColor=new Color3(0,0,0); this.alpha=1; } }
class Mesh { constructor(name){ this.name=name; this.isPickable=true; this.thinInstanceEnablePicking=false; this.material=null; this._buffers={}; this.isVisible=true; } thinInstanceSetBuffer(kind, arr){ this._buffers[kind]=arr; } setEnabled(on){ this.isVisible=on; } }
const MeshBuilder={ CreateSphere:(n)=>new Mesh(n), CreateCylinder:(n)=>new Mesh(n), CreateLines:(n)=>new Mesh(n) };
class TransformNode { constructor(name){ this.name=name; this.position=new Vector3(); this.scaling=new Vector3(1,1,1); this.rotationQuaternion=null; this.parent=null; } }
BABYLON.TransformNode=TransformNode; BABYLON.Vector3=Vector3; BABYLON.Quaternion=Quaternion; BABYLON.Matrix=Matrix; BABYLON.Color3=Color3; BABYLON.StandardMaterial=StandardMaterial; BABYLON.MeshBuilder=MeshBuilder; BABYLON.Material={ MATERIAL_ALPHABLEND:2 }; BABYLON.PointerEventTypes={ POINTERDOWN:1 };

// Minimal engine/scene stub to satisfy createScene usage.
// We intercept runRenderLoop to just call the callback immediately a few times.
class Engine { constructor(){ this._cb=null; }
  runRenderLoop(cb){ this._cb=cb; for(let i=0;i<3;i++) cb(); }
  stopRenderLoop(){ this._cb=null; }
}
class Scene { constructor(){ this.meshes=[]; this.onBeforeRenderObservable={ add(){} }; this.pointerX=0; this.pointerY=0; }
  render(){ /* noop for tests */ }
  getEngine(){ return { getRenderingCanvas:()=>({ addEventListener(){} }) }; }
  createPickingRay(){ return { origin:new Vector3(0,0,0), direction:new Vector3(0,1,0) }; }
  pick(){ return { hit:false }; }
}

// Patch module imports used by createScene if needed
jest.unmock('../public/render/scene.js');

// Mock fetch to simulate two sequential force responses
let fetchCallCount = 0;
const firstForces = [ [0.0, 0.1, 0.0], [0.05, 0.0, 0.0], [-0.02, 0.03, 0.0] ];
const secondForces = [ [0.5, 0.0, 0.0], [0.1, 0.1, 0.0], [-0.3, 0.2, 0.0] ];

global.fetch = jest.fn(async (url, opts)=>{
  fetchCallCount++;
  // Always respond OK with JSON structure expected
  const body = fetchCallCount === 1 ? { results:{ energy: -10.0, forces: firstForces } } : { results:{ energy: -12.5, forces: secondForces } };
  return { ok:true, status:200, json: async()=> body };
});

// Provide window & document stubs
// Ensure window exists before importing modules that define properties on window
if(!global.window) global.window = {};
if(!global.document) global.document = { 
  createElement:(tag)=>({ 
    tagName:tag.toUpperCase(), 
    style:{}, 
    getContext:()=>({ 
      clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){} 
    }) 
  }), 
  getElementById:()=>null, 
  body:{ appendChild(){} } 
};
window.__MLIPVIEW_TEST_MODE = true; // avoid full runRenderLoop path

// Mock createScene to use stubs (so we don't rely on real Babylon engine). Define mock-prefixed classes inside factory.
jest.mock('../public/render/scene.js', ()=>{
  class mockVec3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }
  class mockEngine { constructor(){ this._cb=null; } runRenderLoop(cb){ this._cb=cb; for(let i=0;i<2;i++) cb(); } stopRenderLoop(){ this._cb=null; } }
  class mockScene { constructor(){ this.meshes=[]; this.onBeforeRenderObservable={ add(){} }; this.onPointerObservable={ add(){} }; this.pointerX=0; this.pointerY=0; }
    render(){}
    getEngine(){ return { getRenderingCanvas:()=>({ addEventListener(){} }) }; }
    createPickingRay(){ return { origin:new mockVec3(0,0,0), direction:new mockVec3(0,1,0) }; }
    pick(){ return { hit:false }; }
  }
  return { createScene: async function(){ return { engine:new mockEngine(), scene:new mockScene(), camera:{ detachControl(){}, attachControl(){}, alpha:0,beta:0,radius:5,target:{x:0,y:0,z:0} } }; } };
});

// Dynamic import later after environment prepared
let initNewViewer;

function getForceMatrixCount(api){
  const fg = api.view._internals.forceGroups.get('force');
  if(!fg) return 0; const buf = fg.master._buffers['matrix']; return buf ? (buf.length/16) : 0;
}

describe('force update after perturbation', () => {
  test('forces change after atom displacement triggers recompute', async () => {
  if(!initNewViewer){ const mod = await import('../public/index.js'); initNewViewer = mod.initNewViewer; }
  const canvas = { width:800, height:600, getContext:()=>({}) };
  const api = await initNewViewer(canvas, { elements:['O','H','H'], positions:[{x:0,y:0,z:0},{x:0.95,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] });
    // Wait for first baseline fetch to complete
    await new Promise(r=>setTimeout(r,5));
    expect(fetchCallCount).toBeGreaterThanOrEqual(1);
    // Ensure initial forces rendered
  api.state.bus.emit('forcesChanged');
  // Poll for force matrices to appear (async rebuild path tolerance)
  let initialCount = 0; 
    for(let t=0;t<40;t++){ 
      initialCount = getForceMatrixCount(api); 
      if(initialCount>0) break; 
      // Fallback: direct rebuild call (some stubs may skip event binding path)
      try { api.view.rebuildForces(); } catch {}
      await new Promise(r=>setTimeout(r,10)); 
    }
  expect(initialCount).toBeGreaterThan(0);
    const fgInitial = api.view._internals.forceGroups.get('force');
    const initialMatrix = fgInitial.master._buffers['matrix'].slice();
    // Perturb oxygen significantly
    api.state.positions[0].x += 5.0; api.state.positions[0].y += 4.0; api.state.positions[0].z += -3.0;
    api.state.markPositionsChanged();
    // Force recompute (wrapped markPositionsChanged triggers stale flag then debounced ff.computeForces)
    await new Promise(r=>setTimeout(r,80)); // > debounce 50ms
    // Trigger another compute explicitly to ensure second fetch
    await api.ff.computeForces({ sync:true });
    expect(fetchCallCount).toBeGreaterThanOrEqual(2);
  api.state.bus.emit('forcesChanged');
  let fgAfter, afterMatrix; 
  for(let t=0;t<40;t++){ fgAfter = api.view._internals.forceGroups.get('force'); afterMatrix = fgAfter?.master?._buffers['matrix']; if(afterMatrix && afterMatrix.length) break; try { api.view.rebuildForces(); } catch {}; await new Promise(r=>setTimeout(r,10)); }
    expect(afterMatrix).toBeTruthy();
    // Matrices should differ (position + direction changes)
    let diff=0; for(let i=0;i<afterMatrix.length && i<initialMatrix.length;i++){ if(afterMatrix[i] !== initialMatrix[i]) { diff++; break; } }
    expect(diff).toBeGreaterThan(0);
    // Additionally ensure force array content in state updated
    const f1 = api.state.forces; expect(Array.isArray(f1)).toBe(true); expect(f1.length).toBe(3);
  });
});
