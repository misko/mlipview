// Ensure fetch exists (Node 20 should, but polyfill defensively for Jest env inconsistencies)
try {
  if (typeof fetch === 'undefined') {
    // cross-fetch/polyfill will install globalThis.fetch
    require('cross-fetch/polyfill');
  }
} catch (e) {
  // eslint-disable-next-line no-console
  console.warn('[jest.setup] fetch polyfill not installed', e?.message||e);
}

// Global BABYLON mock for tests needing molecule/cell logic
if (!global.BABYLON) {
  const BABYLON = {
  Engine: class Engine { constructor(canvas){ this._canvas = canvas||{ getBoundingClientRect(){ return { left:0, top:0, width:800, height:600 }; } }; this._loop=null; } runRenderLoop(fn){ this._loop = setInterval(()=>{ try{ fn(); }catch{} }, 16); } stopRenderLoop(){ if(this._loop){ clearInterval(this._loop); this._loop=null; } } getRenderingCanvas(){ return this._canvas; } },
  Scene: class Scene { constructor(){ this.onPointerObservable={ _l:[], add(fn){ this._l.push(fn); }, notify(ev){ this._l.forEach(f=>f(ev)); }, notifyObservers(ev){ this._l.forEach(f=>f(ev)); } }; this.onBeforeRenderObservable={ _l:[], add(fn){ this._l.push(fn); return fn; }, remove(fn){ const i=this._l.indexOf(fn); if(i>=0) this._l.splice(i,1); } }; this.activeCamera=null; this.clearColor=null; } getLightByName(){ return null; } pick(){ return { hit:false }; } createPickingRay(){ return { origin:new BABYLON.Vector3(0,0,0), direction:new BABYLON.Vector3(0,0,1) }; } render(){ /* noop */ } },
  Vector3: class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } add(v){ return new BABYLON.Vector3(this.x+v.x,this.y+v.y,this.z+v.z);} addInPlace(v){ this.x+=v.x; this.y+=v.y; this.z+=v.z; return this;} scale(s){ return new BABYLON.Vector3(this.x*s,this.y*s,this.z*s);} length(){ return Math.hypot(this.x,this.y,this.z);} negate(){ return new BABYLON.Vector3(-this.x,-this.y,-this.z);} clone(){ return new BABYLON.Vector3(this.x,this.y,this.z);} subtract(v){ return new BABYLON.Vector3(this.x-v.x,this.y-v.y,this.z-v.z);} normalize(){ const L=this.length()||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L);} normalizeToNew(){ return this.normalize(); } toString(){ return `{X: ${this.x} Y: ${this.y} Z: ${this.z}}`; } static Zero(){ return new BABYLON.Vector3(0,0,0);} static Up(){ return new BABYLON.Vector3(0,1,0);} static Right(){ return new BABYLON.Vector3(1,0,0);} static Dot(a,b){ return a.x*b.x + a.y*b.y + a.z*b.z; } static Cross(a,b){ return new BABYLON.Vector3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);} },
  Quaternion: class Quaternion { static Identity(){ return new BABYLON.Quaternion(); } static RotationAxis(){ return new BABYLON.Quaternion(); } },
    Matrix: class Matrix { constructor(){ this.m=new Float32Array(16);} static Compose(scale,rot,pos){ const m=new BABYLON.Matrix(); m.m[12]=pos.x; m.m[13]=pos.y; m.m[14]=pos.z; return m;} clone(){ const c=new BABYLON.Matrix(); c.m.set(this.m); return c;} copyToArray(arr,off){ for(let i=0;i<16;i++) arr[off+i]=this.m[i]||0; } static Identity(){ return new BABYLON.Matrix(); } },
    Color3: class Color3 { constructor(r,g,b){ this.r=r; this.g=g; this.b=b; } clone(){ return new BABYLON.Color3(this.r,this.g,this.b);} scale(f){ return new BABYLON.Color3(this.r*f,this.g*f,this.b*f);} },
    Color4: class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a;} },
    StandardMaterial: class StandardMaterial { constructor(){ this.alpha=1; this.diffuseColor=new BABYLON.Color3(1,1,1); this.emissiveColor=new BABYLON.Color3(0,0,0);} },
    MeshBuilder: { CreateSphere(name,opts,scene){ return new BABYLON.Mesh(name); }, CreateCylinder(name,opts,scene){ return new BABYLON.Mesh(name);} },
    Mesh: class Mesh { constructor(name){ this.name=name; this.material=null; this.isPickable=false; this.thinInstanceEnablePicking=false; this.alwaysSelectAsActiveMesh=false; this.isVisible=true; this.scaling=new BABYLON.Vector3(1,1,1); this.position=new BABYLON.Vector3(); } dispose(){} thinInstanceSetBuffer(){} thinInstanceRefreshBoundingInfo(){} setEnabled(v){ this.isVisible=v; } },
  ArcRotateCamera: class ArcRotateCamera { constructor(){ this.alpha=0; this.beta=0; this.radius=25; this.position=new BABYLON.Vector3(); this.target=new BABYLON.Vector3(); } attachControl(){ /* noop */ } detachControl(){ /* noop */ } getFrontPosition(){ return new BABYLON.Vector3(0,0,1); } },
  HemisphericLight: class HemisphericLight { constructor(){ this.intensity=1; } dispose(){} },
  DirectionalLight: class DirectionalLight { constructor(){ this.intensity=1; this.direction=new BABYLON.Vector3(0,0,-1); this.position=new BABYLON.Vector3(0,0,0);} dispose(){} },
  TransformNode: class TransformNode { constructor(name){ this.name=name; this.position=new BABYLON.Vector3(); this.isPickable=false; } setEnabled(){ /* noop */ } },
    Material: { MATERIAL_ALPHABLEND: 2 },
    PointerEventTypes: { POINTERDOWN: 1 }
  };
  global.BABYLON = BABYLON;
}

// ---- Test Log Suppression Layer -------------------------------------------
// Default: reduce noisy console output so test results highlight only warnings/errors.
// Opt-out: set env TEST_VERBOSE=1 (e.g. `TEST_VERBOSE=1 npm test`).
// Strategy:
//  - Buffer frequent identical lines and only show the first occurrence plus a summary count at process exit.
//  - Suppress per-test health logs (moved to debug level) unless verbose.
//  - Preserve console.error and console.warn (still visible) but dedupe repeats.
//  - Provide global helper `__logFlush()` for ad-hoc flushing within a test if needed.
try {
  const verbose = process.env.TEST_VERBOSE === '1';
  if(!verbose){
    const orig = { log: console.log, info: console.info, warn: console.warn, error: console.error };
    const counts = new Map();
    const firstLines = new Set();
    function track(kind, args){
      const key = args.map(a=> (typeof a==='string'? a : JSON.stringify(a))).join(' | ');
      const c = counts.get(key) || { kind, n:0, line:key };
      c.n++; counts.set(key,c);
      if(!firstLines.has(key)){
        firstLines.add(key);
        // Show the very first occurrence (except per-test-health noise)
        if(!/\[per-test-health\]/.test(key)) orig[kind](...args);
      }
    }
    console.log = (...a)=> track('log', a);
    console.info = (...a)=> track('info', a);
    // Keep warnings & errors but still de-duplicate after first display
    console.warn = (...a)=> track('warn', a);
    console.error = (...a)=> track('error', a);
    global.__logFlush = ()=>{
      const summary = [];
      for(const { kind, n, line } of counts.values()){
        if(n>1) summary.push({ kind, occurrences:n, preview: line.slice(0,160) });
      }
      if(summary.length){
        orig.log('\n[TestLogSummary]', summary.length,'collapsed lines');
        for(const s of summary){ orig.log(`  [${s.kind}] x${s.occurrences}: ${s.preview}`); }
      }
    };
    process.on('exit', ()=>{ try { global.__logFlush(); } catch{} });
  }
} catch {}
// ---------------------------------------------------------------------------
