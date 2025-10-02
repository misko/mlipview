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
  Vector3: class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } add(v){ return new BABYLON.Vector3(this.x+v.x,this.y+v.y,this.z+v.z);} addInPlace(v){ this.x+=v.x; this.y+=v.y; this.z+=v.z; return this;} scale(s){ return new BABYLON.Vector3(this.x*s,this.y*s,this.z*s);} length(){ return Math.hypot(this.x,this.y,this.z);} negate(){ return new BABYLON.Vector3(-this.x,-this.y,-this.z);} clone(){ return new BABYLON.Vector3(this.x,this.y,this.z);} subtract(v){ return new BABYLON.Vector3(this.x-v.x,this.y-v.y,this.z-v.z);} normalize(){ const L=this.length()||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L);} normalizeToNew(){ return this.normalize(); } toString(){ return `{X: ${this.x} Y: ${this.y} Z: ${this.z}}`; } static Zero(){ return new BABYLON.Vector3(0,0,0);} static Up(){ return new BABYLON.Vector3(0,1,0);} static Right(){ return new BABYLON.Vector3(1,0,0);} static Dot(a,b){ return a.x*b.x + a.y*b.y + a.z*b.z; } static Cross(a,b){ return new BABYLON.Vector3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);} },
  Quaternion: class Quaternion { static Identity(){ return new BABYLON.Quaternion(); } static RotationAxis(){ return new BABYLON.Quaternion(); } },
    Matrix: class Matrix { constructor(){ this.m=new Float32Array(16);} static Compose(scale,rot,pos){ const m=new BABYLON.Matrix(); m.m[12]=pos.x; m.m[13]=pos.y; m.m[14]=pos.z; return m;} clone(){ const c=new BABYLON.Matrix(); c.m.set(this.m); return c;} copyToArray(arr,off){ for(let i=0;i<16;i++) arr[off+i]=this.m[i]||0; } },
    Color3: class Color3 { constructor(r,g,b){ this.r=r; this.g=g; this.b=b; } clone(){ return new BABYLON.Color3(this.r,this.g,this.b);} scale(f){ return new BABYLON.Color3(this.r*f,this.g*f,this.b*f);} },
    Color4: class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a;} },
    StandardMaterial: class StandardMaterial { constructor(){ this.alpha=1; this.diffuseColor=new BABYLON.Color3(1,1,1); this.emissiveColor=new BABYLON.Color3(0,0,0);} },
    MeshBuilder: { CreateSphere(name,opts,scene){ return new BABYLON.Mesh(name); }, CreateCylinder(name,opts,scene){ return new BABYLON.Mesh(name);} },
    Mesh: class Mesh { constructor(name){ this.name=name; this.material=null; this.isPickable=false; this.thinInstanceEnablePicking=false; this.alwaysSelectAsActiveMesh=false; this.isVisible=true; this.scaling=new BABYLON.Vector3(1,1,1); this.position=new BABYLON.Vector3(); } dispose(){} thinInstanceSetBuffer(){} thinInstanceRefreshBoundingInfo(){} setEnabled(v){ this.isVisible=v; } },
  TransformNode: class TransformNode { constructor(name){ this.name=name; this.position=new BABYLON.Vector3(); this.isPickable=false; } setEnabled(){ /* noop */ } },
    Material: { MATERIAL_ALPHABLEND: 2 }
  };
  global.BABYLON = BABYLON;
}
