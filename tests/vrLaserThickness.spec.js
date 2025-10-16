import { createVRSupport } from '../public/vr/setup.js';

// Minimal Babylon stubs (only what the feature/test will touch). If real Babylon exists these are ignored.
if(!global.BABYLON){
  class Vec3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } clone(){ return new Vec3(this.x,this.y,this.z); } subtract(v){ return new Vec3(this.x-v.x,this.y-v.y,this.z-v.z); } add(v){ return new Vec3(this.x+v.x,this.y+v.y,this.z+v.z); } scale(f){ return new Vec3(this.x*f,this.y*f,this.z*f); } length(){ return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z); } normalize(){ const l=this.length()||1; this.x/=l; this.y/=l; this.z/=l; return this; } };
  class Ray { constructor(o,d,l){ this.origin=o; this.direction=d; this.length=l; } }
  const Axis={ X:new Vec3(1,0,0), Y:new Vec3(0,1,0), Z:new Vec3(0,0,1) };
  global.BABYLON={ Vector3:Vec3, Ray, Axis, MeshBuilder:{ CreateCylinder:(name, opts)=>({ name, opts, scaling:{x:1,y:1,z:1}, position:new Vec3(), parent:null, material:{}, dispose(){} }) }, Quaternion: class { static Identity(){ return {x:0,y:0,z:0,w:1}; } }, Color4: class{ constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } } };
}

// Build a simple onBeforeRenderObservable harness
function makeScene(pickSequence){
  const observers=[];
  let frame=0;
  // Provide a fake masters accessor used indirectly (simulate presence of masters to fully run logic)
  const fakeMaster = { name:'molecule_root', rotationQuaternion:null, scaling:{x:1,y:1,z:1}, position:{x:0,y:0,z:0}, refreshBoundingInfo(){} };
  const scene = {
    meshes:[],
    onBeforeRenderObservable:{ add(fn){ observers.push(fn); } },
    pickWithRay:(ray)=>{
      if(frame < pickSequence.length){ const res = pickSequence[frame]; frame++; return res; }
      frame++; return { hit:false };
    },
    // Simulate a method used by vr-utils to retrieve masters (simplified)
    getMeshByName(name){ return name==='molecule_root'?fakeMaster:null; },
    step(n=1){ for(let i=0;i<n;i++){ observers.slice().forEach(fn=>fn()); } },
    createDefaultXRExperienceAsync: async ()=>({ input:{ controllers:[], onControllerAddedObservable:{ add(){ /* noop in this test */ } } }, baseExperience:{ sessionManager:{} } })
  };
  return scene;
}

// Utility to fabricate a controller after init shim path
function fabricateController(vr){
  const controllers = vr.controllers();
  if(controllers.length) return controllers[0];
  const fake = { uniqueId:'shimControllerLaser', pointer:{ getAbsolutePosition:()=>new BABYLON.Vector3(0,1,0), getDirection:()=>new BABYLON.Vector3(0,0,1) } };
  controllers.push(fake);
  return fake;
}

// Sequence: first a few misses (no hover), then a bond-like hit, then miss again.
// We'll simulate a hit by returning { hit:true, pickedMesh:{} } which our resolver maps to an atom.
const atomPick = { hit:true, pickedMesh:{} };

// A view resolver that returns an atom for any hit pick result we pass in (simplified)
const view = { resolveAtomPick:(p)=> p===atomPick? { kind:'atom', index:1 }: null, resolveBondPick:()=> null };

const selectionService = { clickAtom(){}, clickBond(){}, clear(){} };

// Test: laser diameter increases on hover and relaxes after

test.skip('VR laser thickness increases on atom hover', async () => {
  // pick sequence: 4 frames miss, 6 frames hit, 6 frames miss
  const pickSeq = [ {hit:false}, {hit:false}, {hit:false}, {hit:false}, atomPick, atomPick, atomPick, atomPick, atomPick, atomPick, {hit:false}, {hit:false}, {hit:false}, {hit:false}, {hit:false}, {hit:false} ];
  const scene = makeScene(pickSeq);
  const picking = { view, selectionService };
  const vr = createVRSupport(scene, { picking });
  const initRes = await vr.init();
  expect(initRes.supported).toBe(true);
  // Fabricate controller and run a frame so laser state initializes
  fabricateController(vr);
  scene.step(1);
  // We expect getLaserInfo to exist after feature implementation; for now it may be undefined.
  // Step frames to allow laser creation and hover detection.
  scene.step(4); // warm-up frames (no hover)
  const infoPre = vr.getLaserInfo ? vr.getLaserInfo() : [];
  const baseDiam = infoPre[0]?.diameter;

  scene.step(6); // hover frames
  const infoHover = vr.getLaserInfo ? vr.getLaserInfo() : [];
  const hoverDiam = infoHover[0]?.diameter;

  scene.step(6); // post-hover frames
  const infoPost = vr.getLaserInfo ? vr.getLaserInfo() : [];
  const postDiam = infoPost[0]?.diameter;

  // Expectations (will fail until implemented)
  expect(baseDiam).toBeDefined();
  expect(hoverDiam).toBeDefined();
  // Hover diameter should be larger than base
  expect(hoverDiam).toBeGreaterThan(baseDiam);
  // After leaving hover it should decrease (allow small tolerance for smoothing)
  expect(postDiam).toBeLessThan(hoverDiam);
});
