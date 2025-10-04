import { attachConsistentLighting } from '../public/render/lighting.js';

// Minimal BABYLON stubs used by lighting helper
if(!global.BABYLON) global.BABYLON = {};
if(!BABYLON.Vector3){
  BABYLON.Vector3 = class Vector3 {
    constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; }
    static Zero(){ return new BABYLON.Vector3(0,0,0); }
    copyFrom(v){ this.x=v.x; this.y=v.y; this.z=v.z; return this; }
    subtractInPlace(v){ this.x-=v.x; this.y-=v.y; this.z-=v.z; return this; }
    normalize(){ const l=Math.hypot(this.x,this.y,this.z)||1; this.x/=l; this.y/=l; this.z/=l; return this; }
    set(x,y,z){ this.x=x; this.y=y; this.z=z; return this; }
  };
}
if(!BABYLON.HemisphericLight){ BABYLON.HemisphericLight = class { constructor(name,dir){ this.name=name; this.direction=dir; this.intensity=1; } dispose(){} }; }
if(!BABYLON.DirectionalLight){ BABYLON.DirectionalLight = class { constructor(name,dir){ this.name=name; this.direction=dir; this.intensity=1; this.position=new BABYLON.Vector3(); } dispose(){} }; }
if(!BABYLON.Color4){ BABYLON.Color4 = class { constructor(){} }; }

function makeSceneStub(){
  const observers=[];
  return {
    onBeforeRenderObservable:{ add(fn){ observers.push(fn); return fn; }, remove(fn){ const i=observers.indexOf(fn); if(i>=0) observers.splice(i,1); } },
    tick(){ observers.forEach(fn=>fn()); },
    getLightByName(){ return null; },
    activeCamera:null
  };
}

describe('lighting consistency', () => {
  test('creates ambient and headlight with configured intensities', () => {
    const scene = makeSceneStub();
    const camera = { alpha:0, beta:Math.PI/3, position:new BABYLON.Vector3(0,0,10) };
    scene.activeCamera = camera;
    const handle = attachConsistentLighting(scene, camera, { ambientIntensity:0.15, directionalIntensity:0.85 });
    expect(handle).toBeTruthy();
    expect(handle.ambient).toBeTruthy();
    expect(handle.head).toBeTruthy();
    expect(handle.ambient.intensity).toBeCloseTo(0.15, 5);
    expect(handle.head.intensity).toBeCloseTo(0.85, 5);
  });
});
