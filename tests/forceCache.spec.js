import { jest } from '@jest/globals';

// Tests deterministic versioned force cache: repeated force computations without geometry change
// should perform only a single network /simple_calculate call; after geometry mutation & markPositionsChanged,
// the next compute should perform a new network call.

describe('force cache: repeated force computations without geometry change hit cache', () => {
  test('second computeForces no network; after mutation triggers network again', async () => {
    const calls = [];
    // Simple deterministic response; energy increments slightly each request so we can detect if fetch happened.
    let energyCounter = -12.0;
    function nextEnergy(){ energyCounter += 0.1; return Number(energyCounter.toFixed(4)); }
    global.fetch = jest.fn(async (url, opts) => {
      calls.push(url);
  if(url.endsWith('/serve/simple')){
        const E = nextEnergy();
        return { ok:true, status:200, json: async ()=> ({ results:{ energy: E, forces:[[0.2,0,0]], stress:null } }) };
      }
      if(url.endsWith('/relax')){
        // Provide a minimal relax response in case code inadvertently calls it
        return { ok:true, status:200, json: async ()=> ({ initial_energy: nextEnergy(), final_energy: nextEnergy(), positions:[[0,0,0]], forces:[[0.2,0,0]], steps_completed:1 }) };
      }
      throw new Error('Unexpected URL '+url);
    });

    // Minimal DOM/canvas stubs required by index.js energy canvas logic
    global.document = { getElementById: ()=> null };
    global.window = Object.assign(global.window||{}, {});
    const canvas = { addEventListener: ()=>{}, getBoundingClientRect: ()=>({ left:0, top:0 }) };
    class DummyEngine { runRenderLoop(fn){ /* ignore render loop for tests */ } }
    class DummyScene { constructor(engine){ this._engine=engine; this.meshes=[]; this.onBeforeRenderObservable={ add:()=>{} }; } render(){} }
    global.BABYLON = {
      Engine: DummyEngine,
      Scene: DummyScene,
      TransformNode: function(){},
      MeshBuilder: {
        CreateCylinder: ()=> ({ dispose(){}, position:{}, rotationQuaternion:null, scaling:{}, setEnabled(){}, material:null, thinInstanceSetBuffer:()=>{} }),
        CreateSphere: ()=> ({ dispose(){}, position:{}, rotationQuaternion:null, scaling:{}, setEnabled(){}, material:null, thinInstanceSetBuffer:()=>{} })
      },
      Matrix: { Compose: (s,q,t)=> ({ m: new Float32Array(16) }) },
      StandardMaterial: function(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; },
      Color3: function(r=0,g=0,b=0){ this.r=r; this.g=g; this.b=b; this.clone=()=>new global.BABYLON.Color3(this.r,this.g,this.b); this.scale=(s)=>new global.BABYLON.Color3(this.r*s,this.g*s,this.b*s); },
      Color4: function(){},
      ArcRotateCamera: function(){ this.attachControl=()=>{}; },
      HemisphericLight: function(){},
      Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new this(0,1,0);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new this(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);} normalize(){ const m=Math.hypot(this.x,this.y,this.z)||1; this.x/=m; this.y/=m; this.z/=m; return this;} },
      Quaternion: { Identity: ()=>({}), RotationAxis: ()=>({}) }
    };

    // Mock modules used by index.js
    jest.unstable_mockModule('../public/render/moleculeView.js', () => ({ createMoleculeView: () => ({ rebuildBonds: ()=>{} }) }));
    jest.unstable_mockModule('../public/domain/moleculeState.js', () => ({ createMoleculeState: ({ elements, positions }) => ({ elements, positions: positions.map(p=>({x:p[0],y:p[1],z:p[2]})), bonds:[], bus:{ on:()=>{} }, markPositionsChanged(){}, dynamics:{} }) }));
    jest.unstable_mockModule('../public/domain/bondService.js', () => ({ createBondService: () => ({ recomputeAndStore: ()=>[] }) }));
    jest.unstable_mockModule('../public/domain/selectionService.js', () => ({ createSelectionService: () => ({}) }));
    jest.unstable_mockModule('../public/core/pickingService.js', () => ({ createPickingService: ()=> ({}) }));
    jest.unstable_mockModule('../public/domain/manipulationService.js', () => ({ createManipulationService: ()=> ({ beginDrag:()=>{}, updateDrag:()=>{}, endDrag:()=>{}, setDragPlane:()=>{}, rotateBond:()=>{} }) }));
    jest.unstable_mockModule('../public/vr/setup.js', () => ({ createVRSupport: ()=> ({ init: async ()=> ({ supported:false }) }) }));
    jest.unstable_mockModule('../public/vr/vr-picker.js', () => ({ createVRPicker: ()=> ({}) }));
    jest.unstable_mockModule('../public/fairchem_provider.js', () => ({ createFairChemForcefield: ()=> ({}) }));
    jest.unstable_mockModule('../public/util/funcCount.js', () => ({ __count: ()=>{} }));

    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(canvas, { elements:[8], positions:[[0,0,0]], bonds:[] });
    // Clear any initial baseline /simple_calculate call(s)
    calls.length = 0;

    // Because initNewViewer performs a baseline force fetch, cache is now primed.
    // First explicit compute should be a cache hit: expect zero network calls.
    await api.ff.computeForces({ sync:true });
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(0);

    // Mutate geometry to invalidate cache
    api.state.positions[0].x += 0.5; api.state.markPositionsChanged();
    calls.length = 0;
    await api.ff.computeForces({ sync:true });
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(1);

    // Second call after invalidation (no further mutation) should be cache hit again
    calls.length = 0;
    await api.ff.computeForces({ sync:true });
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(0);
    api.shutdown && api.shutdown();
  });
});
