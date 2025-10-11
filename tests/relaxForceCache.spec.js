import { jest } from '@jest/globals';

// Integration: relaxStep should seed force cache; subsequent computeForces no network until geometry changes.

describe.skip('relax + force cache integration [cache removed]', () => {
  test('relaxStep seeds cache; next force compute hits cache; mutation triggers new fetch', async () => {
    const calls = [];
    global.fetch = jest.fn(async (url, opts) => {
      calls.push(url);
  if(url.endsWith('/serve/relax')){
        return { ok:true, status:200, json: async ()=> ({ initial_energy:-5, final_energy:-4.9, positions:[[0,0,0]], forces:[[0.05,0,0]], steps_completed:1 }) };
      }
  if(url.endsWith('/serve/simple')){
        return { ok:true, status:200, json: async ()=> ({ results:{ energy:-4.9, forces:[[0.05,0,0]] } }) };
      }
      throw new Error('Unexpected URL '+url);
    });

    global.document = { getElementById: ()=> null };
    global.window = Object.assign(global.window||{}, {});
    const canvas = { addEventListener: ()=>{}, getBoundingClientRect: ()=>({ left:0, top:0 }) };
    class DummyEngine { runRenderLoop(fn){} stopRenderLoop(){} }
    class DummyScene { constructor(engine){ this._engine=engine; this.meshes=[]; this.onBeforeRenderObservable={ add:()=>{} }; } render(){} }
    global.BABYLON = {
      Engine: DummyEngine,
      Scene: DummyScene,
      TransformNode: function(){},
      MeshBuilder: { CreateCylinder: ()=>({ dispose(){}, position:{}, rotationQuaternion:null, scaling:{}, setEnabled(){}, material:null, thinInstanceSetBuffer:()=>{} }), CreateSphere: ()=>({ dispose(){}, position:{}, rotationQuaternion:null, scaling:{}, setEnabled(){}, material:null, thinInstanceSetBuffer:()=>{} }) },
      Matrix: { Compose: ()=> ({ m:new Float32Array(16) }) },
      StandardMaterial: function(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; },
      Color3: function(r=0,g=0,b=0){ this.r=r; this.g=g; this.b=b; this.clone=()=>new global.BABYLON.Color3(this.r,this.g,this.b); this.scale=(s)=>new global.BABYLON.Color3(this.r*s,this.g*s,this.b*s); },
      Color4: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; }, HemisphericLight: function(){},
      Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new this(0,1,0);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new this(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);} normalize(){ const m=Math.hypot(this.x,this.y,this.z)||1; this.x/=m; this.y/=m; this.z/=m; return this;} },
      Quaternion: { Identity: ()=>({}), RotationAxis: ()=>({}) }
    };

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
    // Clear baseline call(s) (there will have been an initial compute)
    calls.length = 0;

    // Run relax step: expect one /relax call only
    await api.relaxStep();
  const relaxCalls = calls.filter(u=>u.endsWith('/serve/relax'));
    expect(relaxCalls.length).toBe(1);
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(0);

    // Force compute should hit cache (seeded by relax forces) => no simple_calculate
    calls.length = 0;
    await api.ff.computeForces({ sync:true });
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(0);

    // Mutate geometry to invalidate
    api.state.positions[0].x += 0.25; api.state.markPositionsChanged();
    calls.length = 0;
    await api.ff.computeForces({ sync:true });
  expect(calls.filter(u=>u.endsWith('/serve/simple')).length).toBe(1);

    api.shutdown && api.shutdown();
  });
});
