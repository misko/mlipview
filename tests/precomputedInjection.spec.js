import { jest } from '@jest/globals';

// Verifies that when a fresh force cache exists, relax & md requests include `precomputed`.

describe('precomputed injection (relax + md)', () => {
  test('includes precomputed energy/forces/stress when cache fresh; skips after geometry change', async () => {
    const calls = [];
    const bodies = [];
    // Mock fetch: first simple -> seeds cache; relax/md echo; second relax after mutation should have no precomputed.
    global.fetch = jest.fn(async (url, opts) => {
      calls.push(url);
      if (opts && opts.body) {
        try { bodies.push({ url, body: JSON.parse(opts.body) }); } catch {}
      }
      if (url.endsWith('/serve/simple')) {
        return { ok:true, status:200, json: async ()=> ({ results:{ energy: 1.234, forces:[[0.1,0.0,0.0]], stress:[0.1,0.2,0.3,0.01,0.02,0.03] } }) };
      }
      if (url.endsWith('/serve/relax')) {
        return { ok:true, status:200, json: async ()=> ({ initial_energy:1.234, final_energy:1.111, positions:[[0,0,0]], forces:[[0.05,0,0]], steps_completed:1 }) };
      }
      if (url.endsWith('/serve/md')) {
        return { ok:true, status:200, json: async ()=> ({ initial_energy:1.234, final_energy:1.200, positions:[[0,0,0]], velocities:[[0,0,0]], forces:[[0.02,0,0]], steps_completed:1, temperature:300 }) };
      }
      throw new Error('Unexpected URL '+url);
    });

    // Minimal DOM / Babylon mocks
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

    // Module mocks
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

    // Seed forces (simple_calculate) explicitly
    await api.ff.computeForces({ sync:true });

    // Perform relax step - should include precomputed
    calls.length = 0; bodies.length = 0;
    await api.relaxStep();
    const relaxReq = bodies.find(b=> b.url.endsWith('/serve/relax'));
    expect(relaxReq).toBeTruthy();
    expect(relaxReq.body.precomputed).toBeTruthy();
    expect(relaxReq.body.precomputed.energy).toBeCloseTo(1.234, 6);
    expect(relaxReq.body.precomputed.forces.length).toBe(1);
    expect(relaxReq.body.precomputed.stress.length).toBe(6);

    // Perform md step - should include precomputed again (cache still valid)
    bodies.length = 0;
    await api.mdStep();
    const mdReq = bodies.find(b=> b.url.endsWith('/serve/md'));
    expect(mdReq).toBeTruthy();
    expect(mdReq.body.precomputed).toBeTruthy();
  // After relax step the force cache energy updates to final_energy (1.111) and that is what should be injected.
  expect(mdReq.body.precomputed.energy).toBeCloseTo(1.111, 6);

    // Invalidate geometry and ensure no precomputed now
    api.state.positions[0].x += 0.5; api.state.markPositionsChanged();
    bodies.length = 0;
    await api.relaxStep();
    const relaxAfterMove = bodies.find(b=> b.url.endsWith('/serve/relax'));
    expect(relaxAfterMove).toBeTruthy();
    expect(relaxAfterMove.body.precomputed).toBeFalsy();

    api.shutdown && api.shutdown();
  });
});
