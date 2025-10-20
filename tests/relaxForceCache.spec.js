import { getWS } from '../public/fairchem_ws_client.js';

// WS-only integration: a relax step seeds the force cache with returned forces/energy.
// Subsequent computeForces without geometry change should keep the cache version stable.
// After a geometry mutation, the next compute should trigger an idle compute and bump version.

describe('relax + force cache integration (WS)', () => {
  test('relaxStep seeds cache; compute stable; mutation -> idle compute bumps version', async () => {
    if (!global.document) global.document = { getElementById: ()=> null };
    if (!global.window) global.window = {};
    if (!window.location) window.location = { protocol:'http:', host:'localhost:8000', origin:'http://localhost:8000' };
    class DummyEngine { runRenderLoop(fn){} stopRenderLoop(){} }
    class DummyScene { constructor(engine){ this._engine=engine; this.meshes=[]; this.onBeforeRenderObservable={ add:()=>{} }; } render(){} }
    global.BABYLON = {
      Engine: DummyEngine,
      Scene: DummyScene,
      TransformNode: function(){},
      MeshBuilder: { CreateCylinder: ()=>({ dispose(){}, setEnabled(){}, thinInstanceSetBuffer(){}, position:{}, rotationQuaternion:null, scaling:{} }), CreateSphere: ()=>({ dispose(){}, setEnabled(){}, thinInstanceSetBuffer(){}, position:{}, rotationQuaternion:null, scaling:{} }) },
      Matrix: { Compose: ()=> ({ m:new Float32Array(16) }) },
      Quaternion: { Identity: ()=>({}), RotationAxis: ()=>({}) },
      StandardMaterial: function(){ this.diffuseColor=null; this.emissiveColor=null; this.specularColor=null; },
      Color3: function(r=1,g=1,b=1){ this.r=r; this.g=g; this.b=b; this.clone=()=>new global.BABYLON.Color3(this.r,this.g,this.b); this.scale=(s)=>new global.BABYLON.Color3(this.r*s,this.g*s,this.b*s); },
      Color4: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; }, HemisphericLight: function(){},
      Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }
    };

    const origWS = global.WebSocket;
    class FakeWS { constructor(){ this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0);} send(){} close(){} onopen(){} onmessage(){} onerror(){} }
    global.WebSocket = FakeWS;
    const ws = getWS();
    // Inject frames for relax and idle computes
    ws.setTestHook((msg)=>{
      try {
        // For relax one-shot step: when START_SIMULATION sent with params, inject a relax frame
        if (msg && msg.simulationParams) {
          setTimeout(()=>{
            ws.injectTestResult({ positions:[[0,0,0]], forces:[[0.05,0,0]], energy:-4.9, temperature:300 });
          }, 0);
        }
        // For subsequent idle compute requests, inject an idle energy frame
        if (msg && msg.type && !msg.simulationParams) {
          setTimeout(()=>{ ws.injectTestResult({ energy:-4.9, forces:[[0.05,0,0]] }); }, 0);
        }
      } catch {}
    });

    const { initNewViewer } = await import('../public/index.js');
    const canvas = { addEventListener: ()=>{}, getBoundingClientRect: ()=>({ left:0, top:0 }) };
    const api = await initNewViewer(canvas, { elements:[8], positions:[[0,0,0]], bonds:[] });

    // Do a relax step: seeds cache with forces/energy
    await api.relaxStep();
    const v0 = api.getForceCacheVersion();

    // computeForces without geometry change should keep version stable
    await api.ff.computeForces({ sync:true });
    const v1 = api.getForceCacheVersion();
    expect(v1).toBe(v0);

    // mutate -> mark change
    api.state.positions[0].x += 0.1; api.state.markPositionsChanged();
    // next compute should accept idle frame and bump version
    await api.ff.computeForces({ sync:true });
    const v2 = api.getForceCacheVersion();
    expect(v2).toBeGreaterThan(v1);

    ws.setTestHook(null); global.WebSocket = origWS; api.shutdown && api.shutdown();
  });
});
