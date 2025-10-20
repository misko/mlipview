import { getWS } from '../public/fairchem_ws_client.js';

// Force cache is now client-managed (versioned) but transport is WS-only.
// This test verifies: (1) initial baseline seeds cache with an idle energy frame,
// (2) repeated computeForces without geometry changes do not change the cache version,
// (3) after geometry mutation and markPositionsChanged, a new idle compute occurs and cache version increments.

describe('force cache: version increments only on geometry changes (WS idle compute)', () => {
  test('baseline -> computeForces (no change) -> mutate -> computeForces triggers idle compute', async () => {
    // Minimal DOM setup for energy canvas (optional, not asserting visuals here)
    if (!global.document) global.document = { getElementById: ()=> null };
    if (!global.window) global.window = {};
    if (!window.location) window.location = { protocol:'http:', host:'localhost:8000', origin:'http://localhost:8000' };
    // Babylon minimal mocks
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

    // WS stub and hooks: inject idle energy/forces on demand when computeForces requests it.
    const origWS = global.WebSocket;
    class FakeWS { constructor(){ this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0);} send(){} close(){} onopen(){} onmessage(){} onerror(){} }
    global.WebSocket = FakeWS;
    const ws = getWS();
    // For this test, we inject an energy frame when the app does a USER_INTERACTION (idle compute), and do nothing otherwise
    ws.setTestHook((msg)=>{
      try {
        if (msg && msg.type) {
          // Defer injection to allow waitForEnergy to attach its listener
          setTimeout(()=>{
            const energy = -12.0; const forces = [[0.2,0,0]]; const positions = [[0,0,0]];
            ws.injectTestResult({ energy, forces, positions });
          }, 0);
        }
      } catch {}
    });

    const { initNewViewer } = await import('../public/index.js');
    const canvas = { addEventListener: ()=>{}, getBoundingClientRect: ()=>({ left:0, top:0 }) };
    const api = await initNewViewer(canvas, { elements:[8], positions:[[0,0,0]], bonds:[] });

    // Initial baseline called by initNewViewer may have already seeded cache
    const v0 = api.getForceCacheVersion();
    // Explicit compute without geometry change should not change the cache version
    await api.ff.computeForces({ sync:true });
    const v1 = api.getForceCacheVersion();
    expect(v1).toBe(v0);

    // Mutate geometry to invalidate and mark change
    api.state.positions[0].x += 0.25; api.state.markPositionsChanged();
    // Next compute should accept a new idle frame and bump the cache version
    await api.ff.computeForces({ sync:true });
    const v2 = api.getForceCacheVersion();
    expect(v2).toBeGreaterThan(v1);

    // Another compute without changes should keep version stable
    await api.ff.computeForces({ sync:true });
    const v3 = api.getForceCacheVersion();
    expect(v3).toBe(v2);

    ws.setTestHook(null); global.WebSocket = origWS; api.shutdown && api.shutdown();
  });
});
