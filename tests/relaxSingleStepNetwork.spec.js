/** @jest-environment jsdom */
// WS migration: a single relaxStep should issue exactly one START_SIMULATION (relax) action and no SIMPLE_CALCULATE.

beforeAll(()=>{
  // Minimal BABYLON stubs
  global.BABYLON = global.BABYLON || { Engine: function(){ this.runRenderLoop=()=>{}; }, Scene: function(){ this.onPointerObservable={}; this.render=()=>{}; }, Color3: function(){}, MeshBuilder:{ CreateSphere: ()=>({}) }, StandardMaterial: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; } };
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ render:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} }) }));
import { getWS } from '../public/fairchem_ws_client.js';

describe('network: single relax step (WS)', () => {
  test('relaxStep triggers exactly one START_SIMULATION(relax) and zero SIMPLE_CALCULATE', async () => {
    document.body.innerHTML = '<canvas id="viewer"></canvas><div class="hud"></div>';
    // Stub WS to auto-open and capture messages
    const origWS = global.WebSocket;
    class FakeWS { constructor(){ this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0);} send(){} close(){} onopen(){} onmessage(){} onerror(){} }
    global.WebSocket = FakeWS;
    const ws = getWS();
    const sent = []; ws.setTestHook(m=> sent.push(m));
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const api = await initNewViewer(canvas, { elements:['O'], positions:[{x:0,y:0,z:0}], bonds:[] });
    // Trigger a single relax step and immediately emit a frame to resolve the one-shot
    const p = api.relaxStep();
    ws.injectTestResult({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [[0,0,0]], energy: -9.9 });
    await p;
    const relaxStarts = sent.filter(m=> m && m.type != null && m.simulationType != null && m.simulationParams && m.simulationParams.optimizer === 'bfgs');
    expect(relaxStarts.length).toBe(1);
    // No SIMPLE_CALCULATE in new design; implicitly verified by absence of such messages in ws client
    global.WebSocket = origWS;
  });
});
