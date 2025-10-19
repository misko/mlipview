/** @jest-environment jsdom */
// WS migration: a single relaxStep should issue exactly one START_SIMULATION (relax) action and no SIMPLE_CALCULATE.

beforeAll(()=>{
  // Minimal BABYLON stubs
  global.BABYLON = global.BABYLON || { Engine: function(){ this.runRenderLoop=()=>{}; }, Scene: function(){ this.onPointerObservable={}; this.render=()=>{}; }, Color3: function(){}, MeshBuilder:{ CreateSphere: ()=>({}) }, StandardMaterial: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; } };
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ render:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} }) }));
import { stubWebSocketAndHook } from './utils/wsTestStub.js';

describe('network: single relax step (WS)', () => {
  test('relaxStep triggers exactly one START_SIMULATION(relax) and zero SIMPLE_CALCULATE', async () => {
    document.body.innerHTML = '<canvas id="viewer"></canvas><div class="hud"></div>';
    const ws = stubWebSocketAndHook();
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const api = await initNewViewer(canvas, { elements:['O'], positions:[{x:0,y:0,z:0}], bonds:[] });
    // Trigger a single relax step and immediately emit a frame to resolve the one-shot
    const p = api.relaxStep();
    ws.emit({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [[0,0,0]], energy: -9.9 });
    await p;
    const sent = ws.sent || [];
    const relaxStarts = sent.filter(m=> m && m.type != null && m.simulationType != null && m.simulationParams && m.simulationParams.optimizer === 'bfgs');
    const simpleMsgs = sent.filter(m=> m && m.type != null && m.type === 'SIMPLE_CALCULATE');
    expect(relaxStarts.length).toBe(1);
    expect(simpleMsgs.length).toBe(0);
  });
});
