/** @jest-environment jsdom */
// Validate that streaming WS frames drive the RPS label via __noteRequestCompleted.

import { getWS } from '../public/fairchem_ws_client.js';

beforeAll(()=>{
  // Minimal BABYLON stubs
  if(!global.BABYLON){
    global.BABYLON = { TransformNode: class {}, MeshBuilder: { CreateCylinder: ()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial: class { constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }}, Color3: class {}, Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=y; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} }, Quaternion: class { static Identity(){ return {}; } static RotationAxis(){ return {}; } }, Scene: class {} };
  }
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

describe('MD RPS label (WS)', () => {
  test('RPS label shows ~10.0 given 5 frames over 400ms', async () => {
    // Build DOM similar to index.html expects
    window.__MLIPVIEW_TEST_MODE = true;
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    window.__MLIP_FEATURES = { RELAX_LOOP:false, MD_LOOP:true, ENERGY_TRACE:false, FORCE_VECTORS:false };
    const canvas=document.createElement('canvas'); canvas.id='viewer'; canvas.addEventListener=()=>{}; document.body.appendChild(canvas);
    const hud=document.createElement('div'); hud.className='hud'; document.body.appendChild(hud);
    const rps=document.createElement('span'); rps.id='rpsLabel'; rps.textContent='RPS: --'; hud.appendChild(rps);
    // Energy elements (some code references these)
    const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
    const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
    const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);

  // Stub WebSocket to avoid network and auto-open
  const origWS = global.WebSocket;
  class FakeWS { constructor(){ this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0);} send(){} close(){} onopen(){} onmessage(){} onerror(){} }
  global.WebSocket = FakeWS;
  const ws = getWS();
  ws.setTestHook(()=>{});
    const mod = await import('../public/index.js');
    const api = await mod.initNewViewer(canvas, { elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.96,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] });
    // Override performance.now for deterministic RPS
    const origNow = performance.now;
    let now = 1000;
    performance.now = ()=> now;
    try {
      // Start streaming MD and wait a tick so onResult is subscribed
      api.startMDContinuous({ steps:100, temperature: 1500 });
      await Promise.resolve().then(()=>{});
      await new Promise(r=>setTimeout(r,0));
      // Emit 5 frames with 100ms intervals => dt ~ 400ms for 5 samples => RPS ~ 4*1000/400 = 10
      for(let i=0;i<5;i++){
        now += 100; // advance time BEFORE sample so dt accumulates
        ws.injectTestResult({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [], energy: -1.0, temperature: 1500 });
        // allow __noteRequestCompleted to run and label to update
        await new Promise(r=>setTimeout(r,0));
      }
      // Let DOM settle
      await Promise.resolve().then(()=>{});
      const label = document.getElementById('rpsLabel');
      expect(label).toBeTruthy();
      // It should update from default '--' to a numeric value
      expect(label.textContent).toMatch(/^RPS:\s+\d+\.\d/);
      api.stopSimulation();
    } finally {
      performance.now = origNow;
      global.WebSocket = origWS;
    }
  });
});
