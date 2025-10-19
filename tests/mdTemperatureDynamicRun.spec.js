/** @jest-environment jsdom */
// Test: During a continuous MD run, moving the temperature slider changes subsequent WS temperatures.

// BABYLON stubs
beforeAll(()=>{ if(!global.BABYLON){ global.BABYLON={ TransformNode:class{}, MeshBuilder:{ CreateCylinder:()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial:class{ constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }}, Color3:class{}, Vector3:class{ constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} }, Quaternion:class{ static Identity(){return{};} }, Scene:class{} }; } });

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

function stubWebSocketAndHook(){
  const sent = [];
  let wsOnMessage = null;
  global.WebSocket = class {
    constructor(){ this.readyState=1; setTimeout(()=> this.onopen && this.onopen(),0); }
    set binaryType(_){}
    send(_){}
    close(){}
    onopen(){}
    onerror(){}
    onclose(){}
    set onmessage(fn){ wsOnMessage = fn; }
    get onmessage(){ return wsOnMessage; }
  };
  window.__WS_TEST_HOOK__ = (msg)=>{ sent.push(msg); };
  return { sent, emit: (obj)=>{ if (typeof window.__ON_WS_RESULT__ === 'function') window.__ON_WS_RESULT__(obj); } };
}

async function setup(){
  const ws = stubWebSocketAndHook();
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  window.__MLIP_FEATURES = { RELAX_LOOP:false, MD_LOOP:true, ENERGY_TRACE:false, FORCE_VECTORS:false };
  window.__MLIPVIEW_TEST_MODE = true;
  const hud=document.createElement('div'); hud.className='hud'; document.body.appendChild(hud);
  hud.innerHTML = `<span></span><button id="btnRelax"></button><button id="btnRelaxRun"></button><span></span><button id="btnMD"></button><button id="btnMDRun"></button><select id="forceProviderSel"></select><button id="btnCell"></button><button id="btnGhosts"></button><button id="btnToggleForces"></button><span id="status">Ready</span><select id="xrModeSelect"></select>`;
  const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
  const canvas=document.createElement('canvas'); canvas.id='viewer'; document.body.appendChild(canvas);
  const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
  const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);
  const mod = await import('../public/index.js');
  const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
  const api = await mod.initNewViewer(canvas, { elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.95,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] });
  window.viewerApi = api;
  initTemperatureSlider({ hudEl:hud, getViewer:()=>api });
  return { api, ws };
}

describe('MD run dynamic temperature', () => {
  test('temperature updates propagate mid-run', async () => {
    const { api, ws } = await setup();
  api.setMinStepInterval(1);
  api.startMDContinuous({ steps:50, temperature: window.__MLIP_TARGET_TEMPERATURE });
  // Allow async ws.ensureConnected + init + startSimulation to run
  await Promise.resolve().then(()=>{});
  await new Promise(r=>setTimeout(r,0));
    // Simulate a few frames at initial temperature
    for(let i=0;i<5;i++){ ws.emit({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [], energy: 0, temperature: window.__MLIP_TARGET_TEMPERATURE }); }
    const firstMsgs = ws.sent.slice();
    const initialTemps = firstMsgs.filter(m=> m.simulationParams).map(m=> m.simulationParams.temperature).filter(t=> typeof t==='number');
    expect(initialTemps.length).toBeGreaterThan(0);
    const uniqueFirst = new Set(initialTemps);
    expect(uniqueFirst.size).toBe(1);
  // Change slider to max; restart MD so a new START_SIMULATION reflects the updated temperature
    const slider = document.getElementById('mdTempSlider');
    slider.value = slider.max; slider.dispatchEvent(new Event('input'));
    const newTarget = window.__MLIP_TARGET_TEMPERATURE;
  api.stopSimulation();
  api.startMDContinuous({ steps: 20, temperature: newTarget });
  await Promise.resolve().then(()=>{});
  await new Promise(r=>setTimeout(r,0));
  for(let i=0;i<5;i++){ ws.emit({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [], energy: 0, temperature: newTarget }); }
  const later = ws.sent.slice(firstMsgs.length);
  const startMsgs = later.filter(m=> m && m.simulationParams && typeof m.simulationParams.temperature==='number');
  const laterTemps = startMsgs.map(m=> m.simulationParams.temperature);
  expect(laterTemps.some(t=> t === newTarget)).toBe(true);
  }, 10000);
});
