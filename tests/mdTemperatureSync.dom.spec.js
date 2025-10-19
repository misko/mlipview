/** @jest-environment jsdom */
// Test: Temperature slider and MD parameter stay in sync on initial load and after updates (WS-only).

beforeAll(()=>{
  // Minimal BABYLON stubs
  global.BABYLON = global.BABYLON || { Engine: function(){ this.runRenderLoop=()=>{}; }, Scene: function(){ this.onPointerObservable={}; this.render=()=>{}; }, Color3: function(){}, MeshBuilder:{ CreateSphere: ()=>({}) }, StandardMaterial: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; } };
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ render:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} }) }));

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
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"></div>`;
  const ws = stubWebSocketAndHook();
  const { initNewViewer } = await import('../public/index.js');
  const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
  const canvas = document.getElementById('viewer');
  const hud = document.querySelector('.hud');
  const api = await initNewViewer(canvas, { elements:['O','H','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0},{x:0,y:1,z:0}], bonds:[] });
  window.viewerApi = api;
  initTemperatureSlider({ hudEl: hud, getViewer: ()=> api });
  return { api, ws };
}

describe('Temperature slider sync', () => {
  test('initial load uses 1500K and mdStep payload matches', async () => {
    const { api, ws } = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    expect(slider).toBeTruthy();
    expect(label.textContent).toContain('1500');
    // Trigger an MD step; verify outgoing START_SIMULATION temperature
    await Promise.resolve().then(()=>{}); // tick for ws connect
    const p = api.mdStep({ temperature: window.__MLIP_TARGET_TEMPERATURE });
    // Simulate a response frame to resolve mdStep
    ws.emit({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [], energy: -1.0, temperature: window.__MLIP_TARGET_TEMPERATURE });
    await p;
    const last = (window.__WS_TEST_HOOK__ && typeof window.__WS_TEST_HOOK__==='function') ? null : null;
    const sent = ws.sent || [];
    const lastStart = sent.reverse().find(m=> m && m.type != null);
    expect(sent.length).toBeGreaterThan(0);
    // Find a START_SIMULATION with simulationParams.temperature
    const startMsg = sent.find(m=> m.simulationParams && typeof m.simulationParams.temperature === 'number');
    expect(startMsg).toBeTruthy();
    expect(startMsg.simulationParams.temperature).toBe(window.__MLIP_TARGET_TEMPERATURE);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(1500);
  });

  test('moving slider updates __MLIP_TARGET_TEMPERATURE and mdStep payload follows', async () => {
    const { api, ws } = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    // Move slider to max (near 3000)
    slider.value = slider.max;
    slider.dispatchEvent(new Event('input'));
    const target = window.__MLIP_TARGET_TEMPERATURE;
    expect(label.textContent).toContain(String(target));
    const p = api.mdStep({ temperature: target });
    ws.emit({ positions: api.state.positions.map(p=>[p.x,p.y,p.z]), forces: [], energy: -1.0, temperature: target });
    await p;
    const sent = ws.sent || [];
    const startMsg = sent.find(m=> m.simulationParams && typeof m.simulationParams.temperature === 'number');
    expect(startMsg).toBeTruthy();
    expect(startMsg.simulationParams.temperature).toBe(target);
  });
});
