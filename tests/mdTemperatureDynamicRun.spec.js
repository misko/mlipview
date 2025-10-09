/** @jest-environment jsdom */
// Test: During a continuous MD run, moving the temperature slider changes subsequent /serve/md temperatures.
import http from 'http';

// BABYLON stubs
beforeAll(()=>{ if(!global.BABYLON){ global.BABYLON={ TransformNode:class{}, MeshBuilder:{ CreateCylinder:()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial:class{ constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }}, Color3:class{}, Vector3:class{ constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} }, Quaternion:class{ static Identity(){return{};} }, Scene:class{} }; } });

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

const mdBodies=[];
global.fetch = async (url, opts={}) => {
  if(typeof url === 'string' && (/\/serve\/md$/.test(url) || /\/serve\/md_from_cache$/.test(url))){
    try { mdBodies.push(JSON.parse(opts.body||'{}')); } catch { mdBodies.push({ parseError:true }); }
    const lastT = mdBodies[mdBodies.length-1].temperature;
    const body = JSON.stringify({ positions:[[0,0,0],[0.9,0,0],[-0.2,0.9,0]], velocities:[[0,0,0],[0,0,0],[0,0,0]], forces:[[0,0,0],[0,0,0],[0,0,0]], final_energy:-1, steps_completed:1, temperature:lastT });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/simple$/.test(url)){
    const body = JSON.stringify({ energy:-1.1, forces:[[0,0,0]], positions:[[0,0,0]] });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/relax$/.test(url)){
    const body = JSON.stringify({ positions:[[0,0,0],[0.9,0,0],[-0.2,0.9,0]], forces:[[0,0,0],[0,0,0],[0,0,0]], final_energy:-1.2 });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/health$/.test(url)){
    const body = JSON.stringify({ ok:true });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  throw new Error('Unexpected fetch '+url);
};

async function setup(){
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
  return api;
}

describe('MD run dynamic temperature', () => {
  test('temperature updates propagate mid-run', async () => {
    const api = await setup();
    mdBodies.length = 0;
    // Start a long MD run (many steps) but with short pacing to collect multiple bodies quickly
    api.setMinStepInterval(1);
    api.startMDContinuous({ steps:50, temperature: window.__MLIP_TARGET_TEMPERATURE });
    // Wait for a few requests to accumulate
    await new Promise(r=>setTimeout(r, 50));
    const initialBodies = mdBodies.slice();
    expect(initialBodies.length).toBeGreaterThan(1);
    const firstTemps = initialBodies.map(b=>b.temperature);
    const uniqueFirst = new Set(firstTemps);
    // All initial temps should be identical (old behavior baseline)
    expect(uniqueFirst.size).toBe(1);
    // Change slider to max
    const slider = document.getElementById('mdTempSlider');
    slider.value = slider.max; slider.dispatchEvent(new Event('input'));
    const newTarget = window.__MLIP_TARGET_TEMPERATURE;
    expect(newTarget).not.toBe(firstTemps[0]);
    // Allow more steps
    await new Promise(r=>setTimeout(r, 80));
    const later = mdBodies.slice(initialBodies.length);
    expect(later.length).toBeGreaterThan(2);
    const laterTemps = later.map(b=>b.temperature);
    // Expect at least one body to reflect new target temperature
    expect(laterTemps.some(t=> t === newTarget)).toBe(true);
  }, 10000);
});
