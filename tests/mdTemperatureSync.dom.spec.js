/** @jest-environment jsdom */
// Test: Temperature slider and MD parameter stay in sync on initial load and after updates.

const mdBodies = [];

beforeAll(()=>{
  // Minimal BABYLON stubs
  global.BABYLON = global.BABYLON || { Engine: function(){ this.runRenderLoop=()=>{}; }, Scene: function(){ this.onPointerObservable={}; this.render=()=>{}; }, Color3: function(){}, MeshBuilder:{ CreateSphere: ()=>({}) }, StandardMaterial: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; } };
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ render:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} }) }));

// Capture /serve/md payloads
global.fetch = async (url, opts={})=>{
  if(typeof url === 'string' && /\/serve\/md(_from_cache)?$/.test(url)){
    const body = JSON.parse(opts.body||'{}');
    mdBodies.push(body);
    const resp = { positions:[[0,0,0],[1,0,0],[0,1,0]], velocities:[[0,0,0],[0,0,0],[0,0,0]], forces:[[0,0,0],[0,0,0],[0,0,0]], final_energy:-1.0, temperature: body.temperature||0 };
    return { ok:true, status:200, json: async ()=> resp, text: async ()=> JSON.stringify(resp) };
  }
  if(typeof url === 'string' && /\/serve\/simple$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ results:{ energy:-1.2, forces:[[0,0,0]] }, cache_key:'x'}) };
  }
  if(typeof url === 'string' && /\/serve\/relax$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ positions:[[0,0,0]], forces:[[0,0,0]], final_energy:-1.1, cache_key:'y' }) };
  }
  if(typeof url === 'string' && /\/serve\/health$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ ok:true }) };
  }
  throw new Error('Unexpected fetch '+url);
};

async function setup(){
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"></div>`;
  const { initNewViewer } = await import('../public/index.js');
  const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
  const canvas = document.getElementById('viewer');
  const hud = document.querySelector('.hud');
  const api = await initNewViewer(canvas, { elements:['O','H','H'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0},{x:0,y:1,z:0}], bonds:[] });
  window.viewerApi = api;
  initTemperatureSlider({ hudEl: hud, getViewer: ()=> api });
  return api;
}

describe('Temperature slider sync', () => {
  test('initial load uses 1500K and mdStep payload matches', async () => {
    mdBodies.length = 0;
    const api = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    expect(slider).toBeTruthy();
    expect(label.textContent).toContain('1500');
  // Trigger an MD step at target temperature to verify payload
  await api.mdStep({ temperature: window.__MLIP_TARGET_TEMPERATURE });
    expect(mdBodies.length).toBe(1);
    expect(mdBodies[0].temperature).toBe(window.__MLIP_TARGET_TEMPERATURE);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(1500);
  });

  test('moving slider updates __MLIP_TARGET_TEMPERATURE and mdStep payload follows', async () => {
    mdBodies.length = 0;
    const api = await setup();
    const slider = document.getElementById('mdTempSlider');
    const label = document.getElementById('tempLabel');
    // Move slider to max (near 3000)
    slider.value = slider.max;
    slider.dispatchEvent(new Event('input'));
    const target = window.__MLIP_TARGET_TEMPERATURE;
    expect(label.textContent).toContain(String(target));
  await api.mdStep({ temperature: target });
  expect(mdBodies.length).toBe(1);
    expect(mdBodies[0].temperature).toBe(target);
  });
});
