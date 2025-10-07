/** @jest-environment jsdom */
// Verifies that the instantaneous temperature returned by /serve/md updates the HUD element #instTemp.

beforeAll(()=>{
  if(!global.BABYLON){
    global.BABYLON = { TransformNode: class {}, MeshBuilder:{ CreateCylinder:()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial: class {}, Color3: class {}, Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }, Quaternion: class {}, Scene: class {} };
  }
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

// Capture last temperature sent and simulate backend returning a slightly different instantaneous temperature
let lastSentTemp = null;

global.fetch = async (url, opts={}) => {
  if(/\/serve\/md$/.test(url)){
    let reqBody={}; try { reqBody = JSON.parse(opts.body||'{}'); } catch{}
    lastSentTemp = reqBody.temperature;
    const inst = (reqBody.temperature||0) + 12.34; // backend instantaneous temp
    const body = JSON.stringify({ positions:[[0,0,0]], velocities:[[0,0,0]], forces:[[0,0,0]], final_energy:-1.0, steps_completed:1, temperature: inst });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/simple$/.test(url)){
    const body = JSON.stringify({ results:{ energy:-5, forces:[[0,0,0]], stress:[0,0,0,0,0,0] } });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/relax$/.test(url)){
    const body = JSON.stringify({ initial_energy:-5, final_energy:-5, positions:[[0,0,0]], forces:[[0,0,0]], steps_completed:1 });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  if(/\/serve\/health$/.test(url)){
    const body = JSON.stringify({ status:'ok' });
    return { ok:true, status:200, json: async ()=> JSON.parse(body), text: async ()=> body };
  }
  throw new Error('Unexpected fetch '+url);
};

async function setup(){
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"><button id="btnMD"></button><button id="btnMDRun"></button><select id="forceProviderSel"></select><button id="btnCell"></button><button id="btnGhosts"></button><button id="btnToggleForces"></button><span id="status"></span><span id="instTemp">T: --.- K</span><select id="xrModeSelect"></select></div><div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>`;
  const mod = await import('../public/index.js');
  const viewer = await mod.initNewViewer(document.getElementById('viewer'), { elements:[{Z:8}], positions:[{x:0,y:0,z:0}], bonds:[] });
  return viewer;
}

describe('Instantaneous MD temperature HUD', () => {
  test('updates #instTemp after mdStep', async () => {
    const viewer = await setup();
    // Pre condition
    const el = document.getElementById('instTemp');
    expect(el).toBeTruthy();
    expect(el.textContent).toMatch(/--/);
    // Perform an mdStep
    await viewer.mdStep({ temperature: 500 });
    expect(lastSentTemp).toBe(500);
    // HUD should reflect instantaneous value (500 + 12.34)
    expect(el.textContent).toMatch(/T: 512\.3/);
  });
});
