/** @jest-environment jsdom */
// Verify: During running MD, while dragging, positions for the dragged atom are excluded; after endDrag,
// subsequent MD responses should resume updating that atom.

const { JSDOM } = require('jsdom');

beforeAll(()=>{
  global.BABYLON = global.BABYLON || {
    Engine: function(){ this.runRenderLoop=()=>{}; },
    Scene: function(){ this.onPointerObservable={ add:()=>{} }; this.render=()=>{}; },
    Color3: function(){}, MeshBuilder:{ CreateSphere: ()=>({}) }, StandardMaterial: function(){}, ArcRotateCamera: function(){ this.attachControl=()=>{}; }
  };
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ render:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} }) }));

function mkResp(positions){ return { positions, forces: positions.map(()=>[0,0,0]), final_energy: -1.0, velocities: positions.map(()=>[0,0,0]), temperature: 1000 }; }

const queue = [];

global.fetch = async (url, opts={})=>{
  if(typeof url === 'string' && /\/serve\/md$/.test(url)){
    const resp = queue.length ? queue.shift() : mkResp([[0,0,0],[1,0,0],[0,1,0]]);
    return { ok:true, status:200, json: async ()=> resp, text: async ()=> JSON.stringify(resp) };
  }
  if(typeof url === 'string' && /\/serve\/simple$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ results:{ energy:-1.2, forces:[[0,0,0],[0,0,0],[0,0,0]] } }) };
  }
  if(typeof url === 'string' && /\/serve\/health$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ ok:true }) };
  }
  if(typeof url === 'string' && /\/serve\/relax$/.test(url)){
    return { ok:true, status:200, json: async ()=> ({ positions:[[0,0,0],[1,0,0],[0,1,0]], forces:[[0,0,0],[0,0,0],[0,0,0]], final_energy:-1.1 }) };
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
  window.viewerApi = api; initTemperatureSlider({ hudEl: hud, getViewer: ()=> api });
  return api;
}

function almostEqual(a,b,eps=1e-9){ return Math.abs(a-b) <= eps; }

describe('drag end resumes sim application', () => {
  test('after endDrag, MD updates start applying to previously dragged atom', async () => {
    const api = await setup();
    // Begin drag on atom 1
    api.selection.clickAtom(1);
    const intersector = ()=>({ x:1, y:0, z:0 });
    expect(api.manipulation.beginDrag(intersector)).toBe(true);

    // Queue an MD response that shifts all atoms by +10 (should not apply to dragged atom 1)
    const before = api.state.positions.map(p=>[p.x,p.y,p.z]);
    queue.push(mkResp(before.map(p=> [p[0]+10, p[1]+10, p[2]+10])));
    await api.mdStep({ temperature: window.__MLIP_TARGET_TEMPERATURE });
    const afterDragStep = api.state.positions.map(p=>[p.x,p.y,p.z]);
    // Atom 1 unchanged; others moved
    expect(almostEqual(afterDragStep[1][0], before[1][0])).toBe(true);
    expect(almostEqual(afterDragStep[0][0], before[0][0]+10)).toBe(true);

    // End the drag
    api.manipulation.endDrag();

    // Next MD response: shift all atoms by +1; now atom 1 should move as well
    queue.push(mkResp(afterDragStep.map(p=> [p[0]+1, p[1]+1, p[2]+1])));
    await api.mdStep({ temperature: window.__MLIP_TARGET_TEMPERATURE });
    const afterEndStep = api.state.positions.map(p=>[p.x,p.y,p.z]);
    expect(almostEqual(afterEndStep[1][0], afterDragStep[1][0]+1)).toBe(true);
    expect(almostEqual(afterEndStep[0][0], afterDragStep[0][0]+1)).toBe(true);
  });
});
