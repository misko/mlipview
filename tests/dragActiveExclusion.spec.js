/*
Verify: While an atom is actively being dragged (beginDrag without endDrag), incoming relax/MD responses
must not overwrite the dragged atom's position even if userInteractionVersion did not change (very still drag).
*/

const { JSDOM } = require('jsdom');

async function initViewer() {
  const html = `<!DOCTYPE html><html><body><canvas id="c" width="100" height="100"></canvas></body></html>`;
  const dom = new JSDOM(html, { runScripts: 'outside-only', resources: 'usable', url: 'http://localhost/' });
  global.window = dom.window; global.document = dom.window.document; global.performance = { now: ()=> Date.now() };
  window.__MLIPVIEW_TEST_MODE = true;
  // Stub minimal BABYLON API used by createScene
  function Color3(r,g,b){ this.r=r; this.g=g; this.b=b; }
  Color3.prototype.clone = function(){ return new Color3(this.r,this.g,this.b); };
  Color3.prototype.scale = function(s){ return new Color3(this.r*s,this.g*s,this.b*s); };
  Color3.FromHexString = ()=> new Color3(1,1,1);
  global.BABYLON = {
    Engine: function(){ this.stopRenderLoop=()=>{}; },
    Scene: function(){ this.dispose=()=>{}; this.onPointerObservable={}; },
    Color3,
    Color4: function(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; },
    Vector3: class { constructor(x,y,z){ this.x=x; this.y=y; this.z=z; } static Zero(){ return new this(0,0,0); } },
    MeshBuilder: { CreateSphere: ()=>({ material:{}, position:{} }) },
    StandardMaterial: function(){},
    HemisphericLight: function(){},
    ArcRotateCamera: function(){ this.attachControl=()=>{}; },
  };
  // Mock heavy rendering modules BEFORE requiring index.js
  jest.mock('../public/render/scene.js', () => ({
    createScene: async () => ({ engine:{ stopRenderLoop:()=>{} }, scene:{ render:()=>{}, dispose:()=>{}, onPointerObservable:{ add:()=>{} } }, camera:{} })
  }));
  jest.mock('../public/core/pickingService.js', () => ({
    createPickingService: () => ({})
  }));
  jest.mock('../public/render/moleculeView.js', () => ({
    createMoleculeView: () => ({ rebuildBonds:()=>{}, rebuildForces:()=>{} })
  }));
  const { initNewViewer } = require('../public/index.js');
  const canvas = document.getElementById('c');
  const elements = ['O','H','H'];
  const positions = [ {x:0,y:0,z:0},{x:0.95,y:0,z:0},{x:-0.24,y:0.93,z:0} ];
  const bonds = [ { i:0, j:1, order:1, opacity:1 }, { i:0, j:2, order:1, opacity:1 } ];
  const api = await initNewViewer(canvas, { elements, positions, bonds });
  return api;
}

function delay(ms){ return new Promise(r=>setTimeout(r, ms)); }

// Basic Response polyfill for jsdom environment if not present
if (typeof Response === 'undefined') {
  global.Response = class {
    constructor(body, init){ this._body = body; this.status = init.status; this.headers = new Map(Object.entries(init.headers||{})); }
    async json(){ return JSON.parse(this._body); }
    get ok(){ return this.status>=200 && this.status<300; }
    async text(){ return this._body; }
  };
}

describe('drag exclusion from incoming sim updates', () => {
  jest.setTimeout(20000);

  test('relax response does not overwrite actively dragged atom', async () => {
    const api = await initViewer();
    // Select atom 1 and begin a drag
    api.selection.clickAtom(1);
    const intersector = ()=>({ x:0.95, y:0, z:0 });
    expect(api.manipulation.beginDrag(intersector)).toBe(true);
    const before = api.state.positions.map(p=>[p.x,p.y,p.z]);

    const origFetch = global.fetch || window.fetch;
    // Patch relax to quickly return shifted positions (which should be applied to others, but not atom 1)
    global.fetch = async (url, opts) => {
      if (typeof url === 'string' && url.includes('/relax')) {
        const shifted = before.map((p,i)=> [p[0]+1, p[1]+2, p[2]+3]);
        return new Response(JSON.stringify({ positions: shifted, forces: [], final_energy: -1 }), { status:200, headers:{'Content-Type':'application/json'} });
      }
      return origFetch(url, opts);
    };

    const res = await api.relaxStep();
    global.fetch = origFetch;

    expect(res.applied).toBe(true);
    expect(res.partial || true).toBe(true);

    const after = api.state.positions.map(p=>[p.x,p.y,p.z]);
    // Atom 1 should remain at its pre-update value (no +1/+2/+3 applied)
    expect(after[1][0]).toBeCloseTo(before[1][0], 12);
    expect(after[1][1]).toBeCloseTo(before[1][1], 12);
    expect(after[1][2]).toBeCloseTo(before[1][2], 12);
    // Other atoms receive the update
    expect(after[0][0]).toBeCloseTo(before[0][0] + 1, 12);
    expect(after[2][1]).toBeCloseTo(before[2][1] + 2, 12);

    // End the drag (cleanup)
    api.manipulation.endDrag();
  });
});
