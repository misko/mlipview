/** @jest-environment jsdom */
// Validate that the forces button toggles label and state, and emits visibility changes

beforeAll(()=>{
  if(!global.BABYLON){
    global.BABYLON = { TransformNode: class {}, MeshBuilder: { CreateCylinder: ()=>({ dispose(){}, setEnabled(){}, thinInstanceSetBuffer(){}, thinInstanceRefreshBoundingInfo(){}, isPickable:false, visibility:1 }) }, StandardMaterial: class { constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }}, Color3: class {}, Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new global.BABYLON.Vector3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);} normalize(){ const m=Math.hypot(this.x,this.y,this.z)||1; this.x/=m; this.y/=m; this.z/=m; return this;} rotateByQuaternionToRef(_, out){ out.x=this.x; out.y=this.y; out.z=this.z; return out; } }, Quaternion: class { static Identity(){ return {}; } static RotationAxis(){ return {}; } }, Scene: class {}, Matrix: class { static Compose(){ return {}; } } };
  }
});

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

// Minimal fetch stubs for endpoints invoked during viewer init
global.fetch = async function(url, opts={}){
  if(/\/serve\/health$/.test(String(url))) return { ok:true, status:200, json: async ()=> ({ ok:true }) };
  if(/\/serve\/(simple|relax|md)$/.test(String(url))) return { ok:true, status:200, json: async ()=> ({ energy:-1, positions:[[0,0,0]], forces:[[0,0,0]] }) };
  return { ok:true, status:200, json: async ()=> ({}) };
};

async function setup(){
  document.body.innerHTML = `<canvas id="viewer"></canvas><div class="hud"><button id="btnRelax"></button><button id="btnRelaxRun"></button><button id="btnMD"></button><button id="btnMDRun"></button><select id="forceProviderSel"></select><button id="btnCell"></button><button id="btnGhosts"></button><button id="btnToggleForces"></button><span id="status">Ready</span><select id="xrModeSelect"></select></div><div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>`;
  const { initNewViewer } = await import('../public/index.js');
  const { initForcesToggle } = await import('../public/ui/forcesToggle.js');
  const canvas = document.getElementById('viewer');
  const viewer = await initNewViewer(canvas, { elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.96,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] });
  window.viewerApi = viewer;
  initForcesToggle({ getViewer: ()=> viewer });
  return viewer;
}

describe('forces toggle button', () => {
  test('label toggles between on/off and state flips', async () => {
    const viewer = await setup();
    const btn = document.getElementById('btnToggleForces');
    expect(btn).toBeTruthy();
  // Initial default is OFF -> label should invite turning them on
  expect(btn.textContent).toMatch(/forces on/i);
  const initial = viewer.state.showForces;
  expect(initial).toBe(false);
  // Click to turn ON
  btn.click();
  expect(viewer.state.showForces).toBe(true);
  expect(btn.textContent).toMatch(/forces off/i);
  // Click to turn OFF
  btn.click();
  expect(viewer.state.showForces).toBe(false);
  expect(btn.textContent).toMatch(/forces on/i);
  });
});
