/**
 * energyPlot.spec.js
 * Reproduces missing energy plot issue: ensure energy time series seeds and grows after steps.
 * @jest-environment jsdom
 */

// Mock the scene creation to avoid real BABYLON dependency
jest.mock('../public/render/scene.js', () => ({
  createScene: async (canvas) => ({
    engine: { runRenderLoop: (fn)=>{/* no render loop in test */} },
    scene: { onPointerObservable:{ _l:[], add(fn){this._l.push(fn);}, notify(){}, notifyObservers(){} } },
    camera: { attachControl: ()=>{} }
  })
}));

// Provide minimal BABYLON global for forcefield or any vector math referencing it
if (!global.BABYLON) {
  global.BABYLON = {
    Vector3: class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } length(){ return Math.hypot(this.x,this.y,this.z);} },
    Color4: class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } },
    ArcRotateCamera: class {}, HemisphericLight: class {}, Engine: class {}
  };
}

describe('energy plot integration', () => {
  let api;
  beforeAll(async () => {
    // Canvas for viewer (even though mocked scene won't use context)
    const viewer = document.createElement('canvas'); viewer.id='viewer'; document.body.appendChild(viewer);
    // Energy plot elements expected by index.js
    const energyWrapper = document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
    const energyCanvas = document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext = () => ({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, fillText(){}, strokeStyle:null, lineWidth:1 }); energyWrapper.appendChild(energyCanvas);
    const energyLabel = document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);
    // Required HUD controls (subset)
    const btnRelax = document.createElement('button'); btnRelax.id='btnRelax'; document.body.appendChild(btnRelax);
    const btnMD = document.createElement('button'); btnMD.id='btnMD'; document.body.appendChild(btnMD);
    const btnStop = document.createElement('button'); btnStop.id='btnStopSim'; document.body.appendChild(btnStop);
    const selProvider = document.createElement('select'); selProvider.id='forceProviderSel'; document.body.appendChild(selProvider);

    const mod = await import('../public/index.js');
    api = await mod.initNewViewer(viewer, { elements:[], positions:[], bonds:[] });
  });

  test('energy series seeded on init', () => {
    expect(api).toBeTruthy();
    const len = api.debugEnergySeriesLength();
    expect(len).toBeGreaterThanOrEqual(1);
  });

  test('energy series grows after relax & md steps', () => {
    const before = api.debugEnergySeriesLength();
    api.relaxStep();
    api.mdStep();
    const after = api.debugEnergySeriesLength();
    expect(after).toBeGreaterThan(before);
  });
});
