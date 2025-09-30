/**
 * energyInteractions.spec.js
 * Verifies that atom drag and bond rotation interactions trigger energy plot updates.
 * @jest-environment jsdom
 */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn)=>{} },
    scene: { onPointerObservable:{ _l:[], add(fn){this._l.push(fn);}, notify(){}, notifyObservers(){} } },
    camera: { attachControl: ()=>{} }
  })
}));

if (!global.BABYLON) {
  global.BABYLON = {
    Vector3: class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } },
    Color4: class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } },
    ArcRotateCamera: class {}, HemisphericLight: class {}, Engine: class {}
  };
}

function setupDOM(){
  const viewer = document.createElement('canvas'); viewer.id='viewer'; document.body.appendChild(viewer);
  const energyWrapper = document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
  const energyCanvas = document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillStyle:null, strokeStyle:null}); energyWrapper.appendChild(energyCanvas);
  const energyLabel = document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);
  return viewer;
}

describe('energy plot updates on interactions', () => {
  test('atom drag triggers energy recompute & record', async () => {
    const viewer = setupDOM();
    const mod = await import('../public/index.js');
    const api = await mod.initNewViewer(viewer, { elements:[6,6], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], bonds:[{i:0,j:1}] });
    const before = api.debugEnergySeriesLength();
    // Simulate selecting atom 1
    api.selection.set?.({ kind:'atom', data:{ index:1 } });
    api.state.selection = { kind:'atom', data:{ index:1 } }; // fallback direct set
    // Begin drag (we supply a fake intersector returning plane point)
    api.manipulation.beginDrag(()=>({ x:1.4,y:0,z:0 }));
    api.manipulation.updateDrag(()=>({ x:1.6,y:0,z:0 })); // move atom
    api.manipulation.endDrag();
    const after = api.debugEnergySeriesLength();
    expect(after).toBeGreaterThan(before);
  });

  test('bond rotation triggers energy update', async () => {
    const viewer = setupDOM();
    const mod = await import('../public/index.js');
    const api = await mod.initNewViewer(viewer, { elements:[6,6,1], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0},{x:2.8,y:0,z:0}], bonds:[{i:0,j:1},{i:1,j:2}] });
    const before = api.debugEnergySeriesLength();
    // Select central bond
    api.state.selection = { kind:'bond', data:{ i:0, j:1, orientation: { axis:'iToJ', side:'j' } } };
    api.manipulation.rotateBond(Math.PI/12);
    const after = api.debugEnergySeriesLength();
    expect(after).toBeGreaterThan(before);
  });
});
