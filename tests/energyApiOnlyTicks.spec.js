/**
 * energyApiOnlyTicks.spec.js
 * Verifies only API energy responses add ticks: relaxStep & mdStep increment; drag/rotation do not.
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

describe('API-only energy ticks', () => {
  test('relaxStep then mdStep add ticks; drag does not', async () => {
    // Provide fetch stubs ONLY for relax & md so drag-triggered simple_calculate calls fail (no tick),
    // but relax/md succeed and return distinct energies to add ticks.
    let seq = 0;
    global.fetch = jest.fn(async (url, opts)=>{
      if(url.includes('/serve/relax')){
        seq += 1;
        return {
          ok: true,
          json: async () => ({
            initial_energy: -10.0 + (seq-1)*0.01,
            final_energy: -10.0 + seq*0.01,
            positions: [ [0,0,0], [1.4,0,0] ],
            forces: [ [0,0,0], [0,0,0] ],
            steps_completed: 1,
            calculator: 'uma'
          })
        };
      }
      if(url.includes('/serve/md')){
        seq += 1;
        return {
          ok: true,
          json: async () => ({
            positions: [ [0,0,0], [1.4,0,0] ],
            forces: [ [0,0,0], [0,0,0] ],
            final_energy: -10.0 + seq*0.01,
            steps_completed: 1,
            calculator: 'uma'
          })
        };
      }
      // Any other endpoint (e.g., /serve/simple) simulates offline network -> rejection
      return Promise.reject(new Error('offline'));
    });
    const viewer = setupDOM();
    const mod = await import('../public/index.js');
    const api = await mod.initNewViewer(viewer, { elements:[6,6], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], bonds:[{i:0,j:1}] });
  // Allow initial async force request (computeForces kicks off automatically); wait a tick
  await new Promise(r=>setTimeout(r,30));
  const start = api.debugEnergySeriesLength();
    // Drag (no tick)
    api.state.selection = { kind:'atom', data:{ index:1 } };
    api.manipulation.beginDrag(()=>({ x:1.4,y:0,z:0 }));
    api.manipulation.updateDrag(()=>({ x:1.6,y:0,z:0 }));
    api.manipulation.endDrag();
    const afterDrag = api.debugEnergySeriesLength();
    expect(afterDrag).toBe(start);
    // Relax step (adds tick)
  await api.relaxStep();
  await new Promise(r=>setTimeout(r,20));
    const afterRelax = api.debugEnergySeriesLength();
    expect(afterRelax).toBeGreaterThan(afterDrag);
    // MD step (adds tick)
  await api.mdStep();
  await new Promise(r=>setTimeout(r,20));
    const afterMd = api.debugEnergySeriesLength();
    expect(afterMd).toBeGreaterThan(afterRelax);
  });
});
