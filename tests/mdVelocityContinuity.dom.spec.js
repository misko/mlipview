/** @jest-environment jsdom */
// jsdom variant: validates velocity continuity using DOM canvas (still mocked network)
import { jest } from '@jest/globals';

let requests = []; let mdCall=0; let simpleCalled=false;
function makeVel(step){ return [[0.1+0.01*step,0,0],[0,0.2+0.01*step,0],[0,0,0.3+0.01*step]]; }

// Provide BABYLON minimal mocks for scene creation under jsdom
if(!global.BABYLON) global.BABYLON={};
class MockEngine { constructor(){ } runRenderLoop(fn){ this._loop=fn; } stopRenderLoop(){} dispose(){} }
class MockScene { constructor(){ this.onBeforeRenderObservable={ add:()=>{} }; this.onPointerObservable={ add:()=>{} }; this.meshes=[]; } render(){} dispose(){} }
class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }
class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } }
class ArcRotateCamera { constructor(){ } attachControl(){} }
Object.assign(global.BABYLON,{ Engine:MockEngine, Scene:MockScene, Vector3, Color4, ArcRotateCamera, PointerEventTypes:{ POINTERDOWN:1 } });

// Mocks
jest.unstable_mockModule('../public/render/moleculeView.js', () => ({ createMoleculeView: () => ({ rebuildBonds: ()=>{}, _internals:{ forceGroups:new Map() } }) }));
jest.unstable_mockModule('../public/domain/moleculeState.js', () => ({ createMoleculeState: ({ elements, positions }) => ({ elements, positions: positions.map(p=>({x:p[0],y:p[1],z:p[2]})), bonds:[], bus:{ on:()=>{}, emit:()=>{} }, markPositionsChanged(){}, dynamics:{} }) }));
jest.unstable_mockModule('../public/domain/bondService.js', () => ({ createBondService: () => ({ recomputeAndStore: ()=>[] }) }));
jest.unstable_mockModule('../public/domain/selectionService.js', () => ({ createSelectionService: () => ({}) }));
jest.unstable_mockModule('../public/core/pickingService.js', () => ({ createPickingService: ()=> ({}) }));
jest.unstable_mockModule('../public/domain/manipulationService.js', () => ({ createManipulationService: ()=> ({ beginDrag:()=>{}, updateDrag:()=>{}, endDrag:()=>{}, setDragPlane:()=>{}, rotateBond:()=>{} }) }));
jest.unstable_mockModule('../public/vr/setup.js', () => ({ createVRSupport: ()=> ({ init: async ()=> ({ supported:false }) }) }));
jest.unstable_mockModule('../public/vr/vr-picker.js', () => ({ createVRPicker: ()=> ({}) }));
jest.unstable_mockModule('../public/fairchem_provider.js', () => ({ createFairChemForcefield: ()=> ({}) }));
jest.unstable_mockModule('../public/util/funcCount.js', () => ({ __count: ()=>{} }));

global.fetch = jest.fn(async (url, opts) => {
  if(/simple_calculate/.test(url) || /\/serve\/simple$/.test(url)){
    simpleCalled = true;
    return { ok:true, status:200, json: async ()=> ({ results:{ energy: -5.0, forces:[[0,0,0],[0,0,0],[0,0,0]], stress:[0,0,0,0,0,0] } }) };
  }
  if(/\/serve\/md$/.test(url) || /(^|\/)md$/.test(url)){
    mdCall++;
    const body = JSON.parse(opts.body);
    requests.push(body);
    const v = makeVel(mdCall-1);
    return { ok:true, status:200, json: async ()=> ({
      initial_energy: body.precomputed? body.precomputed.energy : -5.0,
      final_energy: -5.0 - 0.05*mdCall,
      positions: body.coordinates,
      velocities: v,
      forces: [[0.02,0,0],[0,0.02,0],[0,0,0.02]],
      steps_completed: 1,
      temperature: 300,
    }) };
  }
  throw new Error('Unexpected URL '+url);
});

describe('md velocity continuity (jsdom)', () => {
  test('velocity reuse between md steps', async () => {
    const canvas = document.createElement('canvas');
    canvas.getBoundingClientRect = ()=>({ left:0, top:0 });
    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(canvas, { elements:['O','H','H'], positions:[[0,0,0],[0.95,0,0],[-0.24,0.93,0]], bonds:[] });
    await api.ff.computeForces({ sync:true });
    expect(simpleCalled).toBe(true);
    await api.mdStep();
    await api.mdStep();
    expect(requests.length).toBeGreaterThanOrEqual(2);
    expect(requests[0].velocities).toBeUndefined();
    expect(Array.isArray(requests[1].velocities)).toBe(true);
    const firstV = requests[1].velocities;
    for(let i=0;i<firstV.length;i++) for(let k=0;k<3;k++) expect(firstV[i][k]).toBeCloseTo(makeVel(0)[i][k], 6);
    const stateV = api.state.dynamics.velocities;
    for(let i=0;i<stateV.length;i++) for(let k=0;k<3;k++) expect(stateV[i][k]).toBeCloseTo(makeVel(1)[i][k], 6);
  });
});
