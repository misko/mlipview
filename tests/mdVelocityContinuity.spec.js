// Verifies velocity continuity and precomputed coexistence for MD steps.
import { jest } from '@jest/globals';

const API_BASE = 'http://localhost:8000';

// Minimal mocks similar to existing tests (avoid full Babylon dependency)
if(!global.window) global.window = {}; var window = global.window;
if(!global.document) global.document = { getElementById:()=>null };
if(!global.BABYLON) global.BABYLON = {};
// Minimal Engine/Scene + camera/vector mocks used by render/scene.js
class MockEngine { constructor(){ this._loaders=[]; } runRenderLoop(fn){ /* noop */ } stopRenderLoop(){} dispose(){} }
class MockScene { constructor(){ this.onBeforeRenderObservable={ add:()=>{} }; this.onPointerObservable={ add:()=>{} }; this.meshes=[]; this._engine=new MockEngine(); } render(){} dispose(){} }
class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } }
class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } }
class ArcRotateCamera { constructor(){ } attachControl(){} }
global.BABYLON.Engine = MockEngine;
global.BABYLON.Scene = MockScene;
global.BABYLON.Vector3 = Vector3;
global.BABYLON.Color4 = Color4;
global.BABYLON.ArcRotateCamera = ArcRotateCamera;

// Provide stable elementToZ mapping indirectly via index.js import.

// Sequence storage
const requests = [];
let mdCall = 0;

// Deterministic velocities for test
function makeVel(step){
  // 3 atoms: pattern increments slightly each step
  return [
    [0.1 + 0.01*step, 0.0, 0.0],
    [0.0, 0.2 + 0.01*step, 0.0],
    [0.0, 0.0, 0.3 + 0.01*step],
  ];
}

// Force cache seeding: simple_calculate then MD; emulate existing pattern
let simpleCalled = false;

global.fetch = jest.fn(async (url, opts) => {
  if(/simple_calculate/.test(url) || /\/serve\/simple$/.test(url)){
    simpleCalled = true;
    return { ok:true, status:200, json: async ()=> ({ results:{ energy: -10.0, forces:[[0,0,0],[0,0,0],[0,0,0]], stress:[0,0,0,0,0,0] } }) };
  }
  if(/md$/.test(url) || /\/serve\/md$/.test(url)){
    mdCall++;
    const body = JSON.parse(opts.body);
    requests.push(body);
    const v = makeVel(mdCall-1);
    return { ok:true, status:200, json: async ()=> ({
      initial_energy: body.precomputed? body.precomputed.energy : -10.0,
      final_energy: -10.0 - 0.05*mdCall,
      positions: body.coordinates,
      velocities: v,
      forces: [[0.01,0,0],[0,0.01,0],[0,0,0.01]],
      steps_completed: 1,
      temperature: 300,
    }) };
  }
  throw new Error('Unexpected URL '+url);
});

// Lightweight engine/scene mocks used by index.js internals
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


describe('md velocity continuity', () => {
  test('second MD step sends velocities from first response and includes precomputed', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer({ addEventListener:()=>{}, getBoundingClientRect:()=>({left:0,top:0}) }, { elements:['O','H','H'], positions:[[0,0,0],[0.95,0,0],[-0.24,0.93,0]], bonds:[] });

    // Seed force cache
    await api.ff.computeForces({ sync:true });
    expect(simpleCalled).toBe(true);

    // First MD step (no outgoing velocities expected)
    await api.mdStep();
    expect(requests.length).toBe(1);
    expect(requests[0].velocities).toBeUndefined();
    expect(requests[0].precomputed).toBeTruthy();

    // Second MD step (should include velocities from first response)
    await api.mdStep();
    expect(requests.length).toBe(2);
    const firstResponseV = makeVel(0);
    const sentV = requests[1].velocities;
    expect(Array.isArray(sentV)).toBe(true);
    expect(sentV.length).toBe(firstResponseV.length);
    for(let i=0;i<sentV.length;i++){
      for(let k=0;k<3;k++) expect(sentV[i][k]).toBeCloseTo(firstResponseV[i][k], 6);
    }
    expect(requests[1].precomputed).toBeTruthy();

    // Confirm internal state updated to second response velocities
    const stateV = api.state.dynamics.velocities;
    const secondResponseV = makeVel(1);
    expect(stateV.length).toBe(secondResponseV.length);
    for(let i=0;i<stateV.length;i++){
      for(let k=0;k<3;k++) expect(stateV[i][k]).toBeCloseTo(secondResponseV[i][k], 6);
    }
    api.shutdown && api.shutdown();
  }, 15000);
});
