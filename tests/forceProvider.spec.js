import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createLocalLJProvider, createFairChemProvider } from '../public/physics/force-provider.js';

describe('force provider abstraction', () => {
  test('local provider returns forces and stress', async () => {
    const st = createMoleculeState({ elements:['C','C'], positions:[{x:0,y:0,z:0},{x:1.5,y:0,z:0}], bonds:[{i:0,j:1}] });
    createBondService(st);
    const provider = createLocalLJProvider(st, {});
    const res = await provider.compute({ wantStress:true });
    expect(Array.isArray(res.forces)).toBe(true);
    expect(res.forces.length).toBe(2);
    expect(res.forces[0].length).toBe(3);
    expect(res.stress).toBeTruthy();
    expect(Number.isFinite(res.stress.voigt[0])).toBe(true);
  });

  test('fairchem provider constructs correct request (mocked)', async () => {
    const fetchCalls = [];
    global.fetch = async (url, init) => {
      fetchCalls.push({ url, init });
      return { ok:true, json: async ()=> ({ results: { energy: 1.23, forces: [[0,0,0]], stress: [1,2,3,4,5,6] } }) };
    };
    const provider = createFairChemProvider({ baseUrl:'http://localhost:9999' });
    const res = await provider.compute({ elements:['C'], positions:[[0,0,0]], wantStress:true });
    expect(res.energy).toBeCloseTo(1.23);
    expect(res.forces.length).toBe(1);
    expect(res.stress.voigt[3]).toBe(4);
    expect(fetchCalls.length).toBe(1);
    const body = JSON.parse(fetchCalls[0].init.body);
    expect(body.properties.includes('stress')).toBe(true);
  });

  test('fairchem provider omits stress when not requested', async () => {
    global.fetch = async (url, init) => ({ ok:true, json: async ()=> ({ results: { energy: 0.5, forces: [[0,0,0]] } }) });
    const provider = createFairChemProvider({ baseUrl:'http://localhost:9999' });
    await provider.compute({ elements:['C'], positions:[[0,0,0]], wantStress:false });
    // parse last body from fetch (only call)
    // global.fetch mock above only sees one call per test
    // To get the body, wrap again
    let capturedBody;
    global.fetch = async (url, init) => { capturedBody = JSON.parse(init.body); return { ok:true, json: async ()=> ({ results:{ energy:0.5, forces:[[0,0,0]] } }) }; };
    await provider.compute({ elements:['C'], positions:[[0,0,0]], wantStress:false });
    expect(capturedBody.properties.includes('stress')).toBe(false);
  });
});
