import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

// Minimal bondService stub that records recompute calls
function makeBondService(state){
  const spy = (...args) => { spy.calls.push(args); state.markBondsChanged(); return undefined; };
  spy.calls = [];
  return { recomputeAndStore: spy };
}

describe('atom drag triggers bond recompute on end', () => {
  test('bondsChanged version increments after drag end with movement', () => {
    const st = createMoleculeState({
      elements:['C','H','H'],
      positions:[{x:0,y:0,z:0},{x:1,y:0,z:0},{x:-1,y:0,z:0}],
      bonds:[{i:0,j:1},{i:0,j:2}]
    });
    st.selection = { kind:'atom', data:{ index:0 } };
    const bondService = makeBondService(st);
    const manip = createManipulationService(st, { bondService });
    const bondsVersionStart = st.versions.bonds;
    const posVersionStart = st.versions.positions;
    manip.beginDrag(()=>({x:0,y:0,z:0}));
    // Simulate movement
    manip.updateDrag(()=>({x:0.5,y:0,z:0}));
    manip.endDrag();
    expect(st.versions.positions).toBeGreaterThan(posVersionStart);
    expect(st.versions.bonds).toBeGreaterThan(bondsVersionStart);
  expect(bondService.recomputeAndStore.calls.length).toBe(1);
  });
  test('no recompute if no movement', () => {
    const st = createMoleculeState({ elements:['C'], positions:[{x:0,y:0,z:0}], bonds:[] });
    st.selection = { kind:'atom', data:{ index:0 } };
    const bondService = makeBondService(st);
    const manip = createManipulationService(st, { bondService });
    const bondsVersionStart = st.versions.bonds;
    manip.beginDrag(()=>({x:0,y:0,z:0}));
    manip.endDrag();
    expect(st.versions.bonds).toBe(bondsVersionStart);
  expect(bondService.recomputeAndStore.calls.length).toBe(0);
  });
});
