import { createMoleculeState } from '../../public/domain/state/moleculeState.js';
import { createBondService, computePeriodicBonds } from '../../public/domain/services/bondService.js';

// Minimal BABYLON vector mock (tests rely only on x,y,z fields)
function V(x,y,z){ return { x,y,z }; }

describe('BondService periodic minimal image', () => {
  test('computes bonds without cell gracefully', () => {
    const st = createMoleculeState({ atoms:[{ id:0, element:'C', pos:V(0,0,0) }, { id:1, element:'C', pos:V(0,0,1.2) }] });
    const svc = createBondService(st);
    const bonds = svc.recompute({ periodic:false });
    expect(Array.isArray(bonds)).toBe(true);
    expect(bonds.length).toBeGreaterThanOrEqual(1);
  });

  test('image expansion yields bond across boundary', () => {
    const cell = { a:V(4,0,0), b:V(0,4,0), c:V(0,0,4) };
    const pA = V(0,0,-1.9); // near -c/2 edge
    const pB = V(0,0, 1.9); // near +c/2 edge
    const directDist = Math.hypot(pB.x-pA.x, pB.y-pA.y, pB.z-pA.z); // ~3.8
    expect(directDist).toBeGreaterThan(2.5);
    const st = createMoleculeState({ atoms:[{ id:0, element:'C', pos:pA }, { id:1, element:'C', pos:pB }], cell });
    const bonds = computePeriodicBonds({ atoms: st.atoms, cell: st.cell, periodic:true });
    const cc = bonds.find(b=> (b.i===0 && b.j===1) || (b.i===1 && b.j===0));
    expect(cc).toBeDefined();
    expect(cc.length).toBeLessThan(0.6); // minimal image length
  });
});
