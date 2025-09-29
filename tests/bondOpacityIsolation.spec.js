import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { computeBondsNoState } from '../public/bond_render.js';

// Simple benzene coordinates (approx planar hexagon with C-H outward) for test purposes
// We'll construct a minimal benzene: 6 carbons in ring radius ~1.4, 6 hydrogens radially outward ~2.4
function benzenePositions() {
  const C_R = 1.40; // ring radius for carbons
  const H_R = 2.40; // hydrogen radius
  const positions = [];
  const elements = [];
  for (let k=0;k<6;k++) {
    const th = (Math.PI*2*k)/6;
    positions.push({ x: C_R*Math.cos(th), y: 0, z: C_R*Math.sin(th) });
    elements.push('C');
  }
  for (let k=0;k<6;k++) {
    const th = (Math.PI*2*k)/6;
    positions.push({ x: H_R*Math.cos(th), y: 0, z: H_R*Math.sin(th) });
    elements.push('H');
  }
  return { elements, positions };
}

function asAtoms(state){
  return state.positions.map((p,i)=>({ element: state.elements[i], pos:[p.x,p.y,p.z] }));
}

function opacityMap(bonds){
  const m = new Map();
  for (const b of bonds) m.set(b.i+"-"+b.j, b.opacity);
  return m;
}

describe('bond opacity isolation when displacing a single carbon', () => {
  test('moving one carbon outward should not dim remote C-H bonds', () => {
    const { elements, positions } = benzenePositions();
    const st = createMoleculeState({ elements, positions, bonds:[] });
    const bondSvc = createBondService(st);
    // Initial compute
    let bonds0 = bondSvc.recomputeAndStore();
    // Ensure all C-H and C-C reasonably opaque initially
    const initial = opacityMap(bonds0);
    for (const b of bonds0) {
      const e1 = st.elements[b.i]; const e2 = st.elements[b.j];
      if ((e1==='C' && e2==='C') || (e1==='C' && e2==='H') || (e1==='H' && e2==='C')) {
        expect(b.opacity).toBeGreaterThan(0.85);
      }
    }
    // Move carbon 0 outward along its radial direction
    const p0 = st.positions[0];
    const dir = { x:p0.x, y:0, z:p0.z };
    const mag = Math.hypot(dir.x, dir.z) || 1;
    dir.x/=mag; dir.z/=mag;
    // Displace far
    st.positions[0] = { x:p0.x + dir.x*2.5, y:0, z:p0.z + dir.z*2.5 };
    st.markPositionsChanged();
    const bonds1 = bondSvc.recomputeAndStore();
    const after = opacityMap(bonds1);
    // Define the directly affected C-H bond (carbon 0 to its hydrogen index 6+0)
    const localCHKey = '0-' + (6+0);
    // Collect all C-H keys
    const allCH = bonds1.filter(b=>{
      const e1=st.elements[b.i], e2=st.elements[b.j];
      return (e1==='C'&&e2==='H')||(e1==='H'&&e2==='C');
    }).map(b=>b.i+'-'+b.j);
    // For every other C-H bond, opacity should remain near 1 (not globally dimmed)
    for (const key of allCH) {
      if (key === localCHKey || key === localCHKey.split('-').reverse().join('-')) continue;
      const op = after.get(key) ?? after.get(key.split('-').reverse().join('-'));
      expect(op).toBeGreaterThan(0.75);
    }
  });
});
