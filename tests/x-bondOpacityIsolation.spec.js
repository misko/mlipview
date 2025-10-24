import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';

function benzenePositions() {
  const C_R = 1.4;
  const H_R = 2.4;
  const positions = [];
  const elements = [];
  for (let k = 0; k < 6; k++) {
    const th = (Math.PI * 2 * k) / 6;
    positions.push({ x: C_R * Math.cos(th), y: 0, z: C_R * Math.sin(th) });
    elements.push('C');
  }
  for (let k = 0; k < 6; k++) {
    const th = (Math.PI * 2 * k) / 6;
    positions.push({ x: H_R * Math.cos(th), y: 0, z: H_R * Math.sin(th) });
    elements.push('H');
  }
  return { elements, positions };
}

function opacityMap(bonds) {
  const m = new Map();
  for (const b of bonds) m.set(`${b.i}-${b.j}`, b.opacity);
  return m;
}

describe('x-bond opacity isolation', () => {
  test('moving one carbon outward should not dim remote C-H bonds', () => {
    const { elements, positions } = benzenePositions();
    const st = createMoleculeState({ elements, positions, bonds: [] });
    const bondSvc = createBondService(st);

    const bonds0 = bondSvc.recomputeAndStore();
    for (const b of bonds0) {
      const e1 = st.elements[b.i];
      const e2 = st.elements[b.j];
      if ((e1 === 'C' && e2 === 'C') || e1 === 'C' || e2 === 'C') {
        expect(b.opacity).toBeGreaterThan(0.85);
      }
    }

    const p0 = st.positions[0];
    const dir = { x: p0.x, y: 0, z: p0.z };
    const mag = Math.hypot(dir.x, dir.z) || 1;
    dir.x /= mag;
    dir.z /= mag;
    st.positions[0] = { x: p0.x + dir.x * 2.5, y: 0, z: p0.z + dir.z * 2.5 };
    st.markPositionsChanged();

    const bonds1 = bondSvc.recomputeAndStore();
    const after = opacityMap(bonds1);

    const localKey = `0-${6 + 0}`;
    const getOpacity = (key) => after.get(key) ?? after.get(key.split('-').reverse().join('-'));
    for (const b of bonds1) {
      const e1 = st.elements[b.i];
      const e2 = st.elements[b.j];
      if (!((e1 === 'C' && e2 === 'H') || (e1 === 'H' && e2 === 'C'))) continue;
      const key = `${b.i}-${b.j}`;
      if (key === localKey || key === localKey.split('-').reverse().join('-')) continue;
      expect(getOpacity(key)).toBeGreaterThan(0.75);
    }
  });
});
