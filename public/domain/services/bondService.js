// BondService: wraps existing computeBondsNoState with optional periodic minimal-image extension.
import { computeBondsNoState } from '../../bond_render.js';

// Periodic bond expansion: duplicate atoms into ±a, ±b, ±c images and then map back minimal image.
// For phase 1 we only replicate logic conceptually (no full optimization). Pure function.
// Enhanced periodic logic: expand atom list with 6 neighbor image shifts then collapse to unique minimal pairs.
export function computePeriodicBonds({ atoms, cell, periodic = true }) {
  if (!periodic || !cell || !cell.a || !cell.b || !cell.c) {
    return computeBondsNoState(atoms.map(a => ({ element: a.element, pos: [a.pos.x,a.pos.y,a.pos.z] })));
  }
  const a=cell.a, b=cell.b, c=cell.c;
  const shifts = [
    { x:0,y:0,z:0 }, a, { x:-a.x,y:-a.y,z:-a.z },
    b, { x:-b.x,y:-b.y,z:-b.z },
    c, { x:-c.x,y:-c.y,z:-c.z }
  ];
  // Build expanded atom list
  const expanded = [];
  const meta = []; // parallel meta: { origIndex, shift }
  for (let si=0; si<shifts.length; si++) {
    const s = shifts[si];
    for (let i=0;i<atoms.length;i++) {
      const ap = atoms[i].pos;
      expanded.push({ element: atoms[i].element, pos: [ap.x + s.x, ap.y + s.y, ap.z + s.z] });
      meta.push({ origIndex: i, shift: s });
    }
  }
  const rawBonds = computeBondsNoState(expanded);
  // Collapse by original pair (i<j)
  const bestByKey = new Map();
  function key(i,j){ return i<j ? i+"_"+j : j+"_"+i; }
  for (const rb of rawBonds) {
    const m1 = meta[rb.i], m2 = meta[rb.j];
    const i0 = m1.origIndex, j0 = m2.origIndex;
    if (i0 === j0) continue; // skip self from image overlap
    const k = key(i0,j0);
    // Compute actual vector using shift difference + original positions
    const p1 = atoms[i0].pos; const p2 = atoms[j0].pos;
    const dx = (p2.x + m2.shift.x) - (p1.x + m1.shift.x);
    const dy = (p2.y + m2.shift.y) - (p1.y + m1.shift.y);
    const dz = (p2.z + m2.shift.z) - (p1.z + m1.shift.z);
    const L = Math.hypot(dx,dy,dz);
    const existing = bestByKey.get(k);
    if (!existing || L < existing.length - 1e-9) {
      // Store minimal image candidate adopting weight/opacity from raw bond
      bestByKey.set(k, {
        i: i0,
        j: j0,
        length: L,
        weight: rb.weight,
        opacity: rb.opacity,
        inRing: rb.inRing
      });
    }
  }
  return [...bestByKey.values()];
}

export function createBondService(state) {
  return {
    recompute({ periodic=true }={}) {
      const bonds = computePeriodicBonds({ atoms: state.atoms, cell: state.cell, periodic });
      return bonds;
    }
  };
}
