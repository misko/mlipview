import { computeBondsNoState } from '../bond_render.js';
import { __count } from '../util/funcCount.js';

export function createBondService(molState) {
  __count('bondService#createBondService');
  function baseAtomArray() {
    __count('bondService#baseAtomArray');
    return molState.positions.map((p, i) => ({
      element: molState.elements[i],
      pos: [p.x, p.y, p.z],
    }));
  }
  function hasCell() {
    __count('bondService#hasCell');
    const c = molState.cell;
    return !!(c && c.enabled && c.a && c.b && c.c);
  }
  function computePeriodicBonds() {
    __count('bondService#computePeriodicBonds');
    // Non-periodic path: preserve full opacity shaping from computeBondsNoState (legacy smooth transparency)
    if (!hasCell())
      return computeBondsNoState(baseAtomArray()).map((b) => ({
        i: b.i,
        j: b.j,
        length: b.length,
        opacity: b.opacity ?? 1,
      }));
    const { a, b, c } = molState.cell;
    const BOND_DBG =
      typeof window !== 'undefined' &&
      (window.BOND_DEBUG === true || /[?&]bondDebug=1/.test(window.location?.search || ''));
    function longThreshold() {
      try {
        if (molState?.cell?.enabled) {
          const vx = a.x + b.x + c.x,
            vy = a.y + b.y + c.y,
            vz = a.z + b.z + c.z;
          const diag = Math.hypot(vx, vy, vz) || 0;
          const mul =
            typeof window !== 'undefined' && window.BOND_DEBUG_MULT
              ? Number(window.BOND_DEBUG_MULT)
              : 0.5;
          return diag * (Number.isFinite(mul) ? mul : 0.5);
        }
      } catch {}
      const fallback =
        typeof window !== 'undefined' && window.BOND_DEBUG_MINLEN
          ? Number(window.BOND_DEBUG_MINLEN)
          : 6.0;
      return Number.isFinite(fallback) ? fallback : 6.0;
    }
    const LONG_THR = longThreshold();
    const shifts = [
      { x: 0, y: 0, z: 0 },
      a,
      b,
      c,
      { x: -a.x, y: -a.y, z: -a.z },
      { x: -b.x, y: -b.y, z: -b.z },
      { x: -c.x, y: -c.y, z: -c.z },
    ];
    const aug = [];
    for (let si = 0; si < shifts.length; si++) {
      const S = shifts[si];
      for (let i = 0; i < molState.positions.length; i++) {
        const p = molState.positions[i];
        aug.push({
          element: molState.elements[i],
          pos: [p.x + S.x, p.y + S.y, p.z + S.z],
          baseIndex: i,
          shiftIndex: si,
        });
      }
    }
    const augBonds = computeBondsNoState(aug.map((a) => ({ element: a.element, pos: a.pos })));
    const av = a,
      bv = b,
      cv = c;
    const O =
      molState.cell && molState.cell.originOffset
        ? molState.cell.originOffset
        : { x: 0, y: 0, z: 0 };
    const M = [av.x, bv.x, cv.x, av.y, bv.y, cv.y, av.z, bv.z, cv.z];
    const det =
      M[0] * (M[4] * M[8] - M[5] * M[7]) -
      M[1] * (M[3] * M[8] - M[5] * M[6]) +
      M[2] * (M[3] * M[7] - M[4] * M[6]);
    let inv = null;
    if (Math.abs(det) > 1e-14) {
      const invDet = 1 / det;
      inv = [
        (M[4] * M[8] - M[5] * M[7]) * invDet,
        (M[2] * M[7] - M[1] * M[8]) * invDet,
        (M[1] * M[5] - M[2] * M[4]) * invDet,
        (M[5] * M[6] - M[3] * M[8]) * invDet,
        (M[0] * M[8] - M[2] * M[6]) * invDet,
        (M[2] * M[3] - M[0] * M[5]) * invDet,
        (M[3] * M[7] - M[4] * M[6]) * invDet,
        (M[1] * M[6] - M[0] * M[7]) * invDet,
        (M[0] * M[4] - M[1] * M[3]) * invDet,
      ];
    }
    function frac(p) {
      // Convert Cartesian to fractional with respect to origin offset O and lattice [a b c]
      const dx = p.x - O.x,
        dy = p.y - O.y,
        dz = p.z - O.z;
      if (!inv) return { u: dx, v: dy, w: dz };
      return {
        u: inv[0] * dx + inv[1] * dy + inv[2] * dz,
        v: inv[3] * dx + inv[4] * dy + inv[5] * dz,
        w: inv[6] * dx + inv[7] * dy + inv[8] * dz,
      };
    }
    function wrap(f) {
      return { u: f.u - Math.floor(f.u), v: f.v - Math.floor(f.v), w: f.w - Math.floor(f.w) };
    }
    function cart(f) {
      // Convert fractional back to Cartesian, add origin offset O
      return {
        x: O.x + av.x * f.u + bv.x * f.v + cv.x * f.w,
        y: O.y + av.y * f.u + bv.y * f.v + cv.y * f.w,
        z: O.z + av.z * f.u + bv.z * f.v + cv.z * f.w,
      };
    }
    const seen = new Set();
    const out = [];
    for (const eb of augBonds) {
      const A = aug[eb.i],
        B = aug[eb.j];
      // Skip self-bonds introduced by periodic augmentation (same base index)
      if (A.baseIndex === B.baseIndex) continue;
      if (A.shiftIndex === B.shiftIndex && A.shiftIndex !== 0) continue;
      const pA = { x: A.pos[0], y: A.pos[1], z: A.pos[2] };
      const pB = { x: B.pos[0], y: B.pos[1], z: B.pos[2] };
      // Hard distance pruning: skip if raw distance exceeds a generous covalent upper bound (~2.2 Å) to avoid spurious long bonds
      const dx0 = pA.x - pB.x,
        dy0 = pA.y - pB.y,
        dz0 = pA.z - pB.z;
      const dist0 = Math.sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
      if (dist0 > 2.2) continue;
      // Additional pruning: avoid creating cross-image C-H or H-H bonds; only allow hydrogen bonds from primary image (shiftIndex 0)
      const eA = molState.elements[A.baseIndex];
      const eB = molState.elements[B.baseIndex];
      if ((eA === 'H' || eB === 'H') && (A.shiftIndex !== 0 || B.shiftIndex !== 0)) continue;
      const fA = wrap(frac(pA));
      const fB = wrap(frac(pB));
      const cA = cart(fA);
      const cB = cart(fB);
      const dx = cA.x - cB.x,
        dy = cA.y - cB.y,
        dz = cA.z - cB.z;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      // Conservative lower-bound filter: drop unphysical ultra-short bonds that can appear
      // due to pathological wrapping across very tight synthetic or rotated cells.
      if (dist < 0.65) continue;
      const i = Math.min(A.baseIndex, B.baseIndex);
      const j = Math.max(A.baseIndex, B.baseIndex);
      let du = fB.u - fA.u,
        dv = fB.v - fA.v,
        dw = fB.w - fA.w;
      du -= Math.round(du);
      dv -= Math.round(dv);
      dw -= Math.round(dw);
      const key = i + '_' + j + ':' + du.toFixed(4) + ',' + dv.toFixed(4) + ',' + dw.toFixed(4);
      if (seen.has(key)) continue;
      seen.add(key);
      const fAraw = frac(pA),
        fBraw = frac(pB);
      const ndu = fBraw.u - fAraw.u - du,
        ndv = fBraw.v - fAraw.v - dv,
        ndw = fBraw.w - fAraw.w - dw;
      // Mark as crossing if either endpoint comes from a non-primary image (shiftIndex != 0),
      // OR if raw fractional delta differs from minimum-image delta (numerical crossing heuristic).
      const crossImage = A.shiftIndex !== 0 || B.shiftIndex !== 0;
      const crossing =
        crossImage || Math.abs(ndu) > 1e-6 || Math.abs(ndv) > 1e-6 || Math.abs(ndw) > 1e-6;
      // If a bond crosses the periodic boundary, hide the primary bond (opacity=0)
      // and rely on ghost bonds for visualization. This prevents a long opaque cylinder across the box.
      if (BOND_DBG && dist > LONG_THR) {
        try {
          const elI = molState.elements[i];
          const elJ = molState.elements[j];
          console.log('[BOND-DBG-LONG][compute]', {
            i,
            j,
            elements: [elI, elJ],
            length: Number(dist.toFixed(4)),
            crossing,
            atomA: {
              x: Number(pA.x.toFixed(4)),
              y: Number(pA.y.toFixed(4)),
              z: Number(pA.z.toFixed(4)),
            },
            atomB: {
              x: Number(pB.x.toFixed(4)),
              y: Number(pB.y.toFixed(4)),
              z: Number(pB.z.toFixed(4)),
            },
            cell: molState.cell
              ? {
                  a: molState.cell.a,
                  b: molState.cell.b,
                  c: molState.cell.c,
                  originOffset: molState.cell.originOffset,
                }
              : null,
          });
        } catch {}
      }
      out.push({ i, j, length: dist, opacity: crossing ? 0.0 : 1.0, crossing });
    }
    // Fallback: if periodic expansion produced no bonds, fall back to non-periodic with shaped opacity
    if (!out.length) {
      return computeBondsNoState(baseAtomArray()).map((b) => ({
        i: b.i,
        j: b.j,
        length: b.length,
        opacity: b.opacity ?? 1,
      }));
    }
    return out;
  }
  function recomputeAndStore() {
    __count('bondService#recomputeAndStore');
    const bonds = computePeriodicBonds();
    // Store opacity for renderer so it can apply group alpha (harmless to existing logic using only i,j)
    // Preserve all bonds with per-instance opacity; renderer will hide crossing ones visually.
    const DBG =
      typeof window !== 'undefined' &&
      (window.__MLIP_DEBUG_STRETCH === true || /[?&]bondStretchDebug=1/.test(window.location?.search || ''));
    if (DBG) {
      try {
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        let softish = 0;
        const samples = [];
        for (const b of bonds) {
          const op = typeof b.opacity === 'number' ? b.opacity : 1;
          if (op < min) min = op;
          if (op > max) max = op;
          if (op < 0.99) softish++;
        }
        for (let i = 0; i < Math.min(6, bonds.length); i++) {
          const b = bonds[i];
          samples.push({
            i: b.i,
            j: b.j,
            opacity: typeof b.opacity === 'number' ? Number(b.opacity.toFixed(4)) : null,
            crossing: !!b.crossing,
          });
        }
        console.log('[BondService][recompute]', {
          count: bonds.length,
          minOpacity: Number.isFinite(min) ? Number(min.toFixed(4)) : null,
          maxOpacity: Number.isFinite(max) ? Number(max.toFixed(4)) : null,
          belowThreshold: softish,
          sample: samples,
        });
      } catch (err) {
        console.warn('[BondService][recompute] debug log error', err);
      }
    }
    molState.bonds = bonds.map((b) => ({ i: b.i, j: b.j, opacity: b.opacity }));
    molState.markBondsChanged();
    return bonds;
  }
  return { computePeriodicBonds, recomputeAndStore };
}
