import { getElementInfo } from './data/periodicTable.js';
// Covalent radii pulled from centralized table with reasonable fallbacks
function covRad(sym) {
  const info = getElementInfo(sym);
  return info && info.covRad ? info.covRad : 0.8;
}
const R_SCALE = 1.15;
const SWITCH_IN = 0.9;
const SWITCH_OUT = 1.15;
const MAX_NEIGHBORS = 12;
const RING_MIN = 5,
  RING_MAX = 7;
const RING_W_MIN = 0.45;
const RING_BOOST = 0.35;
const RING_LEN_PENALTY = (L) => (L >= 5 && L <= 7 ? 1.0 : 0.0);
const OPACITY_GAMMA = 0.9;
const OPACITY_FLOOR = 0.03;
const smoothstep = (a, b, x) => {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};
const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const vlen = (v) => Math.hypot(v[0], v[1], v[2]);
const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scale = (v, s) => [v[0] * s, v[1] * s, v[2] * s];
const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
const norm = (v) => {
  const L = vlen(v) || 1;
  return scale(v, 1 / L);
};
const keyIJ = (i, j) => (i < j ? `${i}_${j}` : `${j}_${i}`);
function bondWeight(d, rc) {
  const r0 = SWITCH_IN * rc;
  const r1 = SWITCH_OUT * rc;
  return 1 - smoothstep(r0, r1, d);
}
export function computeBondsNoState(atoms) {
  const n = atoms.length;
  const positions = atoms.map((a) => a.pos);
  const elements = atoms.map((a) => a.element);
  const nbrs = Array.from({ length: n }, () => []);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const rc = R_SCALE * ((covRad(elements[i]) ?? 0.8) + (covRad(elements[j]) ?? 0.8));
      const d = vlen(vsub(positions[i], positions[j]));
      if (d < rc * SWITCH_OUT * 1.1) {
        nbrs[i].push(j);
        nbrs[j].push(i);
      }
    }
  }
  for (let i = 0; i < n; i++) {
    if (nbrs[i].length > MAX_NEIGHBORS) {
      nbrs[i].sort((a, b) => {
        const da = vlen(vsub(positions[i], positions[a]));
        const db = vlen(vsub(positions[i], positions[b]));
        return da - db;
      });
      nbrs[i] = nbrs[i].slice(0, MAX_NEIGHBORS);
    }
  }
  const edges = [];
  const wByKey = {};
  for (let i = 0; i < n; i++) {
    for (const j of nbrs[i]) {
      if (j <= i) continue;
      const rc = R_SCALE * ((covRad(elements[i]) ?? 0.8) + (covRad(elements[j]) ?? 0.8));
      const d = vlen(vsub(positions[i], positions[j]));
      const w = bondWeight(d, rc);
      if (w > 0) {
        edges.push({ i, j, d, w });
        wByKey[keyIJ(i, j)] = w;
      }
    }
  }
  const ccAdj = new Map();
  for (let i = 0; i < n; i++) if (elements[i] === 'C') ccAdj.set(i, []);
  for (const e of edges) {
    const isCC = elements[e.i] === 'C' && elements[e.j] === 'C';
    if (isCC && e.w >= RING_W_MIN) {
      ccAdj.get(e.i).push({ nb: e.j, w: e.w });
      ccAdj.get(e.j).push({ nb: e.i, w: e.w });
    }
  }
  const ringScoreByEdge = {};
  const seenCycleCanon = new Set();
  function planarity(idx) {
    if (idx.length < 3) return 0;
    const pts = idx.map((k) => positions[k]);
    let c = [0, 0, 0];
    for (const p of pts) c = vadd(c, p);
    c = scale(c, 1 / pts.length);
    let sumCross = [0, 0, 0];
    let denom = 0;
    for (let k = 0; k < pts.length; k++) {
      const a = vsub(pts[k], c);
      const b = vsub(pts[(k + 1) % pts.length], c);
      const cr = cross(a, b);
      sumCross = vadd(sumCross, cr);
      denom += vlen(cr);
    }
    if (denom < 1e-9) return 0;
    return Math.min(1, vlen(sumCross) / denom);
  }
  function dfsCycle(start, cur, path, visited, parent) {
    if (path.length > RING_MAX) return;
    for (const { nb, w } of ccAdj.get(cur) || []) {
      if (nb === parent) continue;
      const idxInPath = path.indexOf(nb);
      if (idxInPath !== -1) {
        const cyc = path.slice(idxInPath);
        if (cyc.length >= RING_MIN && cyc.length <= RING_MAX) {
          const canon = [...cyc]
            .slice()
            .sort((a, b) => a - b)
            .join(',');
          if (!seenCycleCanon.has(canon)) {
            seenCycleCanon.add(canon);
            let wsum = 0;
            for (let k = 0; k < cyc.length; k++) {
              const a = cyc[k],
                b = cyc[(k + 1) % cyc.length];
              const kw = keyIJ(a, b);
              wsum += wByKey[kw] ?? 0;
            }
            const wmean = wsum / cyc.length;
            const pScore = planarity(cyc);
            const cycScore = wmean * pScore * RING_LEN_PENALTY(cyc.length);
            for (let k = 0; k < cyc.length; k++) {
              const a = cyc[k],
                b = cyc[(k + 1) % cyc.length];
              const ke = keyIJ(a, b);
              ringScoreByEdge[ke] = Math.max(ringScoreByEdge[ke] ?? 0, cycScore);
            }
          }
        }
        continue;
      }
      if (visited.has(nb)) continue;
      if (nb < start) continue;
      visited.add(nb);
      dfsCycle(start, nb, [...path, nb], visited, cur);
      visited.delete(nb);
    }
  }
  const carbons = [...ccAdj.keys()].sort((a, b) => a - b);
  for (const c of carbons) {
    const vis = new Set([c]);
    dfsCycle(c, c, [c], vis, -1);
  }
  const out = [];
  for (const e of edges) {
    const base = Math.max(0, Math.min(1, e.w));
    const k = keyIJ(e.i, e.j);
    const rScore = ringScoreByEdge[k] ?? 0;
    const boosted = base * (1 + RING_BOOST * rScore);
    const opacity = Math.max(OPACITY_FLOOR, Math.min(1, Math.pow(boosted, OPACITY_GAMMA)));
    out.push({ i: e.i, j: e.j, length: e.d, weight: base, opacity, inRing: rScore > 0.05 });
  }
  return out;
}
