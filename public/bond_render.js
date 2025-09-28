// ---------- Types (JSDoc for clarity) ----------
/** @typedef {[number, number, number]} Vec3 */
/** @typedef {{ element: string, pos: Vec3 }} Atom */
/**
 * @typedef {Object} BondOut
 * @property {number} i
 * @property {number} j
 * @property {number} length
 * @property {number} weight
 * @property {number} opacity
 * @property {boolean} inRing
 */

// ---------- Config (tweak as needed) ----------
/** @type {Record<string, number>} */
const COV_RAD = {
  H: 0.31, C: 0.76, N: 0.71, O: 0.66, S: 1.05, F: 0.57, Cl: 1.02, P: 1.07,
};
const R_SCALE = 1.15;      // slightly expand cutoffs
const SWITCH_IN = 0.90;    // lower edge of soft window (as multiple of rc)
const SWITCH_OUT = 1.15;   // upper edge of soft window
const MAX_NEIGHBORS = 12;  // local pruning for speed

// Ring detection / opacity shaping
const RING_MIN = 5, RING_MAX = 7;
const RING_W_MIN = 0.45;   // only consider reasonably strong CC edges as ring candidates
const RING_BOOST = 0.35;   // multiplicative opacity boost factor (scaled by ring score)
const RING_LEN_PENALTY = (L) => (L >= 5 && L <= 7 ? 1.0 : 0.0); // can taper if you want
const OPACITY_GAMMA = 0.90;  // base opacity = w^gamma (keeps smooth but readable)
const OPACITY_FLOOR = 0.03;  // never fully vanish (helps anti-flicker visually)

// ---------- Small math helpers ----------
const smoothstep = (a, b, x) => {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};
const vsub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const vlen = (v) => Math.hypot(v[0], v[1], v[2]);
const vadd = (a, b) => [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
const scale = (v, s) => [v[0]*s, v[1]*s, v[2]*s];
const cross = (a, b) => [
  a[1]*b[2] - a[2]*b[1],
  a[2]*b[0] - a[0]*b[2],
  a[0]*b[1] - a[1]*b[0],
];
const norm = (v) => {
  const L = vlen(v) || 1;
  return scale(v, 1/L);
};
const keyIJ = (i, j) => (i < j ? `${i}_${j}` : `${j}_${i}`);

// Base bond strength from soft covalent cutoff
function bondWeight(d, rc) {
  const r0 = SWITCH_IN * rc;
  const r1 = SWITCH_OUT * rc;
  return 1 - smoothstep(r0, r1, d); // 1 strong -> 0 weak
}

// ---------- Main (single-frame, stateless) ----------
/**
 * @param {Atom[]} atoms
 * @returns {BondOut[]}
 */
export function computeBondsNoState(atoms) {
  const n = atoms.length;
  const positions = atoms.map(a => a.pos);
  const elements = atoms.map(a => a.element);

  // 1) Neighbor pre-pass (prune by max allowed distance)
  const nbrs = Array.from({ length: n }, () => []);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const rc = R_SCALE * ((COV_RAD[elements[i]] ?? 0.8) + (COV_RAD[elements[j]] ?? 0.8));
      const d = vlen(vsub(positions[i], positions[j]));
      if (d < rc * SWITCH_OUT * 1.1) {
        nbrs[i].push(j);
        nbrs[j].push(i);
      }
    }
  }
  // Degree pruning for speed in dense environments
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

  // 2) Build weighted bond list (continuous in geometry)
  const edges = [];
  const wByKey = {};
  for (let i = 0; i < n; i++) {
    for (const j of nbrs[i]) {
      if (j <= i) continue;
      const rc = R_SCALE * ((COV_RAD[elements[i]] ?? 0.8) + (COV_RAD[elements[j]] ?? 0.8));
      const d = vlen(vsub(positions[i], positions[j]));
      const w = bondWeight(d, rc);
      if (w > 0) {
        edges.push({ i, j, d, w });
        wByKey[keyIJ(i, j)] = w;
      }
    }
  }

  // 3) Ring score (single-frame, continuous)
  // Build adjacency on carbon subgraph with edge weight w (>= RING_W_MIN)
  const ccAdj = new Map();
  for (let i = 0; i < n; i++) if (elements[i] === "C") ccAdj.set(i, []);
  for (const e of edges) {
    const isCC = (elements[e.i] === "C" && elements[e.j] === "C");
    if (isCC && e.w >= RING_W_MIN) {
  ccAdj.get(e.i).push({ nb: e.j, w: e.w });
  ccAdj.get(e.j).push({ nb: e.i, w: e.w });
    }
  }

  // Find small cycles (length 5–7). We’ll accumulate a *continuous* score per edge:
  //   ringScore(edge) = max over cycles containing edge of [ mean_w(cycle) * planarity(cycle) ]
  // The mean_w is based on current edge weights (continuous); planarity is continuous.
  const ringScoreByEdge = {};
  const seenCycleCanon = new Set();

  // Helper: planarity from polygon of points (continuous): |sum cross| / sum |cross|
  function planarity(idx) {
    if (idx.length < 3) return 0;
    const pts = idx.map(k => positions[k]);
    // centroid
  let c = [0,0,0];
    for (const p of pts) c = vadd(c, p);
    c = scale(c, 1/pts.length);
    // crosses
  let sumCross = [0,0,0];
    let denom = 0;
    for (let k = 0; k < pts.length; k++) {
      const a = vsub(pts[k], c);
      const b = vsub(pts[(k+1)%pts.length], c);
      const cr = cross(a, b);
  sumCross = vadd(sumCross, cr);
  denom += vlen(cr);
    }
    if (denom < 1e-9) return 0;
    return Math.min(1, vlen(sumCross) / denom); // 0..1
  }

  // Bounded DFS to collect simple cycles (length 5..7), starting on carbon nodes
  function dfsCycle(start, cur, path, visited, parent) {
    if (path.length > RING_MAX) return;
    for (const { nb, w } of (ccAdj.get(cur) || [])) {
      if (nb === parent) continue;
      const idxInPath = path.indexOf(nb);
      if (idxInPath !== -1) {
        const cyc = path.slice(idxInPath);
        if (cyc.length >= RING_MIN && cyc.length <= RING_MAX) {
          const canon = [...cyc].slice().sort((a,b)=>a-b).join(",");
          if (!seenCycleCanon.has(canon)) {
            seenCycleCanon.add(canon);
            // cycle strength: mean of edge weights in cycle
            let wsum = 0;
            for (let k = 0; k < cyc.length; k++) {
              const a = cyc[k], b = cyc[(k+1)%cyc.length];
              const kw = keyIJ(a, b);
              wsum += (wByKey[kw] ?? 0);
            }
            const wmean = wsum / cyc.length;
            const pScore = planarity(cyc);
            const cycScore = wmean * pScore * RING_LEN_PENALTY(cyc.length); // 0..1
            // write score to all ring edges
            for (let k = 0; k < cyc.length; k++) {
              const a = cyc[k], b = cyc[(k+1)%cyc.length];
              const ke = keyIJ(a, b);
              ringScoreByEdge[ke] = Math.max(ringScoreByEdge[ke] ?? 0, cycScore);
            }
          }
        }
        continue;
      }
      if (visited.has(nb)) continue;
      // order pruning: only expand nodes >= start to suppress permutations
      if (nb < start) continue;
      visited.add(nb);
      dfsCycle(start, nb, [...path, nb], visited, cur);
      visited.delete(nb);
    }
  }

  const carbons = [...ccAdj.keys()].sort((a,b)=>a-b);
  for (const c of carbons) {
    const vis = new Set([c]);
    dfsCycle(c, c, [c], vis, -1);
  }

  // 4) Finalize opacity per edge (continuous, stateless)
  const out = [];
  for (const e of edges) {
    const base = Math.max(0, Math.min(1, e.w));
    const k = keyIJ(e.i, e.j);
    const rScore = ringScoreByEdge[k] ?? 0;          // 0..1
    const boosted = base * (1 + RING_BOOST * rScore);
    const opacity = Math.max(OPACITY_FLOOR, Math.min(1, Math.pow(boosted, OPACITY_GAMMA)));
    out.push({
      i: e.i,
      j: e.j,
      length: e.d,
      weight: base,
      opacity,
      inRing: rScore > 0.05
    });
  }

  return out;
}
