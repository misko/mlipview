// public/state.js
// Persistent molecule state: initial positions, bonds, rotation history.
// Can rebuild coordinates from rotations alone and export XYZ.

import { buildAdjacency, sideAtoms } from "./graph.js";

export function createStateStore(molecule) {
  const atomsAPI = molecule.atoms;   // [{type, index, pos, scale, mesh}]
  const n = atomsAPI.length;

  // Immutable element list
  const elements = atomsAPI.map(a => a.type);

  // Deep copy initial positions as plain vectors
  const initial = atomsAPI.map(a => a.pos.clone());

  // Bonds (indices) from the molecule builder
  const bonds = molecule.bonds.map(({ i, j }) => ({ i, j }));

  // Adjacency for side-set queries
  let adj = buildAdjacency(bonds, n);

  // Unified operation list (ordered)
  // - Rotation op: { type:'rot', i:number, j:number, angleDeg:number }
  //   Meaning: oriented bond (i -> j). i is the fixed pivot; j-side rotates.
  // - Atom translation: { type:'move', k:number, dx:number, dy:number, dz:number }
  const ops = [];

  // Back-compat view for code that reads rotations.length, etc.
  // Each entry normalized to side:'j' so pivot=i, rotate=j
  const rotations = [];

  // Monotonic counter of rotation steps (every user rotation action increments)
  let rotationSteps = 0;

  // Utility: apply a single oriented bond rotation to a positions array (Vector3[])
  // Oriented semantics: pivot = i, rotate-side = j
  function applyRotationToPositions(positions, { i, j, angleDeg }) {
    if (i === j) return;
    const pivotIdx  = i;
    const rotateIdx = j;

    const pivotPos = positions[pivotIdx].clone();
    const otherPos = positions[rotateIdx].clone();

    const axis = otherPos.subtract(pivotPos).normalize();
    if (!isFinite(axis.lengthSquared()) || axis.lengthSquared() < 1e-10) return;

    const q = BABYLON.Quaternion.RotationAxis(axis, angleDeg * Math.PI / 180);

    // Which atoms move? Everything on the chosen side of the bond (component of rotateIdx not crossing pivotIdx)
    const movers = sideAtoms(adj, pivotIdx, rotateIdx);
    for (const k of movers) {
      const local = positions[k].subtract(pivotPos);
      const rotated = local.rotateByQuaternionToRef(q, new BABYLON.Vector3());
      positions[k] = rotated.addInPlace(pivotPos);
    }
  }

  // Rebuild all current positions by replaying ops starting from 'initial'
  function rebuildPositions() {
    const positions = initial.map(v => v.clone());
    for (const op of ops) {
      if (!op) continue;
      if (op.type === 'rot') {
        applyRotationToPositions(positions, op);
      } else if (op.type === 'move') {
        const k = op.k|0;
        if (k >= 0 && k < positions.length) {
          positions[k].x += +op.dx || 0;
          positions[k].y += +op.dy || 0;
          positions[k].z += +op.dz || 0;
        }
      }
    }
    return positions;
  }

  // Commit a fresh positions array to the scene (update thin instances + bonds)
  function commitPositionsToScene(positions) {
    for (let k = 0; k < n; k++) {
      const a = atomsAPI[k];
      a.pos.copyFrom(positions[k]); // keep API positions in sync
      const s = a.scale ?? 1;
      const m = BABYLON.Matrix.Compose(
        new BABYLON.Vector3(s, s, s),
        BABYLON.Quaternion.Identity(),
        positions[k]
      );
      molecule.updateAtomMatrixByElement(a.type, a.index, m);
    }
    molecule.refreshBonds();
  }

  // Public API

  function recordRotation(entry) {
    // Convert possibly (i,j, side) into oriented op: pivot=i, rotate=j
    let i = entry.i|0;
    let j = entry.j|0;
    const angleDeg = +entry.angleDeg || 0;
    let side = (entry.side === 'i' || entry.side === 'j') ? entry.side : 'j';
    if (i === j || angleDeg === 0) return;
    // Count this as a rotation step regardless of compaction later
    rotationSteps++;
    // Oriented pair: if side === 'i', swap so i becomes pivot
    if (side === 'i') {
      const tmp = i; i = j; j = tmp;
    }
    // Normalize helper to compact angles to (-180, 180]
    const norm = (a) => {
      let x = a % 360; // [-360, 360)
      if (x <= -180) x += 360; else if (x > 180) x -= 360;
      if (Math.abs(x) < 1e-9) x = 0;
      return x;
    };
    // Compound with previous if same oriented bond
    const last = ops.length > 0 ? ops[ops.length - 1] : null;
    if (last && last.type === 'rot' && last.i === i && last.j === j) {
      last.angleDeg = norm(last.angleDeg + angleDeg);
      if (last.angleDeg === 0) {
        // Remove no-op rotation
        ops.pop();
      }
    } else {
      const na = norm(angleDeg);
      if (na !== 0) ops.push({ type: 'rot', i, j, angleDeg: na });
    }
    // Maintain back-compat rotations mirror (always side:'j' semantics)
    const lastR = rotations.length > 0 ? rotations[rotations.length - 1] : null;
    if (lastR && lastR.i === i && lastR.j === j && lastR.side === 'j') {
      lastR.angleDeg = norm(lastR.angleDeg + angleDeg);
      if (lastR.angleDeg === 0) rotations.pop();
    } else {
      const na = norm(angleDeg);
      if (na !== 0) rotations.push({ i, j, side: 'j', angleDeg: na });
    }
  }

  function recordAtomTranslation({ k, dx = 0, dy = 0, dz = 0 }) {
    const kk = k|0;
    const vx = +dx || 0, vy = +dy || 0, vz = +dz || 0;
    if (kk < 0 || kk >= n) return;
    if (vx === 0 && vy === 0 && vz === 0) return;
    ops.push({ type: 'move', k: kk, dx: vx, dy: vy, dz: vz });
  }

  function clearOps() {
    ops.length = 0;
    rotations.length = 0;
  }

  function setOps(newOps) {
    ops.length = 0;
    rotations.length = 0;
    for (const op of newOps) {
      if (!op) continue;
      if (op.type === 'rot') recordRotation({ i: op.i, j: op.j, side: 'j', angleDeg: op.angleDeg });
      else if (op.type === 'move') recordAtomTranslation(op);
    }
  }

  // Back-compat setter for pure rotations
  function setRotations(newList) {
    setOps((newList || []).map(r => {
      let i = r.i|0, j = r.j|0; const ang = +r.angleDeg || 0; const s = (r.side === 'i' || r.side === 'j') ? r.side : 'j';
      if (s === 'i') { const tmp = i; i = j; j = tmp; }
      return { type: 'rot', i, j, angleDeg: ang };
    }));
  }

  function recomputeAndCommit() {
    const pos = rebuildPositions();
    commitPositionsToScene(pos);
  }

  function exportXYZ(name = "STATE") {
    const pos = rebuildPositionsFromRotations();
    const lines = [];
    lines.push(String(n));
    lines.push(`${name}  (reconstructed from torsions)`);
    for (let k = 0; k < n; k++) {
      const e = elements[k];
      const p = pos[k];
      lines.push(`${e} ${p.x.toFixed(6)} ${p.y.toFixed(6)} ${p.z.toFixed(6)}`);
    }
    return lines.join("\n");
  }

  function getStateJSON() {
    return JSON.stringify({
      elements,
      bonds,
      ops: ops.map(o => ({ ...o })),
      // back-compat mirror
      rotations: rotations.map(r => ({ ...r }))
    }, null, 2);
  }

  function loadStateJSON(json) {
    try {
      const obj = (typeof json === "string") ? JSON.parse(json) : json;
      // (Optional) validate elements/bonds match current molecule; here we assume same model.
      if (Array.isArray(obj.ops)) setOps(obj.ops);
      else setRotations(obj.rotations || []);
      recomputeAndCommit();
    } catch (e) {
      console.warn("[state] loadStateJSON failed:", e);
    }
  }

  function refreshAdjacency() {
    // Call this if you run molecule.recomputeBonds()
    adj = buildAdjacency(molecule.bonds, n);
  }

  return {
    molecule,               // expose for desktop change detection
    elements,
    bonds,
    initial,                 // Vector3[]
    rotations,               // back-compat view (ordered, side:'j')
    ops,                     // unified operations list
    getRotationSteps: () => rotationSteps,
    getOps: () => ops,       // accessor for debugging/inspection
    debugPrint: (tag = '') => {
      const len = ops.length;
      const last = len ? ops[len - 1] : null;
      const rlen = rotations.length;
      const rlast = rlen ? rotations[rlen - 1] : null;
      try { console.log(`[state] ${tag} ops.len=${len}`, last, '| rotations.len=', rlen, rlast); } catch {}
      return { len, last, rlen, rlast };
    },
    recordRotation,
    recordAtomTranslation,
    clearOps,
    setOps,
    setRotations,
    rebuildPositions,
    recomputeAndCommit,
    exportXYZ,
    getStateJSON,
    loadStateJSON,
    refreshAdjacency
  };
}
