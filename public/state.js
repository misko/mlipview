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

  // Torsion rotation history (ordered list)
  // Each entry: { i, j, side:'i'|'j', angleDeg:number }
  const rotations = [];

  // Utility: apply a single bond rotation to a positions array (Vector3[])
  function applyRotationToPositions(positions, { i, j, side = "j", angleDeg }) {
    if (i === j) return;
    const pivotIdx  = (side === "j") ? i : j;
    const rotateIdx = (side === "j") ? j : i;

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

  // Rebuild all current positions by replaying rotations starting from 'initial'
  function rebuildPositionsFromRotations() {
    const positions = initial.map(v => v.clone());
    for (const r of rotations) applyRotationToPositions(positions, r);
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
    // Normalize entry & push to history
    const e = {
      i: entry.i|0,
      j: entry.j|0,
      side: (entry.side === "i") ? "i" : "j",
      angleDeg: +entry.angleDeg || 0
    };
    if (e.i === e.j || e.angleDeg === 0) return;
    rotations.push(e);
  }

  function clearRotations() {
    rotations.length = 0;
  }

  function setRotations(newList) {
    rotations.length = 0;
    for (const r of newList) {
      recordRotation(r);
    }
  }

  function recomputeAndCommit() {
    const pos = rebuildPositionsFromRotations();
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
      rotations: rotations.map(r => ({ ...r }))
    }, null, 2);
  }

  function loadStateJSON(json) {
    try {
      const obj = (typeof json === "string") ? JSON.parse(json) : json;
      // (Optional) validate elements/bonds match current molecule; here we assume same model.
      setRotations(obj.rotations || []);
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
    elements,
    bonds,
    initial,                 // Vector3[]
    rotations,               // live array (ordered)
    recordRotation,
    clearRotations,
    setRotations,
    rebuildPositionsFromRotations,
    recomputeAndCommit,
    exportXYZ,
    getStateJSON,
    loadStateJSON,
    refreshAdjacency
  };
}
