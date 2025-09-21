import { buildAdjacency, sideAtoms } from "./graph.js";

/**
 * Create a torsion controller for a molecule built with buildMolecule().
 * - Computes adjacency from molecule.bonds
 * - Exposes rotateAroundBond({ i, j, side, angleDeg, recompute=false })
 */
export function createTorsionController(molecule) {
  const atomsAPI = molecule.atoms; // [{ mesh, type, index, pos, scale }]
  let adj = buildAdjacency(molecule.bonds, atomsAPI.length);

  function refreshGraph() {
    adj = buildAdjacency(molecule.bonds, atomsAPI.length);
  }

  function rotateAroundBond({ i, j, side = "j", angleDeg = 5, recompute = false }) {
    if (i === j) return;
    if (i < 0 || j < 0 || i >= atomsAPI.length || j >= atomsAPI.length) return;

    const pivotIdx = side === "j" ? i : j;
    const movedFrom = side === "j" ? j : i;

    const rotatingSet = sideAtoms(adj, pivotIdx, movedFrom); // set of global atom indices

    const pivotPos = atomsAPI[pivotIdx].pos.clone();
    const otherPos = atomsAPI[movedFrom].pos.clone();
    const axis = otherPos.subtract(pivotPos).normalize();
    if (!isFinite(axis.lengthSquared()) || axis.lengthSquared() < 1e-8) return;

    const angleRad = angleDeg * Math.PI / 180;
    const q = BABYLON.Quaternion.RotationAxis(axis, angleRad);

    for (const k of rotatingSet) {
      const a = atomsAPI[k];
      const local = a.pos.subtract(pivotPos);
      const rotated = local.rotateByQuaternionToRef(q, new BABYLON.Vector3());
      const newPos = rotated.addInPlace(pivotPos);

      const s = a.scale ?? 1;
      const mat = BABYLON.Matrix.Compose(
        new BABYLON.Vector3(s, s, s),
        BABYLON.Quaternion.Identity(),
        newPos
      );
      molecule.updateAtomMatrixByElement(a.type, a.index, mat);
      a.pos.copyFrom(newPos);
    }

    // Bulk-refresh all bond buffers so cylinders follow immediately
    molecule.refreshBonds();

    if (recompute) {
      molecule.recomputeBonds({ hysteresis: 1.03 });
      refreshGraph();
    }
  }

  return { rotateAroundBond, refreshGraph };
}
