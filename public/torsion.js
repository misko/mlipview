// public/torsion.js
import { buildAdjacency, sideAtoms } from "./graph.js";

/**
 * Create a torsion controller for a molecule built with buildMolecule().
 * If a stateStore is provided, rotations are recorded and can be replayed to rebuild positions.
 */
export function createTorsionController(molecule, stateStore = null) {
  const atomsAPI = molecule.atoms; // [{ mesh, type, index, pos, scale }]
  let adj = buildAdjacency(molecule.bonds, atomsAPI.length);

  function refreshGraph() {
    adj = buildAdjacency(molecule.bonds, atomsAPI.length);
    if (stateStore) stateStore.refreshAdjacency();
  }

  /**
   * Rotate around bond i-j.
   * Preferred: pass orientation (0 -> original ordering i->j, 1 -> reversed) instead of side.
   * Backward compatible: if side provided explicitly it wins. Otherwise orientation maps:
   *   orientation 0 => side 'j' (rotate atoms on j side, pivot i)
   *   orientation 1 => side 'i' (rotate atoms on i side, pivot j)
   */
  function rotateAroundBond({ i, j, side, orientation, angleDeg = 5, recompute = false }) {
    if (i === j) return;
    if (i < 0 || j < 0 || i >= atomsAPI.length || j >= atomsAPI.length) return;

    // Resolve side from orientation if not provided
    if (!side) {
      if (orientation === 0) side = 'j';
      else if (orientation === 1) side = 'i';
      else side = 'j'; // default
    }

    const pivotIdx  = (side === "j") ? i : j;
    const movedFrom = (side === "j") ? j : i;

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

    molecule.refreshBonds();

    // Mark molecule as changed for physics cache invalidation
    if (molecule.markChanged) {
      molecule.markChanged();
    }

    // Rotation history recording removed with simplified state store

    if (recompute) {
      molecule.recomputeBonds({ hysteresis: 1.03 });
      refreshGraph();
    }
  }

  return { rotateAroundBond, refreshGraph };
}
