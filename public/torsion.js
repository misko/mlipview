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

  function rotateAroundBond({ i, j, side = "j", angleDeg = 5, recompute = false, record = true }) {
    if (i === j) return;
    if (i < 0 || j < 0 || i >= atomsAPI.length || j >= atomsAPI.length) return;

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

    if (record && stateStore) {
      stateStore.recordRotation({ i, j, side, angleDeg });
    }

    if (recompute) {
      molecule.recomputeBonds({ hysteresis: 1.03 });
      refreshGraph();
    }
  }

  return { rotateAroundBond, refreshGraph };
}
