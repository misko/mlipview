// public/selection.js
// Shared picking helpers for atoms and bonds across desktop and VR.

export function makePicker(scene, molecule) {
  if (!scene || !molecule) throw new Error("makePicker requires scene and molecule");

  // Build a fast lookup from (master mesh, per-element index) -> global atom index
  const atomIndexByMesh = new Map(); // Mesh -> Map(perIndex -> globalIndex)
  if (Array.isArray(molecule.atoms)) {
    for (let g = 0; g < molecule.atoms.length; g++) {
      const a = molecule.atoms[g];
      if (!a || !a.mesh) continue;
      let perMap = atomIndexByMesh.get(a.mesh);
      if (!perMap) {
        perMap = new Map();
        atomIndexByMesh.set(a.mesh, perMap);
      }
      perMap.set(a.index, g);
    }
  }

  const bondGroups = molecule.bondGroups; // Map(pairKey -> { indices: [{i,j}], master, ... })

  const isBondMesh = (m) => !!m?.name && m.name.startsWith("bond_") && !m.name.startsWith("ghost_bond_");
  const isAtomMesh = (m) => !!m?.name && m.name.startsWith("base_") && !m.name.startsWith("ghost_atom_");

  function mapPickedBondFromPick(pick) {
    if (!(pick?.hit) || !isBondMesh(pick.pickedMesh)) return null;
    const mesh = pick.pickedMesh;
    const idx = pick.thinInstanceIndex;
    if (idx == null || idx < 0) return null;
    const key = mesh.name.replace(/^bond_/, "");
    const group = bondGroups?.get ? bondGroups.get(key) : null;
    const pair = group?.indices ? group.indices[idx] : null;
    if (!pair) return null;
    return { key, index: idx, i: pair.i, j: pair.j, mesh };
  }

  function mapPickedAtomFromPick(pick) {
    if (!(pick?.hit) || !isAtomMesh(pick.pickedMesh)) return null;
    const mesh = pick.pickedMesh;
    const perIndex = pick.thinInstanceIndex;
    if (perIndex == null || perIndex < 0) return null;
    const perMap = atomIndexByMesh.get(mesh);
    const globalIndex = perMap ? perMap.get(perIndex) : undefined;
    if (globalIndex == null) return null;
    const atom = molecule.atoms[globalIndex];
    return { idx: globalIndex, type: atom.type, perIndex, mesh, atom };
  }

  function pickBondAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY, isBondMesh);
    return mapPickedBondFromPick(pick);
  }
  function pickAtomAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY, isAtomMesh);
    return mapPickedAtomFromPick(pick);
  }

  function pickBondWithRay(ray) {
    const pick = scene.pickWithRay(ray, isBondMesh);
    return mapPickedBondFromPick(pick);
  }
  function pickAtomWithRay(ray) {
    const pick = scene.pickWithRay(ray, isAtomMesh);
    return mapPickedAtomFromPick(pick);
  }

  return {
    isBondMesh,
    isAtomMesh,
    pickBondAtPointer,
    pickAtomAtPointer,
    pickBondWithRay,
    pickAtomWithRay,
    mapPickedBondFromPick,
    mapPickedAtomFromPick
  };
}
