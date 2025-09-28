import { elInfo } from "./elements.js";
import { quatYto, setThinInstanceMatrices } from "./helpers.js";

/**
 * Generic molecule builder using thin instances.
 * atoms: [{ element, pos: Vector3 }]
 * bonds?: [{ i, j }]  (optional; auto-detect if omitted)
 * bondScale: number   (cutoff multiplier for auto bonds; default 1.15)
 * mode: 'ballstick' | 'spacefill'
 * bondStyles: per element pair: { "C-C": { radius, color }, "C-H": { ... } }
 */
export function buildMolecule(scene, opts) {
  const {
    atoms = [],
    bonds = null,
    bondScale = 1.15,
    mode = "ballstick",
    debugAlwaysActive = true,
    bondStyles = {}
  } = opts || {};

  // --- group atoms by element; masters are UNIT spheres, per-instance scale carries size
  const groups = new Map(); // elem -> { master, mats: [], indices:[], mat, info }
  const perAtomScale = new Array(atoms.length);

  for (let i = 0; i < atoms.length; i++) {
    const { element, pos } = atoms[i];
    const info = elInfo(element);

    if (!groups.has(element)) {
      const m = new BABYLON.StandardMaterial(`mat_${element}`, scene);
      m.diffuseColor = info.color.clone();
      m.emissiveColor = info.color.scale(0.06);

      const base = BABYLON.MeshBuilder.CreateSphere(`base_${element}`, {
        diameter: 1, segments: 24
      }, scene);
      base.material = m;
      base.isPickable = true;
      base.thinInstanceEnablePicking = true;
      if (debugAlwaysActive) base.alwaysSelectAsActiveMesh = true;

      groups.set(element, { master: base, mats: [], indices: [], mat: m, info });
    }

    const g = groups.get(element);
    const d = (mode === "spacefill") ? info.vdw * 0.7 : info.scale; // visual diameter
    perAtomScale[i] = d;

    const mat = BABYLON.Matrix.Compose(
      new BABYLON.Vector3(d, d, d),
      BABYLON.Quaternion.Identity(),
      pos
    );
    g.mats.push(mat);
    g.indices.push(i);
  }

  for (const [, g] of groups) setThinInstanceMatrices(g.master, g.mats);

  // --- bond style helpers
  const defaultBondStyles = (pairKey) => {
    if (pairKey === "C-C") return { radius: 0.12, color: new BABYLON.Color3(0.35, 0.35, 0.38) };
    if (pairKey === "C-H") return { radius: 0.07, color: new BABYLON.Color3(0.78, 0.82, 0.90) };
    return { radius: 0.09, color: new BABYLON.Color3(0.70, 0.73, 0.78) };
  };
  function pairKeyOf(i, j) {
    const a = atoms[i].element, b = atoms[j].element;
    return [a, b].sort().join("-");
  }
  function bondMatrix(aPos, bPos, radius) {
    const mid = aPos.add(bPos).scale(0.5);
    const v = bPos.subtract(aPos);
    const len = v.length();
    const rot = quatYto(v);
    const scale = new BABYLON.Vector3(radius * 2, len, radius * 2);
    return BABYLON.Matrix.Compose(scale, rot, mid);
  }

  // --- bond groups (per element pair)
  const bondGroups = new Map(); // pairKey -> { master, mats: [], style, indices: [] }

  function getOrCreateBondGroup(pairKey) {
    if (bondGroups.has(pairKey)) return bondGroups.get(pairKey);
    const style = { ...defaultBondStyles(pairKey), ...(bondStyles[pairKey] || {}) };
    const mat = new BABYLON.StandardMaterial(`matBond_${pairKey}`, scene);
    mat.diffuseColor = style.color.clone();
    mat.emissiveColor = style.color.scale(0.05);

    const master = BABYLON.MeshBuilder.CreateCylinder(`bond_${pairKey}`, {
      height: 1, diameter: 1, tessellation: 20
    }, scene);
    master.material = mat;
    master.isPickable = true;                 // pick bonds
    master.thinInstanceEnablePicking = true;  // pick specific instances
    if (debugAlwaysActive) master.alwaysSelectAsActiveMesh = true;

    const group = { master, mats: [], style, indices: [] };
    bondGroups.set(pairKey, group);
    return group;
  }

  // --- (A) initial bond detection or use provided bonds
  let bondList = bonds || [];
  if (!bonds) {
    for (let i = 0; i < atoms.length; i++) {
      const ai = atoms[i]; const infoi = elInfo(ai.element);
      for (let j = i + 1; j < atoms.length; j++) {
        const aj = atoms[j]; const infoj = elInfo(aj.element);
        const cutoff = (infoi.covRad + infoj.covRad) * bondScale;
        if (BABYLON.Vector3.Distance(ai.pos, aj.pos) <= cutoff) bondList.push({ i, j });
      }
    }
  }

  function rebuildBondGroupsFromList(list) {
    // reset mats/indices but keep masters/materials
    for (const [, g] of bondGroups) { g.mats = []; g.indices = []; }
    // fill mats per pair
    for (const { i, j } of list) {
      const key = pairKeyOf(i, j);
      const g = getOrCreateBondGroup(key);
      g.indices.push({ i, j });
      g.mats.push(bondMatrix(atoms[i].pos, atoms[j].pos, g.style.radius));
    }
    // bulk upload per group
    for (const [, g] of bondGroups) {
      setThinInstanceMatrices(g.master, g.mats);
    }
  }

  rebuildBondGroupsFromList(bondList);

  // --- per-element index mapping for interaction
  const elementToPerMeshIndex = new Map();
  for (const [elem, g] of groups) {
    const idxMap = new Map();
    g.indices.forEach((globalIdx, k) => idxMap.set(globalIdx, k));
    elementToPerMeshIndex.set(elem, idxMap);
  }

  // atoms presented to interaction layer (include scale)
  const atomsAPI = atoms.map((a, globalIdx) => {
    const elem = a.element;
    const g = groups.get(elem);
    const perMeshIdx = elementToPerMeshIndex.get(elem).get(globalIdx);
    return { mesh: g.master, type: elem, index: perMeshIdx, pos: a.pos, scale: perAtomScale[globalIdx] };
  });

  // --- (B) fast & robust: recompute matrices for all bonds and bulk re-upload
  function refreshBonds() {
    for (const [, g] of bondGroups) {
      // rebuild matrices based on current atom positions
      const newMats = new Array(g.indices.length);
      for (let k = 0; k < g.indices.length; k++) {
        const { i, j } = g.indices[k];
        newMats[k] = bondMatrix(atoms[i].pos, atoms[j].pos, g.style.radius);
      }
      g.mats = newMats;
      setThinInstanceMatrices(g.master, g.mats); // bulk upload + bounds refresh
    }
  }

  // --- (C) slow path: recompute connectivity by cutoff, then rebuild buffers
  function recomputeBonds({ hysteresis = 1.02 } = {}) {
    const newList = [];
    for (let i = 0; i < atoms.length; i++) {
      const ai = atoms[i]; const infoi = elInfo(ai.element);
      for (let j = i + 1; j < atoms.length; j++) {
        const aj = atoms[j]; const infoj = elInfo(aj.element);
        const cutoff = (infoi.covRad + infoj.covRad) * bondScale * hysteresis;
        if (BABYLON.Vector3.Distance(ai.pos, aj.pos) <= cutoff) newList.push({ i, j });
      }
    }
    bondList = newList;
    rebuildBondGroupsFromList(bondList);
  }

  // --- update a single atom's matrix via global index (kept for completeness)
  function updateAtom(globalIdx, newPos) {
    const elem = atoms[globalIdx].element;
    const g = groups.get(elem);
    const perMeshIdx = elementToPerMeshIndex.get(elem).get(globalIdx);
    const d = perAtomScale[globalIdx];

    atoms[globalIdx].pos.copyFrom(newPos);
    const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(d, d, d), BABYLON.Quaternion.Identity(), newPos);

    g.mats[perMeshIdx] = mat;
    setThinInstanceMatrices(g.master, g.mats);
  }

  // --- direct per-element/per-index matrix setter (used by interaction)
  function updateAtomMatrixByElement(elem, perMeshIdx, mat) {
    const g = groups.get(elem);
    if (!g) return;
    g.mats[perMeshIdx] = mat;
    setThinInstanceMatrices(g.master, g.mats);
  }

  // --- refresh all atom positions based on current atom.pos values
  function refreshAtoms() {
    for (const [elem, g] of groups) {
      // Rebuild matrices for all atoms of this element
      for (let k = 0; k < g.indices.length; k++) {
        const globalIdx = g.indices[k];
        const pos = atoms[globalIdx].pos;
        const d = perAtomScale[globalIdx];
        const mat = BABYLON.Matrix.Compose(
          new BABYLON.Vector3(d, d, d),
          BABYLON.Quaternion.Identity(),
          pos
        );
        g.mats[k] = mat;
      }
      setThinInstanceMatrices(g.master, g.mats);
    }
  }

  return {
    groups,
    atoms: atomsAPI,
    bonds: bondList,
    bondGroups,
    refreshBonds,             // bulk updates each group
    refreshAtoms,             // bulk updates all atom positions
    recomputeBonds,
    updateAtom,
    updateAtomMatrixByElement
  };
}
