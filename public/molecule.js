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
    // Use chemistry-aware updating for all bond refreshes
    updateBondsWithChemistry();
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

  // --- Enhanced bond updating with realistic chemistry and transparency
  function updateBondsWithChemistry() {
    // Clear existing bonds and reset materials
    for (const [, g] of bondGroups) { 
      g.mats = []; 
      g.indices = [];
      // Reset material to opaque
      if (g.master && g.master.material) {
        g.master.material.alpha = 1.0;
        g.master.material.transparencyMode = BABYLON.Material.MATERIAL_OPAQUE;
      }
    }
    
    // Group bonds by their transparency level
    const bondsByOpacity = new Map(); // opacity -> bonds[]
    
    for (let i = 0; i < atoms.length; i++) {
      const ai = atoms[i]; const infoi = elInfo(ai.element);
      for (let j = i + 1; j < atoms.length; j++) {
        const aj = atoms[j]; const infoj = elInfo(aj.element);
        
        const distance = BABYLON.Vector3.Distance(ai.pos, aj.pos);
        const expectedBondLength = (infoi.covRad + infoj.covRad) * bondScale;
        
        // Calculate transparency based on deviation from expected length
        const minVisibleLength = expectedBondLength * 0.9; // 90% of expected
        const maxVisibleLength = expectedBondLength * 1.1; // 110% of expected
        
        let opacity = 1.0;
        
        if (distance < minVisibleLength) {
          // Fade out when compressed below 90%
          const compressionFactor = distance / minVisibleLength;
          opacity = compressionFactor; // Fades to 0 as distance approaches 0
        } else if (distance > maxVisibleLength) {
          // Fade out when stretched beyond 110%
          const stretchLimit = expectedBondLength * 1.5; // Complete disappearance at 150%
          if (distance > stretchLimit) {
            opacity = 0.0; // Completely gone
          } else {
            const stretchFactor = (distance - maxVisibleLength) / (stretchLimit - maxVisibleLength);
            opacity = 1.0 - stretchFactor; // Fade from 1.0 to 0.0
          }
        }
        
        // Only show bonds with some visibility
        if (opacity > 0.05) {
          const opacityKey = Math.round(opacity * 20) / 20; // Round to nearest 0.05
          if (!bondsByOpacity.has(opacityKey)) {
            bondsByOpacity.set(opacityKey, []);
          }
          bondsByOpacity.get(opacityKey).push({ i, j, pairKey: pairKeyOf(i, j) });
        }
      }
    }
    
    // Create bonds grouped by opacity for efficiency
    for (const [opacity, bonds] of bondsByOpacity) {
      for (const { i, j, pairKey } of bonds) {
        const group = getOrCreateBondGroup(pairKey);
        
        // Set material transparency for this opacity level
        if (group.master && group.master.material) {
          if (opacity < 1.0) {
            group.master.material.alpha = Math.min(group.master.material.alpha || 1.0, opacity);
            group.master.material.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
          }
        }
        
        const mat = bondMatrix(atoms[i].pos, atoms[j].pos, group.style.radius);
        group.mats.push(mat);
        group.indices.push({ i, j });
      }
    }
    
    // Update all bond group meshes
    for (const [, g] of bondGroups) {
      setThinInstanceMatrices(g.master, g.mats);
    }
    
    const totalBonds = Array.from(bondsByOpacity.values()).reduce((sum, bonds) => sum + bonds.length, 0);
    console.log(`[Bonds] Updated with ${totalBonds} bonds (${bondsByOpacity.size} opacity levels)`);
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
    updateBondsWithChemistry, // realistic chemistry-based bonds with transparency
    updateAtom,
    updateAtomMatrixByElement
  };
}
