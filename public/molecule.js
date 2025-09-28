import { elInfo } from "./elements.js";
import { quatYto, setThinInstanceMatrices } from "./helpers.js";
import { computeBondsNoState } from "./bond_render.js";

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
    bondStyles = {},
    center = true // center molecule on load (geometric center of positions)
  } = opts || {};

  // Optional centering (currently geometric center; mass-weighting could be added later)
  if (center && atoms.length) {
    let cx = 0, cy = 0, cz = 0;
    for (const a of atoms) { cx += a.pos.x; cy += a.pos.y; cz += a.pos.z; }
    cx /= atoms.length; cy /= atoms.length; cz /= atoms.length;
    if (Math.abs(cx) + Math.abs(cy) + Math.abs(cz) > 1e-9) {
      for (const a of atoms) { a.pos.x -= cx; a.pos.y -= cy; a.pos.z -= cz; }
    }
  }

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
  // Minimize visual artifact: scale master host to epsilon so it never appears
  // as a white cylinder at origin, while keeping it 'visible' for the thin instance system.
  master.scaling = new BABYLON.Vector3(1, 1, 1); // normal scale for proper bond length reference
  master.isVisible = true; // visible when instances present, may be hidden dynamically if empty
  master.isPickable = false;                // picking occurs via thin instances
  master.thinInstanceEnablePicking = true;  // allow per-instance picking
  if (debugAlwaysActive) master.alwaysSelectAsActiveMesh = true;

    const group = { master, mats: [], style, indices: [] };
    bondGroups.set(pairKey, group);
    return group;
  }

  // --- (A) initial bond generation unified with runtime logic.
  // If explicit bonds are provided, we respect them; otherwise we derive
  // an initial list using computeBondsNoState (same logic used later).
  let bondList = [];
  if (bonds) {
    bondList = bonds.slice();
  } else {
    const initial = computeBondsNoState(atoms.map(a => ({ element: a.element, pos: [a.pos.x, a.pos.y, a.pos.z] })));
    bondList = initial.map(b => ({ i: b.i, j: b.j }));
  }
  // Use the normal chemistry update path to build thin instances (ensures opacity treatment)
  // but since updateBondsWithChemistry reads bondList only for historical purposes now,
  // we still perform one direct population pass with full opacity for speed, then immediately
  // call updateBondsWithChemistry to settle final visuals.
  (function initialPopulate() {
    for (const [, g] of bondGroups) { g.mats = []; g.indices = []; }
    for (const { i, j } of bondList) {
      const key = pairKeyOf(i, j);
      const g = getOrCreateBondGroup(key);
      g.indices.push({ i, j });
      g.mats.push(bondMatrix(atoms[i].pos, atoms[j].pos, g.style.radius));
    }
    for (const [, g] of bondGroups) setThinInstanceMatrices(g.master, g.mats);
    // Now immediately apply unified opacity/weight logic
    updateBondsWithChemistry();
  })();

  // --- per-element index mapping for interaction
  const elementToPerMeshIndex = new Map();
  for (const [elem, g] of groups) {
    const idxMap = new Map();
    g.indices.forEach((globalIdx, k) => idxMap.set(globalIdx, k));
    elementToPerMeshIndex.set(elem, idxMap);
  }

  // ring detection now handled inside computeBonds (pure function)

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
    // Re-run lightweight proximity pass to seed a new list, then apply chemistry logic
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
    updateBondsWithChemistry();
  }

  // --- Enhanced bond updating with realistic chemistry and transparency
  function updateBondsWithChemistry() {
    // Reset existing bond group buffers/materials
    for (const [, g] of bondGroups) {
      g.mats = [];
      g.indices = [];
      if (g.master && g.master.material) {
        g.master.material.alpha = 1.0;
        g.master.material.transparencyMode = BABYLON.Material.MATERIAL_OPAQUE;
      }
  // Do not upload an empty buffer here; wait until new mats are populated to avoid flicker.
    }

    // Use stateless bond computation from bond_render.js
    // Adapt atom vector3 -> plain tuple expected by computeBondsNoState
    const bondsRaw = computeBondsNoState(atoms.map(a => ({ element: a.element, pos: [a.pos.x, a.pos.y, a.pos.z] })));

    // Map to existing rendering grouping (pairKey) and feed matrices
    const bonds = bondsRaw.map(b => ({
      i: b.i,
      j: b.j,
      opacity: b.opacity,
      pairKey: [atoms[b.i].element, atoms[b.j].element].sort().join("-")
    }));

    // Track min alpha per group for coarse transparency
    const minAlphaPerGroup = new Map();
    for (const b of bonds) {
      const group = getOrCreateBondGroup(b.pairKey);
      const mat = bondMatrix(atoms[b.i].pos, atoms[b.j].pos, group.style.radius);
      group.mats.push(mat);
      group.indices.push({ i: b.i, j: b.j });
      if (b.opacity < 0.999) {
        const prev = minAlphaPerGroup.get(b.pairKey) ?? 1.0;
        if (b.opacity < prev) minAlphaPerGroup.set(b.pairKey, b.opacity);
      }
    }

    for (const [pairKey, g] of bondGroups) {
      const matAlpha = minAlphaPerGroup.get(pairKey);
      if (matAlpha !== undefined && g.master.material) {
        g.master.material.alpha = matAlpha;
        g.master.material.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
      }
      setThinInstanceMatrices(g.master, g.mats);
      // Toggle host visibility based on whether we have instances; prevents lone cylinder
      g.master.isVisible = g.mats.length > 0;
    }

    console.log(`[Bonds] Updated with ${bonds.length} bonds (bond_render.js integration)`);
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
    updateAtomMatrixByElement,
    // Cell / tiling API will be attached below
    setCellVectors: (a,b,c)=>cell_setVectors(a,b,c),
    clearCell: ()=>cell_clear(),
    toggleCell: (v)=>cell_toggle(v)
  };
}

// ---------------- Cell / Periodic Tile Support ------------------
// Injected by augmenting buildMolecule's closure AFTER export return object assembly above.
// We rely on function hoisting for helper definitions.

// NOTE: This patch adds three visual components when enabled:
// 1. Three yellow arrow cylinders (cell vectors a,b,c) starting at origin.
// 2. One translucent (50% alpha) replica of the molecule translated by +a, +b, +c (no combinations like a+b).
// 3. Matching translucent bond replicas for each translation.
// Ghost atoms/bonds are non-pickable and excluded from makePicker logic by using distinct mesh name prefixes.

// We attach the API to each molecule instance by mutating the returned object before it's used externally.

// To preserve original file structure minimalism, we monkey-patch buildMolecule after declaration if not already augmented.
if (!buildMolecule.__cellAugmented) {
  const original = buildMolecule;
  buildMolecule = function(scene, opts) {
    const mol = original(scene, opts);
    // internal state container
    const cellState = {
      vectors: null, // {a: Vector3, b: Vector3, c: Vector3}
      arrows: [],    // arrow meshes
      ghostAtomGroups: new Map(), // element -> { master, mats: [] }
      ghostBondGroups: new Map(), // pairKey -> { master, mats: [] }
      visible: false
    };

    function cell_disposeGhosts() {
      for (const [, g] of cellState.ghostAtomGroups) { try { g.master.dispose(); } catch {} }
      for (const [, g] of cellState.ghostBondGroups) { try { g.master.dispose(); } catch {} }
      cellState.ghostAtomGroups.clear();
      cellState.ghostBondGroups.clear();
    }
    function cell_disposeArrows() {
      for (const m of cellState.arrows) { try { m.dispose(); } catch {} }
      cellState.arrows = [];
    }
    function cell_makeArrow(name, v) {
      const len = v.length();
      if (len < 1e-6) return null;
      const mat = new BABYLON.StandardMaterial(name+"Mat", scene);
      mat.diffuseColor = new BABYLON.Color3(1,0.9,0.1);
      mat.emissiveColor = new BABYLON.Color3(1,0.9,0.1).scale(0.6);
      mat.alpha = 1.0;
      const shaft = BABYLON.MeshBuilder.CreateCylinder(name, { height: len, diameter: 0.08, tessellation: 20 }, scene);
      shaft.material = mat;
      shaft.isPickable = false;
      // orient along v (reuse quatYto logic from helpers via temporary computation)
      try {
        const rot = (function quatYtoVec(vec){
          const up = new BABYLON.Vector3(0,1,0);
            const axis = BABYLON.Vector3.Cross(up, vec);
            const len2 = axis.length();
            if (len2 < 1e-8) {
              if (vec.y > 0) return BABYLON.Quaternion.Identity();
              return new BABYLON.Quaternion(0,0,1,0); // 180 deg around Z
            }
            axis.normalize();
            const angle = Math.acos(BABYLON.Vector3.Dot(up, vec) / vec.length());
            return BABYLON.Quaternion.RotationAxis(axis, angle);
        })(v);
        shaft.rotationQuaternion = rot;
      } catch {}
      shaft.position = v.scale(0.5); // center cylinder between origin and end
      // Arrow head
      const head = BABYLON.MeshBuilder.CreateCylinder(name+"Head", { height: len*0.15, diameterTop: 0, diameterBottom: 0.20, tessellation: 20 }, scene);
      head.material = mat;
      head.isPickable = false;
      head.position = v.clone().scale(0.92);
      try { head.rotationQuaternion = shaft.rotationQuaternion.clone(); } catch {}
      const parent = new BABYLON.TransformNode(name+"Root", scene);
      shaft.parent = parent;
      head.parent = parent;
      parent.isPickable = false;
      return parent;
    }

    function cell_buildArrows() {
      cell_disposeArrows();
      if (!cellState.vectors) return;
      const { a, b, c } = cellState.vectors;
      for (const [name, v] of [["cell_a", a],["cell_b", b],["cell_c", c]]) {
        const arrow = cell_makeArrow(name, v);
        if (arrow) cellState.arrows.push(arrow);
      }
      cell_updateVisibility();
    }

    function cell_getOrCreateGhostAtomGroup(element) {
      if (cellState.ghostAtomGroups.has(element)) return cellState.ghostAtomGroups.get(element);
      const master = BABYLON.MeshBuilder.CreateSphere("ghost_atom_"+element, { diameter: 1, segments: 18 }, scene);
      const mat = new BABYLON.StandardMaterial("ghost_atomMat_"+element, scene);
      const baseInfo = mol.groups.get(element)?.info;
      const baseColor = baseInfo ? baseInfo.color : new BABYLON.Color3(0.8,0.8,0.8);
      mat.diffuseColor = baseColor.clone();
      mat.emissiveColor = baseColor.scale(0.03);
      mat.alpha = 0.5; // 50% transparency
      mat.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
      master.material = mat;
      master.isPickable = false;
      master.thinInstanceEnablePicking = false;
      const g = { master, mats: [] };
      cellState.ghostAtomGroups.set(element, g);
      return g;
    }
    function cell_getOrCreateGhostBondGroup(pairKey) {
      if (cellState.ghostBondGroups.has(pairKey)) return cellState.ghostBondGroups.get(pairKey);
      const master = BABYLON.MeshBuilder.CreateCylinder("ghost_bond_"+pairKey, { height: 1, diameter: 1, tessellation: 16 }, scene);
      const mat = new BABYLON.StandardMaterial("ghost_bondMat_"+pairKey, scene);
      mat.diffuseColor = new BABYLON.Color3(0.85,0.85,0.85);
      mat.emissiveColor = mat.diffuseColor.scale(0.04);
      mat.alpha = 0.5;
      mat.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
      master.material = mat;
      master.isPickable = false;
      master.thinInstanceEnablePicking = false;
      const g = { master, mats: [] };
      cellState.ghostBondGroups.set(pairKey, g);
      return g;
    }

    function cell_rebuildGhosts() {
      if (!cellState.vectors) return;
      const shifts = [cellState.vectors.a, cellState.vectors.b, cellState.vectors.c].filter(v=>v && v.length()>1e-9);
      if (!shifts.length) return;
      // Reset existing mats
      for (const [,g] of cellState.ghostAtomGroups) { g.mats = []; }
      for (const [,g] of cellState.ghostBondGroups) { g.mats = []; }

      // Atoms
      for (let gi=0; gi<mol.atoms.length; gi++) {
        const a = mol.atoms[gi];
        const scale = a.scale || 1;
        const basePos = a.pos;
        for (const shift of shifts) {
          const pos = basePos.add(shift);
          const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(scale, scale, scale), BABYLON.Quaternion.Identity(), pos);
          const g = cell_getOrCreateGhostAtomGroup(a.type);
          g.mats.push(mat);
        }
      }
      for (const [,g] of cellState.ghostAtomGroups) { setThinInstanceMatrices(g.master, g.mats); }

      // Bonds: reuse existing bondGroups definition from primary molecule
      for (const [pairKey, bg] of mol.bondGroups) {
        if (!bg.indices || !bg.indices.length) continue;
        const gg = cell_getOrCreateGhostBondGroup(pairKey);
        for (const pair of bg.indices) {
          const ai = mol.atoms[pair.i];
          const aj = mol.atoms[pair.j];
          if (!ai || !aj) continue;
          for (const shift of shifts) {
            const aPos = ai.pos.add(shift);
            const bPos = aj.pos.add(shift);
            const mid = aPos.add(bPos).scale(0.5);
            const v = bPos.subtract(aPos);
            const len = v.length();
            if (len < 1e-6) continue;
            const up = new BABYLON.Vector3(0,1,0);
            const axis = BABYLON.Vector3.Cross(up, v);
            const axisLen = axis.length();
            let rot = BABYLON.Quaternion.Identity();
            if (axisLen > 1e-8) {
              axis.normalize();
              const ang = Math.acos(BABYLON.Vector3.Dot(up, v) / len);
              rot = BABYLON.Quaternion.RotationAxis(axis, ang);
            } else if (v.y < 0) {
              rot = new BABYLON.Quaternion(0,0,1,0); // 180 deg
            }
            const radius = bg.style?.radius || 0.1;
            const scaleVec = new BABYLON.Vector3(radius*2, len, radius*2);
            const mat = BABYLON.Matrix.Compose(scaleVec, rot, mid);
            gg.mats.push(mat);
          }
        }
        setThinInstanceMatrices(gg.master, gg.mats);
      }
      cell_updateVisibility();
    }

    function cell_updateVisibility() {
      const vis = cellState.visible && !!cellState.vectors;
      for (const m of cellState.arrows) { m.setEnabled(vis); }
      for (const [,g] of cellState.ghostAtomGroups) { g.master.setEnabled(vis); }
      for (const [,g] of cellState.ghostBondGroups) { g.master.setEnabled(vis); }
    }

    function cell_setVectors(aVec, bVec, cVec) {
      // Accept array [x,y,z] or BABYLON.Vector3
      function toV(v) {
        if (!v) return new BABYLON.Vector3(0,0,0);
        if (v instanceof BABYLON.Vector3) return v.clone();
        if (Array.isArray(v) && v.length===3) return new BABYLON.Vector3(v[0],v[1],v[2]);
        if (typeof v === 'object' && 'x' in v) return new BABYLON.Vector3(v.x,v.y,v.z);
        return new BABYLON.Vector3(0,0,0);
      }
      cellState.vectors = { a: toV(aVec), b: toV(bVec), c: toV(cVec) };
      cell_buildArrows();
      cell_rebuildGhosts();
      cellState.visible = true;
      cell_updateVisibility();
    }
    function cell_clear() {
      cell_disposeArrows();
      cell_disposeGhosts();
      cellState.vectors = null;
      cellState.visible = false;
    }
    function cell_toggle(v) {
      if (typeof v === 'boolean') cellState.visible = v; else cellState.visible = !cellState.visible;
      if (cellState.visible && cellState.vectors) {
        cell_rebuildGhosts();
      }
      cell_updateVisibility();
      return cellState.visible;
    }
    function cell_defaultFromBounds() {
      // Compute axis-aligned bounds of original molecule atoms
      if (!mol.atoms || !mol.atoms.length) return;
      let minX=Infinity,minY=Infinity,minZ=Infinity,maxX=-Infinity,maxY=-Infinity,maxZ=-Infinity;
      for (const a of mol.atoms) {
        const p=a.pos; if (p.x<minX) minX=p.x; if (p.y<minY) minY=p.y; if (p.z<minZ) minZ=p.z;
        if (p.x>maxX) maxX=p.x; if (p.y>maxY) maxY=p.y; if (p.z>maxZ) maxZ=p.z;
      }
      const pad = 1.5; // Å padding for clarity
      const aVec = new BABYLON.Vector3((maxX-minX)+pad,0,0);
      const bVec = new BABYLON.Vector3(0,(maxY-minY)+pad,0);
      const cVec = new BABYLON.Vector3(0,0,(maxZ-minZ)+pad);
      cell_setVectors(aVec,bVec,cVec);
    }

    // Expose API on molecule instance
    mol.setCellVectors = cell_setVectors;
    mol.clearCell = cell_clear;
    mol.toggleCell = cell_toggle;
    mol.buildDefaultCell = cell_defaultFromBounds;
    mol.__cellState = cellState;

    // Patch refreshAtoms to keep ghosts in sync
    const oldRefreshAtoms = mol.refreshAtoms;
    mol.refreshAtoms = function() {
      oldRefreshAtoms();
      if (cellState.visible && cellState.vectors) {
        cell_rebuildGhosts();
      }
    };

    return mol;
  };
  buildMolecule.__cellAugmented = true;
}

