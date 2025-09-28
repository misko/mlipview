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

  // Will hold the final returned molecule object once assembled so internal helpers
  // (defined earlier) can safely reference cell/bond state without using undefined 'mol'.
  let moleculeRef = null;

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
    // Local self reference (the returned molecule object) is 'self'; ensure we don't
    // reference 'mol' here because buildMolecule wraps/augments later.
    const self = moleculeRef || null;
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

    // If periodic cell active, perform periodic-union bond detection:
    // We generate candidate positions for each atom shifted by ±a,±b,±c (single vectors only)
    // then compute connectivity, project pairs back into primary by minimum-image mapping.
    let bondsRaw;
  const cellState = self && self.__cellState ? self.__cellState : null;
    const usePeriodic = !!(cellState && cellState.visible && cellState.vectors && cellState.vectors.a && cellState.vectors.b && cellState.vectors.c);
    if (usePeriodic) {
      try {
        const { a:avec, b:bvec, c:cvec } = cellState.vectors;
        // Build shift list (±a, ±b, ±c) plus origin
        const shifts = [
          new BABYLON.Vector3(0,0,0),
          avec, bvec, cvec,
          avec.scale(-1), bvec.scale(-1), cvec.scale(-1)
        ];
        // Original atom data
        const baseAtoms = atoms.map(a => ({ element: a.element, pos: [a.pos.x, a.pos.y, a.pos.z] }));
        // Collect augmented atom list with references to original index
        const aug = [];
        for (let si=0; si<shifts.length; si++) {
          const S = shifts[si];
          for (let i=0;i<atoms.length;i++) {
            const p = atoms[i].pos.add(S);
            aug.push({ element: atoms[i].element, pos: [p.x, p.y, p.z], baseIndex: i, shiftIndex: si });
          }
        }
        // Compute bonds on augmented set
        const augBonds = computeBondsNoState(aug.map(a => ({ element: a.element, pos: a.pos })));
        // Helper: build inverse matrix for fractional conversion (to choose minimum-image)
        const a=avec,b=bvec,c=cvec;
        const M=[a.x,b.x,c.x,a.y,b.y,c.y,a.z,b.z,c.z];
        const det = (
          M[0]*(M[4]*M[8]-M[5]*M[7]) -
          M[1]*(M[3]*M[8]-M[5]*M[6]) +
          M[2]*(M[3]*M[7]-M[4]*M[6])
        );
        let inv=null;
        if (Math.abs(det)>1e-14) {
          const invDet=1/det;
          inv=[
            (M[4]*M[8]-M[5]*M[7])*invDet,
            (M[2]*M[7]-M[1]*M[8])*invDet,
            (M[1]*M[5]-M[2]*M[4])*invDet,
            (M[5]*M[6]-M[3]*M[8])*invDet,
            (M[0]*M[8]-M[2]*M[6])*invDet,
            (M[2]*M[3]-M[0]*M[5])*invDet,
            (M[3]*M[7]-M[4]*M[6])*invDet,
            (M[1]*M[6]-M[0]*M[7])*invDet,
            (M[0]*M[4]-M[1]*M[3])*invDet
          ];
        }
        function frac(p) {
          if (!inv) return {u:p.x,v:p.y,w:p.z};
          return {
            u: inv[0]*p.x + inv[1]*p.y + inv[2]*p.z,
            v: inv[3]*p.x + inv[4]*p.y + inv[5]*p.z,
            w: inv[6]*p.x + inv[7]*p.y + inv[8]*p.z
          };
        }
        function wrapFrac(f) {
          return { u: f.u-Math.floor(f.u), v: f.v-Math.floor(f.v), w: f.w-Math.floor(f.w) };
        }
        function cart(f) {
          return new BABYLON.Vector3(
            a.x*f.u + b.x*f.v + c.x*f.w,
            a.y*f.u + b.y*f.v + c.y*f.w,
            a.z*f.u + b.z*f.v + c.z*f.w
          );
        }
        const seen = new Set();
        bondsRaw = [];
        for (const eb of augBonds) {
          const A = aug[eb.i];
            const B = aug[eb.j];
          // Skip bonds entirely inside the same shift that's not origin to reduce duplicates; keep those where at least one endpoint is origin OR different shifts.
          if (A.shiftIndex === B.shiftIndex && A.shiftIndex !== 0) continue;
          // Map endpoints to fractional, wrap into primary to get canonical midpoint & distance
          const pA = new BABYLON.Vector3(A.pos[0],A.pos[1],A.pos[2]);
          const pB = new BABYLON.Vector3(B.pos[0],B.pos[1],B.pos[2]);
          const fA_raw = frac(pA); const fB_raw = frac(pB);
          const fA = wrapFrac(fA_raw);
          const fB = wrapFrac(fB_raw);
          // Reconstruct canonical cartesian positions (inside primary) for distance
          const cA = cart(fA); const cB = cart(fB);
          const d = BABYLON.Vector3.Distance(cA, cB);
          // Key independent of ordering: base indices + rounded fractional delta (to prune duplicates)
          const di = Math.min(A.baseIndex, B.baseIndex);
          const dj = Math.max(A.baseIndex, B.baseIndex);
          // fractional separation (minimum image) for keying
          const df = { du: fB.u - fA.u, dv: fB.v - fA.v, dw: fB.w - fA.w };
          // normalize into [-0.5,0.5) for uniqueness
          df.du -= Math.round(df.du); df.dv -= Math.round(df.dv); df.dw -= Math.round(df.dw);
          const key = di+"_"+dj+":"+df.du.toFixed(4)+","+df.dv.toFixed(4)+","+df.dw.toFixed(4);
          if (seen.has(key)) continue;
          seen.add(key);
          // Detect crossing: if raw fractional difference differs from wrapped delta significantly
          let crossing = false;
          const rdu = fB_raw.u - fA_raw.u;
          const rdv = fB_raw.v - fA_raw.v;
          const rdw = fB_raw.w - fA_raw.w;
          function normComp(x){ return x - Math.round(x); }
          const ndu = normComp(rdu), ndv = normComp(rdv), ndw = normComp(rdw);
          if (Math.abs(ndu - rdu) > 1e-6 || Math.abs(ndv - rdv) > 1e-6 || Math.abs(ndw - rdw) > 1e-6) crossing = true;
          const baseOpacity = 1.0;
          const finalOpacity = crossing ? Math.min(0.5, baseOpacity) : baseOpacity;
          bondsRaw.push({ i: di, j: dj, length: d, opacity: finalOpacity, crossing });
        }
        // Fallback: if periodic path produced no bonds (edge case), revert
        if (!bondsRaw.length) {
          bondsRaw = computeBondsNoState(baseAtoms).map(b=>({ i:b.i,j:b.j,length:b.length, opacity:b.opacity }));
        }
      } catch (e) {
        console.warn('[Cell] periodic bond build failed, fallback to non-periodic', e.message);
        bondsRaw = computeBondsNoState(atoms.map(a => ({ element: a.element, pos: [a.pos.x, a.pos.y, a.pos.z] })));
      }
    } else {
      // Non-periodic: standard computation
      bondsRaw = computeBondsNoState(atoms.map(a => ({ element: a.element, pos: [a.pos.x, a.pos.y, a.pos.z] })));
    }

    // Map to existing rendering grouping (pairKey) and feed matrices
    const bonds = bondsRaw.map(b => ({
      i: b.i,
      j: b.j,
      opacity: b.opacity ?? 1.0,
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

  const _moleculeObject = {
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
  moleculeRef = _moleculeObject;
  return _moleculeObject;
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
      // Always have a vectors object; will be initialized from current bounds if null
      vectors: null, // {a: Vector3, b: Vector3, c: Vector3}
      arrows: [],    // arrow meshes
      ghostAtomGroups: new Map(), // element -> { master, mats: [] }
      ghostBondGroups: new Map(), // pairKey -> { master, mats: [] }
      visible: false, // visibility toggle (render + backend inclusion)
      _didInitialCenter: false,
      // originOffset: Cartesian shift applied to visual cell constructs (arrows, ghosts) instead of mutating atom positions.
      originOffset: new BABYLON.Vector3(0,0,0)
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
      const off = cellState.originOffset || BABYLON.Vector3.Zero();
      for (const [name, v] of [["cell_a", a],["cell_b", b],["cell_c", c]]) {
        const arrow = cell_makeArrow(name, v);
        if (arrow) {
          try { arrow.position.addInPlace(off); } catch {}
          cellState.arrows.push(arrow);
        }
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
  mat.alpha = 0.3; // 30% transparency (updated requirement)
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
  mat.alpha = 0.3; // 30% transparency (updated requirement)
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
      const { a, b, c } = cellState.vectors;
      const off = cellState.originOffset || BABYLON.Vector3.Zero();
      // Exactly 6 neighbor shifts (no origin, no diagonals)
      const shifts = [a, a.scale(-1), b, b.scale(-1), c, c.scale(-1)].filter(v => v && v.length() > 1e-9);
      if (!shifts.length) return;
      // Reset existing mats
      for (const [,g] of cellState.ghostAtomGroups) { g.mats = []; }
      for (const [,g] of cellState.ghostBondGroups) { g.mats = []; }

      // Atoms
      for (let gi=0; gi<mol.atoms.length; gi++) {
        const a = mol.atoms[gi];
        const scale = a.scale || 1;
        const basePos = a.pos.add(off);
        for (const shift of shifts) {
          const pos = basePos.add(shift);
          const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(scale, scale, scale), BABYLON.Quaternion.Identity(), pos);
          const g = cell_getOrCreateGhostAtomGroup(a.type);
          g.mats.push(mat);
        }
      }
      for (const [,g] of cellState.ghostAtomGroups) { setThinInstanceMatrices(g.master, g.mats); }

      // Bonds: translate existing primary bond matrices instead of recompute geometry
      for (const [pairKey, bg] of mol.bondGroups) {
        if (!bg.mats || !bg.mats.length) continue;
        const gg = cell_getOrCreateGhostBondGroup(pairKey);
        for (const mat of bg.mats) {
          // Extract original translation
          const m = mat.clone();
          const origPos = new BABYLON.Vector3(m.m[12], m.m[13], m.m[14]).add(off);
          for (const shift of shifts) {
            const newPos = origPos.add(shift);
            const clone = m.clone();
            clone.m[12] = newPos.x; clone.m[13] = newPos.y; clone.m[14] = newPos.z;
            gg.mats.push(clone);
          }
        }
        setThinInstanceMatrices(gg.master, gg.mats);
        try {
          if (gg.master && gg.master.material) {
            gg.master.material.alpha = Math.min(gg.master.material.alpha ?? 0.3, 0.5);
            gg.master.material.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
          }
        } catch {}
      }
      cell_updateVisibility();
    }

    function cell_updateVisibility() {
      const vis = cellState.visible && !!cellState.vectors;
      for (const m of cellState.arrows) { m.setEnabled(vis); }
      for (const [,g] of cellState.ghostAtomGroups) { g.master.setEnabled(vis); }
      for (const [,g] of cellState.ghostBondGroups) { g.master.setEnabled(vis); }
    }

  function cell_setVectors(aVec, bVec, cVec, { recenter = true } = {}) {
      // Accept array [x,y,z] or BABYLON.Vector3
      function toV(v) {
        if (!v) return new BABYLON.Vector3(0,0,0);
        if (v instanceof BABYLON.Vector3) return v.clone();
        if (Array.isArray(v) && v.length===3) return new BABYLON.Vector3(v[0],v[1],v[2]);
        if (typeof v === 'object' && 'x' in v) return new BABYLON.Vector3(v.x,v.y,v.z);
        return new BABYLON.Vector3(0,0,0);
      }
      cellState.vectors = { a: toV(aVec), b: toV(bVec), c: toV(cVec) };
      console.log('[Cell][setVectors]', {
        a: cellState.vectors.a.toString(),
        b: cellState.vectors.b.toString(),
        c: cellState.vectors.c.toString(),
        recenter,
        alreadyCentered: cellState._didInitialCenter
      });
      // Compute origin offset once; do not move atoms (prevents visible jump when toggling cell)
      if (!cellState._didInitialCenter && recenter) try {
        const { a, b, c } = cellState.vectors;
        const M = [a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z];
        const det = (
          M[0]*(M[4]*M[8]-M[5]*M[7]) -
          M[1]*(M[3]*M[8]-M[5]*M[6]) +
          M[2]*(M[3]*M[7]-M[4]*M[6])
        );
        if (Math.abs(det) > 1e-14 && mol.atoms && mol.atoms.length) {
          const invDet = 1/det;
          const inv = [
            (M[4]*M[8]-M[5]*M[7])*invDet,
            (M[2]*M[7]-M[1]*M[8])*invDet,
            (M[1]*M[5]-M[2]*M[4])*invDet,
            (M[5]*M[6]-M[3]*M[8])*invDet,
            (M[0]*M[8]-M[2]*M[6])*invDet,
            (M[2]*M[3]-M[0]*M[5])*invDet,
            (M[3]*M[7]-M[4]*M[6])*invDet,
            (M[1]*M[6]-M[0]*M[7])*invDet,
            (M[0]*M[4]-M[1]*M[3])*invDet
          ];
          function toFrac(pos) {
            const x=pos.x,y=pos.y,z=pos.z;
            return {
              u: inv[0]*x + inv[1]*y + inv[2]*z,
              v: inv[3]*x + inv[4]*y + inv[5]*z,
              w: inv[6]*x + inv[7]*y + inv[8]*z
            };
          }
          let uc=0, vc=0, wc=0;
            for (const atom of mol.atoms) {
              const f = toFrac(atom.pos);
              uc += f.u; vc += f.v; wc += f.w;
            }
          uc /= mol.atoms.length; vc /= mol.atoms.length; wc /= mol.atoms.length;
          // Desired center (0.5,0.5,0.5); shift delta in fractional then convert to cartesian shift
          const du = 0.5 - uc, dv = 0.5 - vc, dw = 0.5 - wc;
          if (Math.abs(du)+Math.abs(dv)+Math.abs(dw) > 1e-9) {
            const shift = new BABYLON.Vector3(
              a.x*du + b.x*dv + c.x*dw,
              a.y*du + b.y*dv + c.y*dw,
              a.z*du + b.z*dv + c.z*dw
            );
            console.log('[Cell][centering] computed originOffset (atoms unchanged)', { uc, vc, wc, du, dv, dw, shift: shift.toString() });
            cellState.originOffset.addInPlace(shift);
            cellState._didInitialCenter = true;
          }
        }
      } catch(e) { console.warn('[Cell] centering failed', e.message); }
      cell_buildArrows();
      cell_rebuildGhosts();
      // Do not auto-force visibility; user can toggle. Cell semantics always exist logically.
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
      console.log('[Cell][toggle]', { visible: cellState.visible, didInitialCenter: cellState._didInitialCenter });
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
      cell_setVectors(aVec,bVec,cVec,{ recenter:true });
    }

    // Expose API on molecule instance
  mol.setCellVectors = cell_setVectors;
  mol.clearCell = cell_clear; // retains vectors? For always-present semantics could noop; left as-is.
  mol.toggleCell = cell_toggle;
    mol.buildDefaultCell = cell_defaultFromBounds;
    mol.__cellState = cellState;

    // Wrap all atom positions into the primary cell (origin-based parallelepiped spanned by a,b,c)
    // Uses fractional coordinate mapping and modulo 1 wrapping. Returns true if any atom wrapped.
    mol.wrapIntoCell = function({ refresh = true, bondRefresh = false, log = true } = {}) {
      const cs = cellState;
      if (!cs || !cs.vectors || !cs.vectors.a || !cs.vectors.b || !cs.vectors.c) return false;
      const { a, b, c } = cs.vectors;
      // Build inverse matrix once
      const M = [a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z];
      const det = (
        M[0]*(M[4]*M[8]-M[5]*M[7]) -
        M[1]*(M[3]*M[8]-M[5]*M[6]) +
        M[2]*(M[3]*M[7]-M[4]*M[6])
      );
      if (Math.abs(det) < 1e-14) return false; // degenerate cell
      const invDet = 1/det;
      const inv = [
        (M[4]*M[8]-M[5]*M[7])*invDet,
        (M[2]*M[7]-M[1]*M[8])*invDet,
        (M[1]*M[5]-M[2]*M[4])*invDet,
        (M[5]*M[6]-M[3]*M[8])*invDet,
        (M[0]*M[8]-M[2]*M[6])*invDet,
        (M[2]*M[3]-M[0]*M[5])*invDet,
        (M[3]*M[7]-M[4]*M[6])*invDet,
        (M[1]*M[6]-M[0]*M[7])*invDet,
        (M[0]*M[4]-M[1]*M[3])*invDet
      ];
      function toFrac(pos) {
        const x=pos.x,y=pos.y,z=pos.z;
        return {
          u: inv[0]*x + inv[1]*y + inv[2]*z,
          v: inv[3]*x + inv[4]*y + inv[5]*z,
          w: inv[6]*x + inv[7]*y + inv[8]*z
        };
      }
      let wrappedCount = 0;
      const tol = 1e-9;
      for (const atom of mol.atoms) {
        const f = toFrac(atom.pos);
        let u=f.u, v=f.v, w=f.w;
        // Bring into [0,1)
        let changed=false;
        const u0=u; const v0=v; const w0=w;
        u = u - Math.floor(u);
        v = v - Math.floor(v);
        w = w - Math.floor(w);
        if (Math.abs(u-u0)>tol || Math.abs(v-v0)>tol || Math.abs(w-w0)>tol) changed=true;
        if (!changed) continue;
        // back to cartesian
        const nx = a.x*u + b.x*v + c.x*w;
        const ny = a.y*u + b.y*v + c.y*w;
        const nz = a.z*u + b.z*v + c.z*w;
        atom.pos.x = nx; atom.pos.y = ny; atom.pos.z = nz;
        wrappedCount++;
      }
      if (wrappedCount>0) {
        if (mol.changeCounter !== undefined) mol.changeCounter++;
        if (refresh && typeof mol.refreshAtoms === 'function') mol.refreshAtoms();
        if (bondRefresh && typeof mol.refreshBonds === 'function') mol.refreshBonds();
        // Retile ghosts immediately if cell visible. If refreshAtoms() was just called it
        // already triggers a rebuild; otherwise do it here to keep replicas in sync.
        if (cellState.visible && cellState.vectors) {
          if (!refresh) {
            try { cell_rebuildGhosts(); } catch(e) { console.warn('[Cell] ghost rebuild after wrap failed', e.message); }
          }
        }
        if (log) console.log('[Cell] wrapped', wrappedCount, 'atoms into primary cell' + (cellState.visible ? ' (ghosts retiled)' : ''));
        return true;
      }
      return false;
    };

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

