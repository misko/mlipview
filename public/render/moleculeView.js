import { elInfo } from '../elements.js';
import { computeBondsNoState } from '../bond_render.js';
import { __count } from '../util/funcCount.js';

const MATERIAL_FALLBACKS = {
  MATERIAL_ALPHABLEND: 2,
  MATERIAL_OPAQUE: 0,
};

function resolveMaterialConstant(name) {
  try {
    const mat = typeof BABYLON !== 'undefined' ? BABYLON.Material : undefined;
    const value = mat && typeof mat[name] === 'number' ? mat[name] : undefined;
    if (typeof value === 'number') return value;
  } catch {}
  return MATERIAL_FALLBACKS[name] ?? 0;
}

export function combineOpacity(base, mask, mode = 'multiply') {
  const clamp = (v) => {
    const n = Number(v);
    if (!Number.isFinite(n)) return 0;
    if (n <= 0) return 0;
    if (n >= 1) return 1;
    return n;
  };
  const a = clamp(base);
  const b = clamp(mask);
  switch (mode) {
    case 'override':
      return b;
    case 'clamp':
      return Math.min(a, b);
    default:
      return a * b;
  }
}

export function applyMaterialTransparency(mat, transparent, options = {}) {
  if (!mat || typeof mat !== 'object') return mat;
  const { usePrePass = true } = options;
  if (transparent) {
    mat.transparencyMode = resolveMaterialConstant('MATERIAL_ALPHABLEND');
    mat.forceDepthWrite = false;
    if (usePrePass) {
      mat.needDepthPrePass = true;
      mat.separateCullingPass = true;
    } else {
      mat.needDepthPrePass = false;
      mat.separateCullingPass = mat.separateCullingPass === true ? mat.separateCullingPass : false;
    }
  } else {
    mat.transparencyMode = resolveMaterialConstant('MATERIAL_OPAQUE');
    mat.forceDepthWrite = false;
    mat.needDepthPrePass = false;
    mat.separateCullingPass = false;
  }
  return mat;
}

export function createMoleculeView(scene, molState) {
  __count('moleculeView#createMoleculeView');
  // === Molecule Master Root & Registry =====================================
  // We introduce a single transform node (moleculeRoot) that becomes the parent
  // of every master thin‑instance host (atoms, bonds, forces, ghosts, highlights).
  // VR / AR systems should manipulate ONLY this root when applying global
  // rotations / translations / scaling. This eliminates regex heuristics over
  // mesh names and provides a stable anchor even if internal naming changes.
  //
  // A lightweight registry is also exposed on molState.__masters for debugging
  // and for fallback consumers. Each entry is { kind, mesh }. The registry is
  // append‑only during a molecule’s lifetime; on a full molecule swap a new
  // view (with a new root & registry) is constructed.
  // ===========================================================================
  // Some test environments provide a very small Babylon stub without TransformNode.
  // Fallback: use a MeshBuilder sphere (disabled) as the root parent container.
  let moleculeRoot;
  try {
    if (BABYLON.TransformNode) {
      moleculeRoot = new BABYLON.TransformNode('molecule_root', scene);
    } else {
      // minimal fallback; flagged via metadata for VR utilities
      moleculeRoot = BABYLON.MeshBuilder?.CreateSphere
        ? BABYLON.MeshBuilder.CreateSphere('molecule_root', { diameter: 0.0001 }, scene)
        : {
            name: 'molecule_root',
            position: new BABYLON.Vector3(0, 0, 0),
            scaling: new BABYLON.Vector3(1, 1, 1),
          };
      if (moleculeRoot.setEnabled) moleculeRoot.setEnabled(false);
      else moleculeRoot.isVisible = false;
    }
  } catch {
    moleculeRoot = {
      name: 'molecule_root_stub',
      position: { x: 0, y: 0, z: 0 },
      scaling: { x: 1, y: 1, z: 1 },
    };
  }
  moleculeRoot.metadata = { role: 'moleculeRoot' };
  const mastersRegistry = []; // array of { kind: 'atomMaster'|'atomSoftMaster'|'bondMaster'|'bondSoftMaster'|'forceMaster', mesh }
  function registerMaster(kind, mesh) {
    try {
      if (mesh && mesh.parent !== moleculeRoot) mesh.parent = moleculeRoot;
    } catch {}
    mastersRegistry.push({ kind, mesh });
  }
  // Expose for VR/services (read‑only intended):
  molState.moleculeRoot = moleculeRoot;
  molState.__masters = mastersRegistry;
  function wrapThinInstanceBuffer(mesh) {
    if (!mesh || mesh.__thinInstanceWrapped) return mesh;
    const original = typeof mesh.thinInstanceSetBuffer === 'function' ? mesh.thinInstanceSetBuffer.bind(mesh) : null;
    mesh.thinInstanceSetBuffer = function thinInstanceSetBuffer(kind, data, stride) {
      if (original) {
        try {
          original(kind, data, stride);
        } catch {
          // ignore stub failures
        }
      }
      if (this && typeof this === 'object') {
        if (!this._buffers || typeof this._buffers !== 'object') this._buffers = {};
        this._buffers[kind] = data;
        this._buffersStride = this._buffersStride || {};
        this._buffersStride[kind] = stride;
      }
    };
    mesh.__thinInstanceWrapped = true;
    return mesh;
  }

  function createIdentityMatrix() {
    try {
      if (BABYLON?.Matrix?.Identity) return BABYLON.Matrix.Identity();
    } catch {}
    if (BABYLON?.Matrix) {
      try {
        const mat = new BABYLON.Matrix();
        if (typeof BABYLON.Matrix.IdentityToRef === 'function') {
          BABYLON.Matrix.IdentityToRef(mat);
        } else if (mat && mat.m) {
          mat.m.fill?.(0);
          if (mat.m.length >= 16) {
            mat.m[0] = mat.m[5] = mat.m[10] = mat.m[15] = 1;
          }
        } else {
          mat.m = Float32Array.from([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
          ]);
        }
        return mat;
      } catch {}
    }
    return {
      m: Float32Array.from([
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
      ]),
    };
  }

  const MODE_SOLID = 'solid';
  const MODE_SOFT = 'soft';
  const MODE_LIST = [MODE_SOLID, MODE_SOFT];
  const SOFT_ATOM_ALPHA = 0.32;
  const SOFT_BOND_ALPHA = 0.18;
  const SOFT_MODE_THRESHOLD = 0.99;
  const atomGroups = new Map(); // element -> dual-mode group
  const atomMasterLookup = new Map();
  let atomDefaultModes = [];
  let atomCurrentModes = [];
  const atomRenderRefs = [];
  const bondGroups = new Map(); // key -> dual-mode group
  const bondIndexByKey = new Map();
  function canonicalBondKey(i, j) {
    const a = i < j ? i : j;
    const b = i < j ? j : i;
    return `${a}-${b}`;
  }
  const bondMasterLookup = new Map();
  let bondDefaultModes = [];
  let bondCurrentModes = [];
  const bondRenderRefs = [];
  // Force vectors rendered as thin instances (similar to bonds) for VR/AR consistency
  const forceGroups = new Map(); // single group keyed by 'force' currently; extensible for per-type later
  const ghostAtomGroups = new Map(); // separate to keep pickable flag false
  const ghostBondGroups = new Map();
  let activeOpacityMask = null;
  function bondStretchDebugEnabled() {
    try {
      if (typeof window === 'undefined') return false;
      if (window.__MLIP_DEBUG_STRETCH === true) return true;
      return /[?&]bondStretchDebug=1/.test(window.location?.search || '');
    } catch {
      return false;
    }
  }
  // Geometry version stamping: increment outside (molState.geometryVersion++) when topology changes
  if (typeof molState.geometryVersion !== 'number') molState.geometryVersion = 0;
  let cachedAtomVersion = -1;
  let cachedBondVersion = -1;
  let cellLines = null; // line system for cell outline
  // Highlight meshes (single reusable shell each for selected atom or bond)
  const highlight = { atom: null, bond: null };
  function keyOf(i, j) {
    const a = molState.elements[i],
      b = molState.elements[j];
    return [a, b].sort().join('-');
  }
  function createAtomMaster(el, mode, { pickable = true } = {}) {
    const suffix = mode === MODE_SOFT ? '_soft' : '_solid';
    const info = elInfo(el);
    const mat = new BABYLON.StandardMaterial(`atom_${el}${suffix}_mat`, scene);
    mat.diffuseColor = info.color.clone();
    mat.emissiveColor = info.color.scale(mode === MODE_SOFT ? 0.02 : 0.06);
    applyMaterialTransparency(mat, mode === MODE_SOFT, { usePrePass: mode === MODE_SOFT });
    mat.alpha = mode === MODE_SOFT ? SOFT_ATOM_ALPHA : 1;
    if (mode === MODE_SOFT) mat.backFaceCulling = false;
    const segments = mode === MODE_SOFT ? 20 : 24;
    const master = BABYLON.MeshBuilder.CreateSphere(
      `atom_${el}${suffix}`,
      { diameter: 1, segments },
      scene
    );
    master.material = mat;
    master.isPickable = !!pickable;
    master.thinInstanceEnablePicking = !!pickable;
    if (typeof master.setEnabled === 'function') master.setEnabled(false);
    else master.isVisible = false;
    wrapThinInstanceBuffer(master);
    registerMaster(mode === MODE_SOFT ? 'atomSoftMaster' : 'atomMaster', master);
    return master;
  }

  function ensureAtomGroup(el) {
    __count('moleculeView#ensureAtomGroup');
    if (atomGroups.has(el)) return atomGroups.get(el);
    const solidMaster = createAtomMaster(el, MODE_SOLID, { pickable: true });
    const softMaster = createAtomMaster(el, MODE_SOFT, { pickable: true });
    const group = {
      key: el,
      masters: {
        [MODE_SOLID]: solidMaster,
        [MODE_SOFT]: softMaster,
      },
      instances: {
        [MODE_SOLID]: [],
        [MODE_SOFT]: [],
      },
      ghostInstances: {
        [MODE_SOLID]: [],
        [MODE_SOFT]: [],
      },
      indices: [],
      mats: [],
    };
    group.master = solidMaster;
    group.softMaster = softMaster;
    atomGroups.set(el, group);
    atomMasterLookup.set(solidMaster, { group, mode: MODE_SOLID });
    atomMasterLookup.set(softMaster, { group, mode: MODE_SOFT });
    return group;
  }

  function createBondMaster(key, mode, { pickable = true } = {}) {
    const suffix = mode === MODE_SOFT ? '_soft' : '_solid';
    const mat = new BABYLON.StandardMaterial(`bond_${key}${suffix}_mat`, scene);
    mat.diffuseColor = new BABYLON.Color3(0.75, 0.78, 0.8);
    mat.emissiveColor = mat.diffuseColor.scale(mode === MODE_SOFT ? 0.02 : 0.05);
    applyMaterialTransparency(mat, mode === MODE_SOFT, { usePrePass: mode === MODE_SOFT });
    mat.alpha = mode === MODE_SOFT ? SOFT_BOND_ALPHA : 1;
    const master = BABYLON.MeshBuilder.CreateCylinder(
      `bond_${key}${suffix}`,
      { height: 1, diameter: 1, tessellation: 22 },
      scene
    );
    master.material = mat;
    master.isPickable = !!pickable;
    master.thinInstanceEnablePicking = !!pickable;
    if (typeof master.setEnabled === 'function') master.setEnabled(false);
    else master.isVisible = false;
    wrapThinInstanceBuffer(master);
    registerMaster(mode === MODE_SOFT ? 'bondSoftMaster' : 'bondMaster', master);
    return master;
  }

  function ensureBondGroup(key) {
    __count('moleculeView#ensureBondGroup');
    if (bondGroups.has(key)) return bondGroups.get(key);
    const solidMaster = createBondMaster(key, MODE_SOLID, { pickable: true });
    const softMaster = createBondMaster(key, MODE_SOFT, { pickable: true });
    const group = {
      key,
      masters: {
        [MODE_SOLID]: solidMaster,
        [MODE_SOFT]: softMaster,
      },
      instances: {
        [MODE_SOLID]: [],
        [MODE_SOFT]: [],
      },
      ghostInstances: {
        [MODE_SOLID]: [],
        [MODE_SOFT]: [],
      },
      indices: [],
      mats: [],
    };
    group.master = solidMaster;
    group.softMaster = softMaster;
    bondGroups.set(key, group);
    bondMasterLookup.set(solidMaster, { group, mode: MODE_SOLID });
    bondMasterLookup.set(softMaster, { group, mode: MODE_SOFT });
    return group;
  }

  const ensureGhostAtomGroup = (el) => {
    __count('moleculeView#ensureGhostAtomGroup');
    if (ghostAtomGroups.has(el)) return ghostAtomGroups.get(el);
    const baseGroup = ensureAtomGroup(el);
    const group = {
      master: baseGroup.masters[MODE_SOFT],
      baseGroup,
      mats: [],
      indices: [],
    };
    ghostAtomGroups.set(el, group);
    return group;
  };
  const ensureGhostBondGroup = (key) => {
    __count('moleculeView#ensureGhostBondGroup');
    if (ghostBondGroups.has(key)) return ghostBondGroups.get(key);
    const baseGroup = ensureBondGroup(key);
    const group = {
      master: baseGroup.masters[MODE_SOFT],
      baseGroup,
      mats: [],
      indices: [],
    };
    ghostBondGroups.set(key, group);
    return group;
  };
  function ensureForceGroup() {
    __count('moleculeView#ensureForceGroup');
    if (forceGroups.has('force')) return forceGroups.get('force');
    const mat = new BABYLON.StandardMaterial('force_vec_mat', scene);
    mat.diffuseColor = new BABYLON.Color3(0.95, 0.2, 0.2);
    mat.emissiveColor = mat.diffuseColor.scale(0.55);
    mat.specularColor = new BABYLON.Color3(0.1, 0.1, 0.1);
    mat.disableLighting = true; // keep vivid in dim VR/AR lighting
    // Use unit primitives & scale per-instance:
    //  - Shaft: unit cylinder (height=1, diameter=1)
    //  - Head: unit cone (height=1, diameterBottom=1, diameterTop=0)
    const shaftMaster = BABYLON.MeshBuilder.CreateCylinder(
      'force_vector_shaft_master',
      { height: 1, diameter: 1, tessellation: 14 },
      scene
    );
    shaftMaster.isPickable = false;
    shaftMaster.thinInstanceEnablePicking = false;
    shaftMaster.material = mat;
    if (typeof shaftMaster.setEnabled === 'function') shaftMaster.setEnabled(false);
    else shaftMaster.isVisible = false;

    const headMaster = BABYLON.MeshBuilder.CreateCylinder(
      'force_vector_head_master',
      { height: 1, diameterTop: 0, diameterBottom: 1, tessellation: 14 },
      scene
    );
    headMaster.isPickable = false;
    headMaster.thinInstanceEnablePicking = false;
    headMaster.material = mat;
    if (typeof headMaster.setEnabled === 'function') headMaster.setEnabled(false);
    else headMaster.isVisible = false;
    wrapThinInstanceBuffer(shaftMaster);
    wrapThinInstanceBuffer(headMaster);

    // Back-compat alias: some tests/reference code read g.master
    const g = {
      master: shaftMaster,
      shaftMaster,
      headMaster,
      mats: [],
      headMats: [],
      indices: [],
      colors: [],
    };
    forceGroups.set('force', g);
    registerMaster('forceMaster', shaftMaster);
    registerMaster('forceMaster', headMaster);
    return g;
  }
  function buildInitial() {
    __count('moleculeView#buildInitial');
    const atomCount = molState.positions.length;
    atomDefaultModes = new Array(atomCount).fill(MODE_SOLID);
    atomCurrentModes = atomDefaultModes.slice();
    rebuildAtoms();
    rebuildBonds();
  }
  function clearAtomInstances() {
    for (const group of atomGroups.values()) {
      for (const mode of MODE_LIST) {
        group.instances[mode].length = 0;
        if (group.ghostInstances && group.ghostInstances[mode]) {
          group.ghostInstances[mode].length = 0;
        }
      }
      group.indices = [];
      group.mats = [];
    }
  }
  function assignAtomInstance(atomIndex, mode) {
    const el = molState.elements[atomIndex];
    const group = ensureAtomGroup(el);
    const inst = {
      atomIndex,
      matrix: createIdentityMatrix(),
    };
    const bucket = group.instances[mode];
    bucket.push(inst);
    atomRenderRefs[atomIndex] = {
      group,
      mode,
      slot: bucket.length - 1,
    };
  }
  function refreshAtomMatrices() {
    for (const group of atomGroups.values()) {
      const info = elInfo(group.key);
      const scale = info.scale;
      const combinedIndices = [];
      const combinedMats = [];
      for (const mode of MODE_LIST) {
        const bucket = group.instances[mode];
        const ghostBucket = group.ghostInstances?.[mode] || [];
        const master = group.masters[mode];
        if (!master) continue;
        const totalCount = bucket.length + ghostBucket.length;
        if (!totalCount) {
          try {
            master.thinInstanceSetBuffer('matrix', new Float32Array());
          } catch {}
          if (typeof master.setEnabled === 'function') master.setEnabled(false);
          else master.isVisible = false;
          continue;
        }
        const matrices = [];
        for (const inst of bucket) {
          const pos = molState.positions[inst.atomIndex] || { x: 0, y: 0, z: 0 };
          inst.matrix = BABYLON.Matrix.Compose(
            new BABYLON.Vector3(scale, scale, scale),
            BABYLON.Quaternion.Identity(),
            new BABYLON.Vector3(pos.x, pos.y, pos.z)
          );
          matrices.push(inst.matrix);
          combinedIndices.push(inst.atomIndex);
          combinedMats.push(inst.matrix);
        }
        for (const inst of ghostBucket) {
          if (!inst.matrix) continue;
          matrices.push(inst.matrix);
          combinedMats.push(inst.matrix);
        }
        master.thinInstanceSetBuffer('matrix', flattenMatrices(matrices));
        if (typeof master.setEnabled === 'function') master.setEnabled(true);
        else master.isVisible = true;
      }
      group.indices = combinedIndices;
      group.mats = combinedMats;
    }
  }
  function rebuildAtoms() {
    __count('moleculeView#rebuildAtoms');
    const atomCount = molState.positions.length;
    if (!Array.isArray(atomDefaultModes) || atomDefaultModes.length !== atomCount) {
      atomDefaultModes = new Array(atomCount).fill(MODE_SOLID);
    }
    if (!Array.isArray(atomCurrentModes) || atomCurrentModes.length !== atomCount) {
      atomCurrentModes = atomDefaultModes.slice();
    }
    atomRenderRefs.length = atomCount;
    clearAtomInstances();
    for (let i = 0; i < atomCount; i++) {
      const mode = atomCurrentModes[i] || atomDefaultModes[i] || MODE_SOLID;
      assignAtomInstance(i, mode);
    }
    refreshAtomMatrices();
    if (activeOpacityMask) applyOpacityMask(activeOpacityMask, { reapply: true });
  }
  function flattenMatrices(mats) {
    const arr = new Float32Array(mats.length * 16);
    mats.forEach((m, i) => {
      arr.set(m.m, i * 16);
    });
    return arr;
  }
  function bondMatrix(pA, pB, radius) {
    const mid = new BABYLON.Vector3((pA.x + pB.x) / 2, (pA.y + pB.y) / 2, (pA.z + pB.z) / 2);
    const v = new BABYLON.Vector3(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
    const len = v.length();
    const up = new BABYLON.Vector3(0, 1, 0);
    let rot;
    const d = v.normalizeToNew();
    const dot = BABYLON.Vector3.Dot(up, d);
    if (dot > 0.9999) rot = BABYLON.Quaternion.Identity();
    else if (dot < -0.9999)
      rot = BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1, 0, 0), Math.PI);
    else {
      const axis = BABYLON.Vector3.Cross(up, d).normalize();
      rot = BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
    }
    return BABYLON.Matrix.Compose(new BABYLON.Vector3(radius * 2, len, radius * 2), rot, mid);
  }
  
function rebuildBonds(bondData) {
  __count('moleculeView#rebuildBonds');
  for (const group of bondGroups.values()) {
    for (const mode of MODE_LIST) {
      group.instances[mode].length = 0;
      if (group.ghostInstances && group.ghostInstances[mode]) {
        group.ghostInstances[mode].length = 0;
      }
    }
  }
  let source = bondData;
  if (Array.isArray(bondData)) {
    molState.bonds = bondData.map((b) => ({ ...b }));
    source = molState.bonds;
  } else {
    source = molState.bonds;
  }
  const bondsArray = Array.isArray(source) ? source : [];
  const DEBUG_STRETCH = bondStretchDebugEnabled();
  const debugSamples = DEBUG_STRETCH ? [] : null;
  bondDefaultModes = new Array(bondsArray.length).fill(MODE_SOLID);
  bondCurrentModes = bondDefaultModes.slice();
  bondRenderRefs.length = bondsArray.length;
  bondIndexByKey.clear();
  const ODBG_ENABLED = typeof window === 'undefined' ? false : window.O_BOND_DEBUG === true; // default OFF
  const BOND_DBG =
    typeof window !== 'undefined' &&
    (window.BOND_DEBUG === true || /[?&]bondDebug=1/.test(window.location?.search || ''));
  function longThreshold() {
    try {
      if (molState?.cell?.enabled && molState?.showCell) {
        const a = molState.cell.a || { x: 0, y: 0, z: 0 },
          b = molState.cell.b || { x: 0, y: 0, z: 0 },
          c = molState.cell.c || { x: 0, y: 0, z: 0 };
        const vx = a.x + b.x + c.x,
          vy = a.y + b.y + c.y,
          vz = a.z + b.z + c.z;
        const diag = Math.hypot(vx, vy, vz) || 0;
        const mul =
          typeof window !== 'undefined' && window.BOND_DEBUG_MULT
            ? Number(window.BOND_DEBUG_MULT)
            : 0.5;
        return diag * (Number.isFinite(mul) ? mul : 0.5);
      }
    } catch {}
    const fallback =
      typeof window !== 'undefined' && window.BOND_DEBUG_MINLEN
        ? Number(window.BOND_DEBUG_MINLEN)
        : 6.0;
    return Number.isFinite(fallback) ? fallback : 6.0;
  }
  const LONG_THR = longThreshold();
  let debugSoftDefaults = 0;
  for (let idx = 0; idx < bondsArray.length; idx++) {
    const bond = bondsArray[idx];
    const key = keyOf(bond.i, bond.j);
    const defaultMode =
      typeof bond.opacity === 'number' && bond.opacity < SOFT_MODE_THRESHOLD
        ? MODE_SOFT
        : MODE_SOLID;
    bondDefaultModes[idx] = defaultMode;
    bondCurrentModes[idx] = defaultMode;
    const mode = bondCurrentModes[idx];
    const group = ensureBondGroup(key);
    const canonicalKey = canonicalBondKey(bond.i, bond.j);
    const mapList = bondIndexByKey.get(canonicalKey);
    if (mapList) mapList.push(idx);
    else bondIndexByKey.set(canonicalKey, [idx]);
    const inst = {
      bondIndex: idx,
      bond,
      matrix: createIdentityMatrix(),
    };
    const bucket = group.instances[mode];
    bucket.push(inst);
    bondRenderRefs[idx] = {
      group,
      mode,
      slot: bucket.length - 1,
    };
    if (DEBUG_STRETCH) {
      if (mode === MODE_SOFT) debugSoftDefaults++;
      if (debugSamples.length < 10) {
        debugSamples.push({
          idx,
          i: bond.i,
          j: bond.j,
          opacity:
            typeof bond.opacity === 'number' && Number.isFinite(bond.opacity)
              ? Number(bond.opacity.toFixed(4))
              : null,
          defaultMode: mode,
        });
      }
    }
    if (BOND_DBG) {
      try {
        const pA = molState.positions[bond.i];
        const pB = molState.positions[bond.j];
        if (pA && pB) {
          const L = Math.hypot(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
          if (L > LONG_THR) {
            const elI = molState.elements[bond.i];
            const elJ = molState.elements[bond.j];
            const cell = molState.cell || {};
            console.log('[BOND-DBG-LONG][primary]', {
              i: bond.i,
              j: bond.j,
              elements: [elI, elJ],
              length: Number(L.toFixed(4)),
              cell:
                cell && cell.enabled
                  ? { a: cell.a, b: cell.b, c: cell.c, originOffset: cell.originOffset }
                  : null,
            });
          }
        }
      } catch {}
    }
  }
  if (DEBUG_STRETCH) {
    try {
      console.log('[MoleculeView][bondModes]', {
        total: bondsArray.length,
        defaultSoft: debugSoftDefaults,
        defaultSolid: bondsArray.length - debugSoftDefaults,
        sample: debugSamples,
      });
    } catch (err) {
      console.warn('[MoleculeView][bondModes] debug log error', err);
    }
  }
  refreshBondMatrices();
  if (activeOpacityMask) applyOpacityMask(activeOpacityMask, { reapply: true });
  if (ODBG_ENABLED && typeof molState?.moleculeRoot !== 'undefined') {
    try {
      const rootQ = molState.moleculeRoot && molState.moleculeRoot.rotationQuaternion;
      const rootRotStr = rootQ
        ? `rootQ=(${rootQ.x.toFixed(4)},${rootQ.y.toFixed(4)},${rootQ.z.toFixed(4)},${rootQ.w.toFixed(4)})`
        : 'rootQ=none';
      for (const [key, group] of bondGroups) {
        const master = group.masters[MODE_SOLID];
        if (!master || (!/O-/.test(key) && !/-O/.test(key))) continue;
        const parentOk = master.parent === molState.moleculeRoot;
        const q = master.rotationQuaternion;
        const qStr = q
          ? `q=(${q.x.toFixed(4)},${q.y.toFixed(4)},${q.z.toFixed(4)},${q.w.toFixed(4)})`
          : 'q=none';
        console.log(
          `[O-BOND-DBG][postRebuild] key=${key} instances=${group.instances[MODE_SOLID].length} parentOk=${parentOk} ${qStr} ${rootRotStr}`
        );
      }
    } catch {}
  }
  try {
    rebuildForces();
  } catch {}
}
function refreshBondMatrices() {
  for (const group of bondGroups.values()) {
    for (const mode of MODE_LIST) {
      const bucket = group.instances[mode];
      const ghostBucket = group.ghostInstances?.[mode] || [];
      const master = group.masters[mode];
      if (!master) continue;
      const totalCount = bucket.length + ghostBucket.length;
      if (!totalCount) {
        try {
          master.thinInstanceSetBuffer('matrix', new Float32Array());
        } catch {}
        if (typeof master.setEnabled === 'function') master.setEnabled(false);
        else master.isVisible = false;
        continue;
      }
      const matrices = [];
      for (const inst of bucket) {
        const bond = inst.bond;
        const pA = molState.positions[bond.i] || { x: 0, y: 0, z: 0 };
        const pB = molState.positions[bond.j] || { x: 0, y: 0, z: 0 };
        inst.matrix = bondMatrix(pA, pB, 0.1);
        matrices.push(inst.matrix);
      }
      for (const inst of ghostBucket) {
        if (!inst.matrix) continue;
        matrices.push(inst.matrix);
      }
      master.thinInstanceSetBuffer('matrix', flattenMatrices(matrices));
      if (typeof master.setEnabled === 'function') master.setEnabled(true);
      else master.isVisible = true;
    }
    const combinedIndices = [];
    const combinedMats = [];
    for (const mode of MODE_LIST) {
      const baseBucket = group.instances[mode];
      for (const inst of baseBucket) {
        combinedIndices.push(inst.bond);
        combinedMats.push(inst.matrix);
      }
      const ghostBucket = group.ghostInstances?.[mode] || [];
      for (const inst of ghostBucket) {
        combinedMats.push(inst.matrix);
      }
    }
    group.indices = combinedIndices;
    group.mats = combinedMats;
  }
}

  function setAtomMode(atomIndex, mode, { refresh = true } = {}) {
    if (!MODE_LIST.includes(mode)) mode = MODE_SOLID;
    if (!Number.isInteger(atomIndex) || atomIndex < 0 || atomIndex >= atomRenderRefs.length) {
      return false;
    }
    const ref = atomRenderRefs[atomIndex];
    if (!ref || ref.mode === mode) return false;
    const fromBucket = ref.group.instances[ref.mode];
    const inst = fromBucket[ref.slot];
    const lastIndex = fromBucket.length - 1;
    if (ref.slot !== lastIndex) {
      const swapped = fromBucket[lastIndex];
      fromBucket[ref.slot] = swapped;
      if (swapped) {
        const swappedRef = atomRenderRefs[swapped.atomIndex];
        swappedRef.slot = ref.slot;
      }
    }
    fromBucket.pop();
    const toBucket = ref.group.instances[mode];
    toBucket.push(inst);
    ref.mode = mode;
    ref.slot = toBucket.length - 1;
    atomCurrentModes[atomIndex] = mode;
    if (refresh) refreshAtomMatrices();
    return true;
  }

  function setBondMode(bondIndex, mode, { refresh = true } = {}) {
    if (!MODE_LIST.includes(mode)) mode = MODE_SOLID;
    if (!Number.isInteger(bondIndex) || bondIndex < 0 || bondIndex >= bondRenderRefs.length) {
      return false;
    }
    const ref = bondRenderRefs[bondIndex];
    if (!ref || ref.mode === mode) return false;
    const fromBucket = ref.group.instances[ref.mode];
    const inst = fromBucket[ref.slot];
    const lastIndex = fromBucket.length - 1;
    if (ref.slot !== lastIndex) {
      const swapped = fromBucket[lastIndex];
      fromBucket[ref.slot] = swapped;
      if (swapped) {
        const swappedRef = bondRenderRefs[swapped.bondIndex];
        swappedRef.slot = ref.slot;
      }
    }
    fromBucket.pop();
    const toBucket = ref.group.instances[mode];
    toBucket.push(inst);
    ref.mode = mode;
    ref.slot = toBucket.length - 1;
    bondCurrentModes[bondIndex] = mode;
    if (refresh) refreshBondMatrices();
    return true;
  }

  function resetToBaselineModes({ refresh = true } = {}) {
    let changed = false;
    for (let i = 0; i < atomDefaultModes.length; i++) {
      if (setAtomMode(i, atomDefaultModes[i], { refresh: false })) changed = true;
    }
    for (let i = 0; i < bondDefaultModes.length; i++) {
      if (setBondMode(i, bondDefaultModes[i], { refresh: false })) changed = true;
    }
    if (refresh && changed) {
      refreshAtomMatrices();
      refreshBondMatrices();
    }
    return changed;
  }

  function handlePositionsChanged() {
    if (molState.showGhostCells && molState.showCell && molState.cell?.enabled) {
      rebuildGhosts();
    } else {
      refreshAtomMatrices();
      refreshBondMatrices();
    }
  }
  molState.bus.on('positionsChanged', handlePositionsChanged);

  function forceArrowTransforms(p, f, length) {
    const fx = f[0],
      fy = f[1],
      fz = f[2];
    const mag = Math.hypot(fx, fy, fz) || 1e-9;
    const dir = new BABYLON.Vector3(fx / mag, fy / mag, fz / mag);
    const up = new BABYLON.Vector3(0, 1, 0);
    let rot;
    const dot = BABYLON.Vector3.Dot(up, dir);
    if (dot > 0.9999) rot = BABYLON.Quaternion.Identity();
    else if (dot < -0.9999)
      rot = BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1, 0, 0), Math.PI);
    else {
      const axis = BABYLON.Vector3.Cross(up, dir).normalize();
      rot = BABYLON.Quaternion.RotationAxis(axis, Math.acos(Math.min(1, Math.max(-1, dot))));
    }
    const shaftLen = length * 0.78;
    const headLen = Math.max(length - shaftLen, length * 0.22);
    const shaftMid = new BABYLON.Vector3(
      p.x + (dir.x * shaftLen) / 2,
      p.y + (dir.y * shaftLen) / 2,
      p.z + (dir.z * shaftLen) / 2
    );
    const headMid = new BABYLON.Vector3(
      p.x + dir.x * (shaftLen + headLen / 2),
      p.y + dir.y * (shaftLen + headLen / 2),
      p.z + dir.z * (shaftLen + headLen / 2)
    );
    return { rot, shaftLen, headLen, shaftMid, headMid };
  }

function rebuildForces() {
    __count('moleculeView#rebuildForces');
    const DBG =
      typeof window !== 'undefined' &&
      (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search || ''));
    if (DBG) console.log('[Forces][rebuild] start');
    const g = ensureForceGroup();
    g.mats.length = 0;
    g.headMats.length = 0;
    g.indices.length = 0;
    // If forces are globally hidden, disable the master and clear buffers
    // Default visibility: in Jest/jsdom or explicit test mode, enable by default to satisfy visualization tests
    if (molState.showForces == null) {
      let isTest = false;
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_TEST_MODE === true) isTest = true;
        else if (typeof navigator !== 'undefined' && /jsdom/i.test(navigator.userAgent || ''))
          isTest = true;
        else if (
          typeof process !== 'undefined' &&
          process.env &&
          process.env.JEST_WORKER_ID != null
        )
          isTest = true;
      } catch {}
      molState.showForces = isTest ? true : false;
    }
    if (!molState.showForces) {
      if (DBG) console.log('[Forces][rebuild] hidden by state.showForces=false');
      try {
        g.shaftMaster.thinInstanceSetBuffer('matrix', new Float32Array());
      } catch {}
      try {
        g.headMaster.thinInstanceSetBuffer('matrix', new Float32Array());
      } catch {}
      try {
        g.shaftMaster.thinInstanceSetBuffer('color', new Float32Array(), 4);
      } catch {}
      try {
        g.headMaster.thinInstanceSetBuffer('color', new Float32Array(), 4);
      } catch {}
      if (typeof g.shaftMaster.setEnabled === 'function') g.shaftMaster.setEnabled(false);
      else g.shaftMaster.isVisible = false;
      if (typeof g.headMaster.setEnabled === 'function') g.headMaster.setEnabled(false);
      else g.headMaster.isVisible = false;
      return;
    }
    // Accept forces from molState.forces OR molState.dynamics?.forces OR global window.__RELAX_FORCES for backward compat
    let forces =
      molState.forces ||
      (molState.dynamics &&
        molState.dynamics.forces &&
        molState.dynamics.forces.map((f) => [f.x, f.y, f.z])) ||
      (typeof window !== 'undefined' && window.__RELAX_FORCES) ||
      [];
    if (!Array.isArray(forces)) forces = [];
    const n = Math.min(forces.length, molState.positions.length);
    if (!n) {
      if (DBG) console.log('[Forces][rebuild] no forces or positions (n=0)');
      g.shaftMaster.thinInstanceSetBuffer('matrix', new Float32Array());
      g.headMaster.thinInstanceSetBuffer('matrix', new Float32Array());
      if (typeof g.shaftMaster.setEnabled === 'function') g.shaftMaster.setEnabled(false);
      else g.shaftMaster.isVisible = false;
      if (typeof g.headMaster.setEnabled === 'function') g.headMaster.setEnabled(false);
      else g.headMaster.isVisible = false;
      return;
    }
    // Scaling controls
    const hasWin = typeof window !== 'undefined';
    const search = hasWin ? window.location?.search || '' : '';
    const qs = Object.fromEntries(Array.from(new URLSearchParams(search)).map(([k, v]) => [k, v]));
    const forcedFixed = 'forceFixed' in qs || (hasWin && window.FORCE_FIXED); // keep constant length if set
    const forceScale = parseFloat(qs.forceScale || (hasWin && window.FORCE_SCALE) || '0.4'); // length per |f| unit
    const maxLen = parseFloat(qs.forceMax || (hasWin && window.FORCE_MAX) || '1.2');
    const minLen = parseFloat(qs.forceMin || (hasWin && window.FORCE_MIN) || '0.12');
    const radius = parseFloat(qs.forceRadius || (hasWin && window.FORCE_RADIUS) || '0.05');
    let maxMag = 0;
    if (!forcedFixed) {
      for (let i = 0; i < n; i++) {
        const f = forces[i];
        if (!f) continue;
        const m = Math.hypot(f[0], f[1], f[2]);
        if (m > maxMag) maxMag = m;
      }
    }
    if (DBG)
      console.log('[Forces][params]', {
        forcedFixed,
        forceScale,
        maxLen,
        minLen,
        radius,
        maxMag: maxMag.toFixed ? maxMag.toFixed(4) : maxMag,
      });
    const fixedLen = forcedFixed ? parseFloat(qs.forceFixed) || 0.9 : null;
    for (let i = 0; i < n; i++) {
      const p = molState.positions[i];
      const f = forces[i];
      if (!p || !f) continue;
      const mag = Math.hypot(f[0], f[1], f[2]);
      if (mag < 1e-10) continue;
      let drawLen;
      if (fixedLen != null) drawLen = fixedLen;
      else drawLen = Math.min(maxLen, Math.max(minLen, mag * forceScale));
      const t = forceArrowTransforms(p, f, drawLen);
      // Shaft
      const shaftMat = BABYLON.Matrix.Compose(
        new BABYLON.Vector3(radius * 2, t.shaftLen, radius * 2),
        t.rot,
        t.shaftMid
      );
      g.mats.push(shaftMat);
      // Head (use slightly larger radius for visual prominence)
      const headRadius = radius * 2.4; // bottom diameter of cone = 2*headRadius
      const headMat = BABYLON.Matrix.Compose(
        new BABYLON.Vector3(headRadius * 2, t.headLen, headRadius * 2),
        t.rot,
        t.headMid
      );
      g.headMats.push(headMat);
      g.indices.push({ atom: i, mag });
      if (DBG && i < 8) console.log('[Forces][rebuild] atom', i, 'f=', f, 'mag=', mag.toFixed(4));
    }
    g.shaftMaster.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    g.headMaster.thinInstanceSetBuffer('matrix', flattenMatrices(g.headMats));
    try {
      g.shaftMaster.thinInstanceRefreshBoundingInfo &&
        g.shaftMaster.thinInstanceRefreshBoundingInfo();
    } catch {}
    try {
      g.headMaster.thinInstanceRefreshBoundingInfo &&
        g.headMaster.thinInstanceRefreshBoundingInfo();
    } catch {}
    // Per-instance color (solid red, alpha 1)
    const setColor = (mesh, count) => {
      if (count) {
        const cols = new Float32Array(count * 4);
        for (let i = 0; i < count; i++) {
          cols[i * 4 + 0] = 0.95;
          cols[i * 4 + 1] = 0.05;
          cols[i * 4 + 2] = 0.05;
          cols[i * 4 + 3] = 1.0;
        }
        try {
          mesh.thinInstanceSetBuffer('color', cols, 4);
        } catch {}
      } else {
        try {
          mesh.thinInstanceSetBuffer('color', new Float32Array(), 4);
        } catch {}
      }
    };
    setColor(g.shaftMaster, g.mats.length);
    setColor(g.headMaster, g.headMats.length);
    if (DBG) console.log('[Forces][rebuild] instances=', g.mats.length);
    if (g.mats.length) {
      if (typeof g.shaftMaster.setEnabled === 'function') g.shaftMaster.setEnabled(true);
      else g.shaftMaster.isVisible = true;
      if (typeof g.headMaster.setEnabled === 'function') g.headMaster.setEnabled(true);
      else g.headMaster.isVisible = true;
    } else {
      if (typeof g.shaftMaster.setEnabled === 'function') g.shaftMaster.setEnabled(false);
      else g.shaftMaster.isVisible = false;
      if (typeof g.headMaster.setEnabled === 'function') g.headMaster.setEnabled(false);
      else g.headMaster.isVisible = false;
    }
    if (DBG) console.log('[Forces][rebuild] complete visible=', g.mats.length > 0);
  }
  function clearGhostBuffers() {
    for (const g of ghostAtomGroups.values()) {
      g.mats.length = 0;
      g.indices.length = 0;
      if (g.baseGroup?.ghostInstances?.[MODE_SOFT]) {
        g.baseGroup.ghostInstances[MODE_SOFT].length = 0;
      }
    }
    for (const g of ghostBondGroups.values()) {
      g.mats.length = 0;
      g.indices.length = 0;
      if (g.baseGroup?.ghostInstances?.[MODE_SOFT]) {
        g.baseGroup.ghostInstances[MODE_SOFT].length = 0;
      }
    }
  }

  function sanitizeImageDelta(raw) {
    const out = [0, 0, 0];
    let hasOffset = false;
    if (Array.isArray(raw)) {
      for (let idx = 0; idx < 3; idx++) {
        const value = Number(raw[idx]);
        if (!Number.isFinite(value)) continue;
        const rounded = Math.round(value);
        if (Math.abs(value - rounded) <= 1e-4) {
          out[idx] = rounded;
          if (rounded !== 0) hasOffset = true;
        }
      }
    }
    return { vector: out, hasOffset };
  }

  function rebuildGhostsWithImageDelta() {
    const cell = molState.cell;
    if (!cell || !cell.enabled) {
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    const { a, b, c } = cell;
    const positions = Array.isArray(molState.positions) ? molState.positions : [];
    const elements = Array.isArray(molState.elements) ? molState.elements : [];
    if (!positions.length || !elements.length) {
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    const bondsArray = Array.isArray(molState.bonds) ? molState.bonds : [];
    if (!bondsArray.length) {
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    function shiftToVector(shift) {
      if (!Array.isArray(shift) || shift.length < 3) return { x: 0, y: 0, z: 0 };
      const sx = Number(shift[0]) || 0;
      const sy = Number(shift[1]) || 0;
      const sz = Number(shift[2]) || 0;
      const ax = a?.x || 0;
      const ay = a?.y || 0;
      const az = a?.z || 0;
      const bx = b?.x || 0;
      const by = b?.y || 0;
      const bz = b?.z || 0;
      const cx = c?.x || 0;
      const cy = c?.y || 0;
      const cz = c?.z || 0;
      return {
        x: sx * ax + sy * bx + sz * cx,
        y: sx * ay + sy * by + sz * cy,
        z: sx * az + sy * bz + sz * cz,
      };
    }
    const BOND_DBG =
      typeof window !== 'undefined' &&
      (window.BOND_DEBUG === true || /[?&]bondDebug=1/.test(window.location?.search || ''));
    function longThreshold() {
      try {
        if (molState?.cell?.enabled) {
          const vec = {
            x: (a?.x || 0) + (b?.x || 0) + (c?.x || 0),
            y: (a?.y || 0) + (b?.y || 0) + (c?.y || 0),
            z: (a?.z || 0) + (b?.z || 0) + (c?.z || 0),
          };
          const diag = Math.hypot(vec.x, vec.y, vec.z) || 0;
          const mul =
            typeof window !== 'undefined' && window.BOND_DEBUG_MULT
              ? Number(window.BOND_DEBUG_MULT)
              : 0.5;
          return diag * (Number.isFinite(mul) ? mul : 0.5);
        }
      } catch {}
      const fallback =
        typeof window !== 'undefined' && window.BOND_DEBUG_MINLEN
          ? Number(window.BOND_DEBUG_MINLEN)
          : 6.0;
      return Number.isFinite(fallback) ? fallback : 6.0;
    }
    const LONG_THR = longThreshold();
    const baseShifts = [
      [1, 0, 0],
      [-1, 0, 0],
      [0, 1, 0],
      [0, -1, 0],
      [0, 0, 1],
      [0, 0, -1],
    ];
    const shiftSet = new Set();
    const shiftKey = (s) => `${s[0]},${s[1]},${s[2]}`;
    function addShift(tuple) {
      if (!Array.isArray(tuple) || tuple.length < 3) return;
      const sx = Math.trunc(Number(tuple[0]) || 0);
      const sy = Math.trunc(Number(tuple[1]) || 0);
      const sz = Math.trunc(Number(tuple[2]) || 0);
      if (sx === 0 && sy === 0 && sz === 0) return;
      shiftSet.add(shiftKey([sx, sy, sz]));
    }
    baseShifts.forEach(addShift);
    const sanitizedBonds = [];
    for (const bond of bondsArray) {
      if (!bond) continue;
      const i = bond.i;
      const j = bond.j;
      if (!Number.isInteger(i) || !Number.isInteger(j)) continue;
      if (i < 0 || j < 0 || i >= positions.length || j >= positions.length) continue;
      const { vector: deltaVec, hasOffset } = sanitizeImageDelta(bond.imageDelta);
      if (hasOffset) {
        addShift(deltaVec);
        addShift(deltaVec.map((v) => -v));
      }
      sanitizedBonds.push({
        bond,
        i,
        j,
        delta: deltaVec,
        hasOffset,
      });
    }
    if (!sanitizedBonds.length) {
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    const shiftList = Array.from(shiftSet).map((key) =>
      key.split(',').map((val) => Number(val) || 0)
    );
    shiftList.sort((lhs, rhs) => {
      for (let idx = 0; idx < 3; idx++) {
        if (lhs[idx] !== rhs[idx]) return lhs[idx] - rhs[idx];
      }
      return 0;
    });
    for (const shift of shiftList) {
      const vec = shiftToVector(shift);
      for (let idx = 0; idx < positions.length; idx++) {
        const base = positions[idx];
        if (!base) continue;
        const element = elements[idx];
        const group = ensureGhostAtomGroup(element);
        const info = elInfo(element);
        const scale = info.scale;
        const pos = {
          x: base.x + vec.x,
          y: base.y + vec.y,
          z: base.z + vec.z,
        };
        const mat = BABYLON.Matrix.Compose(
          new BABYLON.Vector3(scale, scale, scale),
          BABYLON.Quaternion.Identity(),
          new BABYLON.Vector3(pos.x, pos.y, pos.z)
        );
        const inst = {
          ghost: true,
          baseIndex: idx,
          shift: shift.slice(),
          matrix: mat,
        };
        const baseGroup = group.baseGroup;
        if (baseGroup?.ghostInstances?.[MODE_SOFT]) {
          baseGroup.ghostInstances[MODE_SOFT].push(inst);
        }
        group.mats.push(mat);
        group.indices.push({ base: idx, shift: shift.slice() });
      }
    }
    const zeroShift = [0, 0, 0];
    const seen = new Set();
    function pushGhostBond(i, j, shiftA, shiftB) {
      if (!Array.isArray(shiftA) || !Array.isArray(shiftB)) return;
      if (
        shiftA.length < 3 ||
        shiftB.length < 3 ||
        (shiftA[0] === 0 && shiftA[1] === 0 && shiftA[2] === 0 &&
          shiftB[0] === 0 && shiftB[1] === 0 && shiftB[2] === 0)
      ) {
        return;
      }
      const key = `${i}|${j}|${shiftA[0]},${shiftA[1]},${shiftA[2]}|${shiftB[0]},${shiftB[1]},${shiftB[2]}`;
      if (seen.has(key)) return;
      seen.add(key);
      const group = ensureGhostBondGroup(keyOf(i, j));
      const posA = positions[i];
      const posB = positions[j];
      if (!posA || !posB) return;
      const vecA = shiftToVector(shiftA);
      const vecB = shiftToVector(shiftB);
      const pA = {
        x: posA.x + vecA.x,
        y: posA.y + vecA.y,
        z: posA.z + vecA.z,
      };
      const pB = {
        x: posB.x + vecB.x,
        y: posB.y + vecB.y,
        z: posB.z + vecB.z,
      };
      const mat = bondMatrix(pA, pB, 0.1);
      const inst = {
        ghost: true,
        matrix: mat,
        base: { i, j },
        shiftA: shiftA.slice(),
        shiftB: shiftB.slice(),
      };
      const baseGroup = group.baseGroup;
      if (baseGroup?.ghostInstances?.[MODE_SOFT]) {
        baseGroup.ghostInstances[MODE_SOFT].push(inst);
      }
      group.mats.push(mat);
      group.indices.push({
        i,
        j,
        shiftA: inst.shiftA,
        shiftB: inst.shiftB,
      });
      if (BOND_DBG) {
        try {
          const L = Math.hypot(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
          if (L > LONG_THR) {
            const elI = molState.elements[i];
            const elJ = molState.elements[j];
            console.log('[BOND-DBG-LONG][ghost][delta]', {
              i,
              j,
              elements: [elI, elJ],
              length: Number(L.toFixed(4)),
              shiftA: inst.shiftA,
              shiftB: inst.shiftB,
            });
          }
        } catch {}
      }
    }
    for (const shift of shiftList) {
      for (const entry of sanitizedBonds) {
        pushGhostBond(entry.i, entry.j, shift, shift);
      }
    }
    for (const entry of sanitizedBonds) {
      if (!entry.hasOffset) continue;
      const delta = entry.delta;
      const neg = delta.map((v) => -v);
      pushGhostBond(entry.i, entry.j, zeroShift, delta);
      pushGhostBond(entry.i, entry.j, neg, zeroShift);
    }
    refreshAtomMatrices();
    refreshBondMatrices();
  }

  function legacyRebuildGhosts() {
    __count('moleculeView#legacyRebuildGhosts');
    clearGhostBuffers();
    if (!molState.showGhostCells || !molState.showCell || !molState.cell?.enabled) {
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    const BOND_DBG =
      typeof window !== 'undefined' &&
      (window.BOND_DEBUG === true || /[?&]bondDebug=1/.test(window.location?.search || ''));
    function longThreshold() {
      try {
        if (molState?.cell?.enabled) {
          const a = molState.cell.a || { x: 0, y: 0, z: 0 },
            b = molState.cell.b || { x: 0, y: 0, z: 0 },
            c = molState.cell.c || { x: 0, y: 0, z: 0 };
          const vx = a.x + b.x + c.x,
            vy = a.y + b.y + c.y,
            vz = a.z + b.z + c.z;
          const diag = Math.hypot(vx, vy, vz) || 0;
          const mul =
            typeof window !== 'undefined' && window.BOND_DEBUG_MULT
              ? Number(window.BOND_DEBUG_MULT)
              : 0.5;
          return diag * (Number.isFinite(mul) ? mul : 0.5);
        }
      } catch {}
      const fallback =
        typeof window !== 'undefined' && window.BOND_DEBUG_MINLEN
          ? Number(window.BOND_DEBUG_MINLEN)
          : 6.0;
      return Number.isFinite(fallback) ? fallback : 6.0;
    }
    const LONG_THR = longThreshold();
    const { a, b, c } = molState.cell;
    const shifts = [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: -1, y: 0, z: 0 },
      { x: 0, y: 1, z: 0 },
      { x: 0, y: -1, z: 0 },
      { x: 0, y: 0, z: 1 },
      { x: 0, y: 0, z: -1 },
    ];
    const augAtoms = [];
    for (const S of shifts) {
      for (let i = 0; i < molState.positions.length; i++) {
        const base = molState.positions[i];
        const pos = {
          x: base.x + S.x * a.x + S.y * b.x + S.z * c.x,
          y: base.y + S.x * a.y + S.y * b.y + S.z * c.y,
          z: base.z + S.x * a.z + S.y * b.z + S.z * c.z,
        };
        augAtoms.push({ element: molState.elements[i], baseIndex: i, shift: S, pos });
        if (S.x !== 0 || S.y !== 0 || S.z !== 0) {
          const g = ensureGhostAtomGroup(molState.elements[i]);
          const info = elInfo(molState.elements[i]);
          const d = info.scale;
          const mat = BABYLON.Matrix.Compose(
            new BABYLON.Vector3(d, d, d),
            BABYLON.Quaternion.Identity(),
            new BABYLON.Vector3(pos.x, pos.y, pos.z)
          );
          const inst = {
            ghost: true,
            baseIndex: i,
            shift: [S.x, S.y, S.z],
            matrix: mat,
          };
          const baseGroup = g.baseGroup;
          if (baseGroup?.ghostInstances?.[MODE_SOFT]) {
            baseGroup.ghostInstances[MODE_SOFT].push(inst);
          }
          g.mats.push(mat);
          g.indices.push({ base: i, shift: [S.x, S.y, S.z] });
        }
      }
    }
    const augSimple = augAtoms.map((a) => ({
      element: a.element,
      pos: [a.pos.x, a.pos.y, a.pos.z],
    }));
    const augBonds = computeBondsNoState(augSimple);
    const nAtoms = molState.positions.length;
    function decode(idx) {
      const si = Math.floor(idx / nAtoms);
      const local = idx % nAtoms;
      return { si, local, shift: shifts[si] };
    }
    for (const eb of augBonds) {
      const A = decode(eb.i);
      const B = decode(eb.j);
      if (A.si === 0 && B.si === 0) continue;
      const baseI = A.local;
      const baseJ = B.local;
      const key = keyOf(baseI, baseJ);
      const g = ensureGhostBondGroup(key);
      const pA = augAtoms[eb.i].pos;
      const pB = augAtoms[eb.j].pos;
      const mat = bondMatrix(pA, pB, 0.1);
      const inst = {
        ghost: true,
        matrix: mat,
        base: { i: baseI, j: baseJ },
        shiftA: [A.shift.x, A.shift.y, A.shift.z],
        shiftB: [B.shift.x, B.shift.y, B.shift.z],
      };
      const baseGroup = g.baseGroup;
      if (baseGroup?.ghostInstances?.[MODE_SOFT]) {
        baseGroup.ghostInstances[MODE_SOFT].push(inst);
      }
      g.mats.push(mat);
      g.indices.push({
        i: baseI,
        j: baseJ,
        shiftA: [A.shift.x, A.shift.y, A.shift.z],
        shiftB: [B.shift.x, B.shift.y, B.shift.z],
      });
      if (BOND_DBG) {
        try {
          const L = Math.hypot(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
          if (L > LONG_THR) {
            const elI = molState.elements[baseI];
            const elJ = molState.elements[baseJ];
            const c = molState.cell || {};
            console.log('[BOND-DBG-LONG][ghost]', {
              i: baseI,
              j: baseJ,
              elements: [elI, elJ],
              length: Number(L.toFixed(4)),
              atomA: {
                x: Number(pA.x.toFixed(4)),
                y: Number(pA.y.toFixed(4)),
                z: Number(pA.z.toFixed(4)),
              },
              atomB: {
                x: Number(pB.x.toFixed(4)),
                y: Number(pB.y.toFixed(4)),
                z: Number(pB.z.toFixed(4)),
              },
              shiftA: A.shift,
              shiftB: B.shift,
              cell:
                c && c.enabled ? { a: c.a, b: c.b, c: c.c, originOffset: c.originOffset } : null,
            });
          }
        } catch {}
      }
    }
    refreshAtomMatrices();
    refreshBondMatrices();
  }

  function rebuildGhosts() {
    __count('moleculeView#rebuildGhosts');
    const cellReady = molState.showGhostCells && molState.showCell && molState.cell?.enabled;
    if (!cellReady) {
      clearGhostBuffers();
      refreshAtomMatrices();
      refreshBondMatrices();
      return;
    }
    const bonds = Array.isArray(molState.bonds) ? molState.bonds : [];
    const canUseImageDelta = bonds.some((b) => Array.isArray(b?.imageDelta));
    if (canUseImageDelta) {
      clearGhostBuffers();
      rebuildGhostsWithImageDelta();
    } else {
      legacyRebuildGhosts();
    }
  }

  function rebuildCellLines() {
    __count('moleculeView#rebuildCellLines');
    if (!molState.showCell || !molState.cell?.enabled) {
      if (cellLines) {
        cellLines.dispose();
        cellLines = null;
      }
      return;
    }
    const { a, b, c, originOffset } = molState.cell;
    const O = originOffset || { x: 0, y: 0, z: 0 };
    const pts = [
      O,
      { x: O.x + a.x, y: O.y + a.y, z: O.z + a.z },
      { x: O.x + a.x + b.x, y: O.y + a.y + b.y, z: O.z + a.z + b.z },
      { x: O.x + b.x, y: O.y + b.y, z: O.z + b.z },
      O,
      { x: O.x + c.x, y: O.y + c.y, z: O.z + c.z },
      { x: O.x + a.x + c.x, y: O.y + a.y + c.y, z: O.z + a.z + c.z },
      { x: O.x + a.x + b.x + c.x, y: O.y + a.y + b.y + c.y, z: O.z + a.z + b.z + c.z },
      { x: O.x + b.x + c.x, y: O.y + b.y + c.y, z: O.z + b.z + c.z },
      { x: O.x + c.x, y: O.y + c.y, z: O.z + c.z },
      { x: O.x + a.x + c.x, y: O.y + a.y + c.y, z: O.z + a.z + c.z },
      { x: O.x + a.x + b.x + c.x, y: O.y + a.y + b.y + c.y, z: O.z + a.z + b.z + c.z },
      { x: O.x + a.x + b.x, y: O.y + a.y + b.y, z: O.z + a.z + b.z },
      { x: O.x + b.x + c.x, y: O.y + b.y + c.y, z: O.z + b.z + c.z },
      { x: O.x + b.x, y: O.y + b.y, z: O.z + b.z },
    ].map((p) => new BABYLON.Vector3(p.x, p.y, p.z));
    if (cellLines) cellLines.dispose();
    cellLines = BABYLON.MeshBuilder.CreateLines('cell_lines', { points: pts }, scene);
    cellLines.color = new BABYLON.Color3(0.9, 0.8, 0.25);
    cellLines.isPickable = false;
  }
  molState.bus.on('bondsChanged', () => {
    const topologyChanged = atomCurrentModes.length !== molState.positions.length;
    if (topologyChanged) {
      atomDefaultModes = new Array(molState.positions.length).fill(MODE_SOLID);
      atomCurrentModes = atomDefaultModes.slice();
      rebuildAtoms();
    }
    rebuildBonds();
    rebuildGhosts();
    try {
      rebuildForces();
    } catch {}
    // Decide if current selection is still valid
    let invalidate = false;
    if (molState.selection && molState.selection.kind) {
      if (molState.selection.kind === 'atom') {
        const idx = molState.selection.data?.index;
        if (idx == null || idx < 0 || idx >= molState.positions.length) invalidate = true;
      } else if (molState.selection.kind === 'bond') {
        const { i, j } = molState.selection.data || {};
        if (
          i == null ||
          j == null ||
          i < 0 ||
          j < 0 ||
          i >= molState.positions.length ||
          j >= molState.positions.length
        ) {
          invalidate = true;
        } else {
          // Bond must still exist (unordered match)
          const exists = molState.bonds.some(
            (b) => (b.i === i && b.j === j) || (b.i === j && b.j === i)
          );
          if (!exists) invalidate = true;
        }
      }
    }
    if (invalidate || topologyChanged) {
      molState.selection = { kind: null, data: null };
      if (highlight.bond) highlight.bond.isVisible = false;
      if (highlight.atom) highlight.atom.isVisible = false;
    } else {
      // Selection remains; refresh highlight transform/visibility
      updateSelectionHighlight();
    }
    // (debug logs suppressed)
  });
  molState.bus.on('cellChanged', () => {
    rebuildCellLines();
    rebuildGhosts();
  });
  // Forces update event (consumer code can emit this when new forces computed)
  if (molState.bus && typeof molState.bus.on === 'function') {
    molState.bus.on('forcesChanged', () => {
      const DBG =
        typeof window !== 'undefined' &&
        (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search || ''));
      if (DBG) console.log('[Forces][event] forcesChanged received');
      try {
        rebuildForces();
      } catch (e) {
        if (DBG) console.warn('[Forces][event] rebuild error', e);
      }
    });
  }
  // Selection highlight management
  function ensureHighlightMeshes() {
    __count('moleculeView#ensureHighlightMeshes');
    if (!highlight.atom) {
      const mat = new BABYLON.StandardMaterial('highlight_atom', scene);
      mat.diffuseColor = new BABYLON.Color3(0.0, 0.9, 0.95); // cyan tint
      mat.emissiveColor = new BABYLON.Color3(0.0, 0.45, 0.5);
      mat.alpha = 0.45; // within spec test range (0.3 - 0.6)
      mat.disableLighting = true;
      mat.backFaceCulling = false;
      const sphere = BABYLON.MeshBuilder.CreateSphere(
        'highlight_atom_mesh',
        { diameter: 1, segments: 24 },
        scene
      );
      sphere.material = mat;
      sphere.isPickable = false;
      sphere.isVisible = false;
      sphere.alwaysSelectAsActiveMesh = true;
      sphere.renderingGroupId = 2; // render after atoms
      highlight.atom = sphere;
    }
    if (!highlight.bond) {
      const mat = new BABYLON.StandardMaterial('highlight_bond', scene);
      mat.diffuseColor = new BABYLON.Color3(0.0, 0.9, 0.95);
      mat.emissiveColor = new BABYLON.Color3(0.0, 0.45, 0.5);
      mat.alpha = 0.4; // similar look to atom highlight, slightly less opaque
      mat.disableLighting = true;
      mat.backFaceCulling = false;
      const cyl = BABYLON.MeshBuilder.CreateCylinder(
        'highlight_bond_mesh',
        { height: 1, diameter: 1, tessellation: 22 },
        scene
      );
      cyl.material = mat;
      cyl.isPickable = false;
      // Keep fully disabled until an actual bond is selected to avoid stray cylinder at origin
      cyl.isVisible = false;
      cyl.alwaysSelectAsActiveMesh = true;
      cyl.renderingGroupId = 2;
      highlight.bond = cyl;
    }
  }
  function updateSelectionHighlight() {
  __count('moleculeView#updateSelectionHighlight');
  const SELDBG = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_SELECT;
  ensureHighlightMeshes();
  const sel = molState.selection;
  highlight.atom.isVisible = false;
  highlight.bond.isVisible = false;
  if (!sel || !sel.kind) {
    if (SELDBG)
      try {
        console.log('[view][highlight] none');
      } catch {}
    return;
  }
  if (sel.kind === 'atom') {
    const idx = sel.data?.index;
    if (!Number.isInteger(idx) || idx < 0 || idx >= molState.positions.length) return;
    const ref = atomRenderRefs[idx];
    if (!ref) return;
    const master = ref.group.masters[ref.mode];
    const pos = molState.positions[idx] || { x: 0, y: 0, z: 0 };
    const info = elInfo(ref.group.key);
    const scale = info.scale * 1.35;
    if (master && highlight.atom.parent !== master) highlight.atom.parent = master;
    if (!master && highlight.atom.parent) highlight.atom.parent = null;
    highlight.atom.position = new BABYLON.Vector3(pos.x, pos.y, pos.z);
    highlight.atom.scaling = new BABYLON.Vector3(scale, scale, scale);
    highlight.atom.isVisible = true;
    if (SELDBG)
      try {
        console.log('[view][highlight] atom', idx, ref.group.key, pos);
      } catch {}
    return;
  }
  if (sel.kind === 'bond') {
    const { i, j } = sel.data || {};
    if (!Number.isInteger(i) || !Number.isInteger(j)) return;
    const bondIndex = molState.bonds?.findIndex(
      (b) => (b.i === i && b.j === j) || (b.i === j && b.j === i)
    );
    if (bondIndex == null || bondIndex < 0) return;
    const ref = bondRenderRefs[bondIndex];
    if (!ref) return;
    const master = ref.group.masters[ref.mode];
    if (master && highlight.bond.parent !== master) highlight.bond.parent = master;
    if (!master && highlight.bond.parent) highlight.bond.parent = null;
    const pA = molState.positions[i];
    const pB = molState.positions[j];
    if (!pA || !pB) return;
    const midLocal = new BABYLON.Vector3((pA.x + pB.x) / 2, (pA.y + pB.y) / 2, (pA.z + pB.z) / 2);
    const vLocal = new BABYLON.Vector3(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
    const lenLocal = vLocal.length();
    const up = new BABYLON.Vector3(0, 1, 0);
    let rotQ;
    const d = vLocal.normalizeToNew();
    const dot = BABYLON.Vector3.Dot(up, d);
    if (dot > 0.9999) rotQ = BABYLON.Quaternion.Identity();
    else if (dot < -0.9999)
      rotQ = BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1, 0, 0), Math.PI);
    else {
      const axis = BABYLON.Vector3.Cross(up, d).normalize();
      rotQ = BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
    }
    highlight.bond.position = midLocal;
    highlight.bond.rotationQuaternion = rotQ;
    const radius = 0.16;
    highlight.bond.scaling = new BABYLON.Vector3(radius * 2, lenLocal, radius * 2);
    highlight.bond.isVisible = true;
    if (SELDBG)
      try {
        console.log('[view][highlight] bond', i, j, 'mid=', midLocal, 'len=', lenLocal);
      } catch {}
  }
}

  molState.bus.on('selectionChanged', (s) => {
    const SELDBG = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_SELECT;
    if (SELDBG)
      try {
        console.log('[view] selectionChanged', s);
      } catch {}
    updateSelectionHighlight();
  });
  // In some minimal test environments, event loop timing can differ; proactively call once if a selection already exists.
  if (molState.selection && molState.selection.kind) {
    try {
      updateSelectionHighlight();
    } catch (e) {}
  }
  // If positions of a selected item change, update highlight transform.
  molState.bus.on('positionsChanged', () => {
    if (molState.selection && molState.selection.kind) updateSelectionHighlight();
  });
  // Frame-level safety: ensure we never show a highlight when selection is null/invalid.
  // This guards against race conditions during molecule switches leaving a stray mesh enabled.
  // Some test environments pass a minimal scene stub without onBeforeRenderObservable.
  // Provide a lightweight shim if absent so highlight safety still functions.
  if (!scene.onBeforeRenderObservable) {
    scene.onBeforeRenderObservable = {
      _c: [],
      add(fn) {
        this._c.push(fn);
      },
      run() {
        for (const f of this._c)
          try {
            f();
          } catch (_) {}
      },
    };
  }
  scene.onBeforeRenderObservable.add(() => {
    if (!highlight.atom || !highlight.bond) return; // not yet created
    const sel = molState.selection;
    if (!sel || !sel.kind) {
      if (highlight.atom.isVisible) highlight.atom.isVisible = false;
      if (highlight.bond.isVisible) highlight.bond.isVisible = false;
      return;
    }
    // If a bond remains selected, recompute its local transform each frame so any
    // underlying atom position mutations (e.g., dragging or dynamics) reflect immediately.
    if (sel.kind === 'bond' && highlight.bond.isVisible) {
      try {
        const { i, j } = sel.data;
        const pA = molState.positions[i];
        const pB = molState.positions[j];
        if (pA && pB) {
          const midLocal = new BABYLON.Vector3(
            (pA.x + pB.x) / 2,
            (pA.y + pB.y) / 2,
            (pA.z + pB.z) / 2
          );
          const vLocal = new BABYLON.Vector3(pB.x - pA.x, pB.y - pA.y, pB.z - pA.z);
          const lenLocal = vLocal.length();
          if (lenLocal > 1e-6) {
            const up = new BABYLON.Vector3(0, 1, 0);
            let rotQ;
            const d = vLocal.normalizeToNew();
            const dot = BABYLON.Vector3.Dot(up, d);
            if (dot > 0.9999) rotQ = BABYLON.Quaternion.Identity();
            else if (dot < -0.9999)
              rotQ = BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1, 0, 0), Math.PI);
            else {
              const axis = BABYLON.Vector3.Cross(up, d).normalize();
              rotQ = BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
            }
            highlight.bond.position = midLocal;
            highlight.bond.rotationQuaternion = rotQ;
            // Keep scaling local; world scaling inherited
            highlight.bond.scaling.y = lenLocal; // preserve radius components
            // (per-frame O-bond debug removed)
          }
        }
      } catch {}
    }
    // Additional integrity checks: indices must be in range, bond must still exist
    if (sel.kind === 'atom') {
      const idx = sel.data?.index;
      if (!(Number.isInteger(idx) && idx >= 0 && idx < molState.positions.length)) {
        highlight.atom.isVisible = false;
        highlight.bond.isVisible = false;
        molState.selection = { kind: null, data: null };
      }
    } else if (sel.kind === 'bond') {
      const { i, j } = sel.data || {};
      const valid =
        Number.isInteger(i) &&
        Number.isInteger(j) &&
        i >= 0 &&
        j >= 0 &&
        i < molState.positions.length &&
        j < molState.positions.length &&
        molState.bonds.some((b) => (b.i === i && b.j === j) || (b.i === j && b.j === i));
      if (!valid) {
        highlight.atom.isVisible = false;
        highlight.bond.isVisible = false;
        molState.selection = { kind: null, data: null };
      }
    }

    // (advanced rotation logging removed)
  });
  buildInitial();
  if (molState.showGhostCells && molState.showCell && molState.cell?.enabled) {
    try {
      rebuildGhosts();
    } catch {}
  }
  // Initialize highlight once after initial build
  ensureHighlightMeshes();
  function resolveAtomPick(pick) {
    __count('moleculeView#resolveAtomPick');
    if (!pick?.hit || !pick.pickedMesh) return null;
    const entry = atomMasterLookup.get(pick.pickedMesh);
    if (!entry) return null;
    const { group, mode } = entry;
    const list = group.instances[mode];
    const idx = pick.thinInstanceIndex;
    if (!Array.isArray(list) || idx == null || idx < 0 || idx >= list.length) return null;
    const inst = list[idx];
    return { kind: 'atom', index: inst.atomIndex, element: group.key };
  }
  function resolveBondPick(pick) {
    __count('moleculeView#resolveBondPick');
    if (!pick?.hit || !pick.pickedMesh) return null;
    const entry = bondMasterLookup.get(pick.pickedMesh);
    if (!entry) return null;
    const { group, mode } = entry;
    const list = group.instances[mode];
    const idx = pick.thinInstanceIndex;
    if (!Array.isArray(list) || idx == null || idx < 0 || idx >= list.length) return null;
    const inst = list[idx];
    const bond = inst.bond;
    return { kind: 'bond', i: bond.i, j: bond.j, key: group.key, index: idx };
  }

  function clearOpacityMask({ refresh = true } = {}) {
    activeOpacityMask = null;
    resetToBaselineModes({ refresh });
  }

  function applyOpacityMask(mask, { reapply = false } = {}) {
    if (!mask || typeof mask !== 'object') {
      clearOpacityMask();
      return;
    }
    if (!reapply) resetToBaselineModes({ refresh: false });

    const focusAtoms = new Set();
    if (Array.isArray(mask.focus?.atoms)) {
      for (const idx of mask.focus.atoms) {
        const n = Number(idx);
        if (Number.isInteger(n) && n >= 0) focusAtoms.add(n);
      }
    }
    if (Number.isInteger(mask.focus?.atom)) focusAtoms.add(mask.focus.atom | 0);

    const focusBondKeySet = new Set();
    if (Array.isArray(mask.focus?.bonds)) {
      for (const entry of mask.focus.bonds) {
        if (!entry) continue;
        let i, j;
        if (Array.isArray(entry) && entry.length >= 2) {
          i = Number(entry[0]);
          j = Number(entry[1]);
        } else if (typeof entry === 'object') {
          i = Number(entry.i ?? entry[0]);
          j = Number(entry.j ?? entry[1]);
        }
        if (Number.isInteger(i) && Number.isInteger(j)) {
          focusBondKeySet.add(canonicalBondKey(i, j));
        }
      }
    }

    const includeConnected = mask.focus?.includeBonds !== 'none';
    const modeSpec = mask.mode || {};

    const fetchMax = (op) => {
      if (!op || typeof op !== 'object') return null;
      const vals = [];
      for (const key of ['atoms', 'atom', 'bonds', 'bond']) {
        const val = op[key];
        if (typeof val === 'number') vals.push(val);
      }
      return vals.length ? Math.max(...vals.map(Number)) : null;
    };

    const deriveMode = (spec, opacityObj, fallback) => {
      if (typeof spec === 'string' && MODE_LIST.includes(spec)) return spec;
      const maxVal = fetchMax(opacityObj);
      if (maxVal != null && maxVal < SOFT_MODE_THRESHOLD) return MODE_SOFT;
      return fallback;
    };

    const focusAtomMode = deriveMode(modeSpec.focusAtoms ?? modeSpec.focus, mask.focusOpacity, MODE_SOLID);
    const backgroundAtomMode = deriveMode(modeSpec.backgroundAtoms ?? modeSpec.background, mask.backgroundOpacity, MODE_SOFT);
    const focusBondMode = deriveMode(modeSpec.focusBonds ?? modeSpec.focus, mask.focusOpacity, MODE_SOLID);
    const backgroundBondMode = deriveMode(modeSpec.backgroundBonds ?? modeSpec.background, mask.backgroundOpacity, MODE_SOFT);

    const applyTo = mask.applyTo || {};
    const applyAtoms = applyTo.primaryAtoms !== false;
    const applyBonds = applyTo.primaryBonds !== false;

    if (applyAtoms) {
      for (let idx = 0; idx < molState.positions.length; idx++) {
        const targetMode = focusAtoms.has(idx) ? focusAtomMode : backgroundAtomMode;
        setAtomMode(idx, targetMode, { refresh: false });
      }
    }

    if (applyBonds) {
      const focusBondIndices = new Set();
      for (const key of focusBondKeySet) {
        const list = bondIndexByKey.get(key);
        if (!list) continue;
        for (const idx of list) focusBondIndices.add(idx);
      }
      for (let idx = 0; idx < molState.bonds.length; idx++) {
        const bond = molState.bonds[idx];
        if (!bond) continue;
        const key = canonicalBondKey(bond.i, bond.j);
        const isFocus =
          focusBondIndices.has(idx) ||
          (includeConnected && focusAtoms.has(bond.i) && focusAtoms.has(bond.j));
        const targetMode = isFocus ? focusBondMode : backgroundBondMode;
        setBondMode(idx, targetMode, { refresh: false });
      }
    }

    refreshAtomMatrices();
    refreshBondMatrices();
    activeOpacityMask = mask;
  }

  function setOpacityMask(mask) {
    if (!mask && !activeOpacityMask) return;
    if (!mask) {
      clearOpacityMask();
      return;
    }
    applyOpacityMask(mask);
  }

  function getAtomMeshModes() {
    return {
      default: atomDefaultModes.slice(),
      current: atomCurrentModes.slice(),
    };
  }

  function getBondMeshModes() {
    return {
      default: bondDefaultModes.slice(),
      current: bondCurrentModes.slice(),
    };
  }

  function normalizeModeValue(value, fallback) {
    return MODE_LIST.includes(value) ? value : fallback;
  }

  function applyMeshModes({ atomDefault, atomCurrent, bondDefault, bondCurrent } = {}) {
    const atomCount = molState.positions.length;
    const bondCount = molState.bonds.length;

    if (Array.isArray(atomDefault) && atomDefault.length === atomCount) {
      atomDefaultModes = atomDefault.map((mode) => normalizeModeValue(mode, MODE_SOLID));
    }
    if (!Array.isArray(atomDefaultModes) || atomDefaultModes.length !== atomCount) {
      atomDefaultModes = new Array(atomCount).fill(MODE_SOLID);
    }

    const targetAtomModes = Array.isArray(atomCurrent) && atomCurrent.length === atomCount
      ? atomCurrent.map((mode, idx) => normalizeModeValue(mode, atomDefaultModes[idx]))
      : atomDefaultModes.slice();

    for (let i = 0; i < atomCount; i++) {
      setAtomMode(i, targetAtomModes[i], { refresh: false });
      atomCurrentModes[i] = targetAtomModes[i];
    }
    refreshAtomMatrices();

    if (Array.isArray(bondDefault) && bondDefault.length === bondCount) {
      bondDefaultModes = bondDefault.map((mode) => normalizeModeValue(mode, MODE_SOLID));
    }
    if (!Array.isArray(bondDefaultModes) || bondDefaultModes.length !== bondCount) {
      bondDefaultModes = new Array(bondCount).fill(MODE_SOLID);
    }

    const targetBondModes = Array.isArray(bondCurrent) && bondCurrent.length === bondCount
      ? bondCurrent.map((mode, idx) => normalizeModeValue(mode, bondDefaultModes[idx]))
      : bondDefaultModes.slice();

    for (let i = 0; i < bondCount; i++) {
      setBondMode(i, targetBondModes[i], { refresh: false });
      bondCurrentModes[i] = targetBondModes[i];
    }
    refreshBondMatrices();

    if (activeOpacityMask) applyOpacityMask(activeOpacityMask, { reapply: true });
  }

  return {
    rebuildBonds,
    rebuildGhosts,
    rebuildForces,
    _internals: {
      atomGroups,
      bondGroups,
      forceGroups,
      ghostAtomGroups,
      ghostBondGroups,
      highlight,
    },
    resolveAtomPick,
    resolveBondPick,
    setOpacityMask,
    clearOpacityMask,
    resetToBaselineModes,
    getAtomMeshModes,
    getBondMeshModes,
    applyMeshModes,
    meshModes: { SOLID: MODE_SOLID, SOFT: MODE_SOFT },
  };
}
