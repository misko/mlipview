// Computes world bounding box for molecule, robust to root-only master list.
// If masters is just root, aggregates bounds from root's child meshes.
// Returns { min, max } or null if no geometry.
export function computeMoleculeWorldBounds(scene, rootOrMasters) {
  let meshes = [];
  if (Array.isArray(rootOrMasters)) {
    // If masters is just root, fallback to children
    if (rootOrMasters.length === 1 && rootOrMasters[0]?.name === 'molecule_root' && rootOrMasters[0].getChildMeshes) {
      meshes = rootOrMasters[0].getChildMeshes();
    } else {
      meshes = rootOrMasters;
    }
  } else if (rootOrMasters && rootOrMasters.name === 'molecule_root' && rootOrMasters.getChildMeshes) {
    meshes = rootOrMasters.getChildMeshes();
  } else if (rootOrMasters) {
    meshes = [rootOrMasters];
  }
  if (!meshes.length) return null;
  let min = null, max = null;
  for (const m of meshes) {
    try {
      m.refreshBoundingInfo?.();
      const bi = m.getBoundingInfo?.();
      if (!bi) continue;
      const mn = bi.boundingBox?.minimumWorld, mx = bi.boundingBox?.maximumWorld;
      if (!mn || !mx) continue;
      if (!min) { min = { x: mn.x, y: mn.y, z: mn.z }; max = { x: mx.x, y: mx.y, z: mx.z }; }
      else {
        min.x = Math.min(min.x, mn.x); min.y = Math.min(min.y, mn.y); min.z = Math.min(min.z, mn.z);
        max.x = Math.max(max.x, mx.x); max.y = Math.max(max.y, mx.y); max.z = Math.max(max.z, mx.z);
      }
    } catch {}
  }
  if (!min || !max) return null;
  return { min, max };
}
// Utility helpers for VR molecule transforms.
// ---------------------------------------------------------------------------
// DESIGN (post-refactor):
// The render layer now creates a single TransformNode named 'molecule_root'
// and parents every master thin‑instance host (atoms, bonds, forces, ghosts,
// highlights) to it. A registry is exposed at molState.__masters plus the
// root node at molState.moleculeRoot.  VR code SHOULD prefer those when
// present. Older regex-scanning master detection is kept as a backward
// compatible fallback (and for tests that construct scenes without the new
// root).  This file therefore attempts (in order):
//   1. Use explicit molecule_root if found.
//   2. Use molState.__masters (if caller passes a scene embedding it globally).
//   3. Fallback to legacy heuristic regex scan (_isMasterCandidate).
// The exported helpers keep the same shape so existing VR modules require no
// changes aside from benefiting from the more stable root anchor.
// ---------------------------------------------------------------------------

let _mastersCache = null;
let _mastersCacheSceneId = null;
let _mastersCacheVersion = 0;
/*
 * REGRESSION GUARD (late master rotation issue):
 * We intentionally collapse the master list to just 'molecule_root' when it exists so that all
 * incremental user rotations are applied once at the root. Previously rotations were pushed per
 * master mesh and any bond/atom master created later missed historical rotations, appearing static.
 * Do NOT reintroduce per-master accumulation without a backfill strategy that reapplies the entire
 * rotation history to newly created masters. Root parenting is the canonical approach.
 */


function _isMasterCandidate(m){
  if(!m || !m.name) return false;
  // Exclude highlight helpers
  if(/^highlight_/i.test(m.name)) return false;
  // We want atom_*, bond_* masters first; force_vector_master should only be included in addition, not alone
  if(/^(atom_|bond_)/i.test(m.name)) return true;
  if(/(base_|sphere_|cyl|mol|molecule|root)/i.test(m.name)) return true;
  if(/force_vector_master/i.test(m.name)) return true; // allowed, but we'll post-filter if it's the only one
  return false;
}

export function refreshMoleculeMasters(scene, { force, molState } = {}){
  if(!scene) return [];
  // Fast path: if a molecule_root exists, treat it as the single authoritative master
  try {
    // Root may be a TransformNode (not in scene.meshes) or a Mesh
    let explicitRoot = null;
    if(scene.transformNodes){
      explicitRoot = scene.transformNodes.find(t=> t && t.name === 'molecule_root');
    }
    if(!explicitRoot){
      explicitRoot = (scene.meshes||[]).find(m=> m && m.name === 'molecule_root');
    }
    if(explicitRoot){
      try { if(typeof window!=='undefined' && window.__XR_ROT_DEBUG) console.log('[XRDBG][vr-utils] explicit molecule_root detected', explicitRoot.name, 'type', explicitRoot.getClassName&&explicitRoot.getClassName()); } catch {}
      _mastersCache = [ explicitRoot ];
      _mastersCacheSceneId = scene.uid;
      _mastersCacheVersion++;
      return _mastersCache;
    }
  } catch {}
  // Alternate path: if caller supplied molState with registry
  if(molState && molState.moleculeRoot){
    _mastersCache = [ molState.moleculeRoot ];
    _mastersCacheSceneId = scene.uid;
    _mastersCacheVersion++;
    return _mastersCache;
  }
  try {
    if(force || !_mastersCache || _mastersCacheSceneId !== scene.uid){
      _mastersCache = (scene.meshes||[]).filter(_isMasterCandidate);
      // If we only captured the force master (e.g. async race where atoms/bonds not yet built), defer until next refresh
      if(_mastersCache.length===1 && /force_vector_master/i.test(_mastersCache[0].name)) {
        // Leave cache empty so next call after atoms appear repopulates
        _mastersCache = [];
      }
      if(_mastersCache.length <= 2){
        const extra = new Set(_mastersCache);
        for(const m of [..._mastersCache]){
          if(m.parent && _isMasterCandidate(m.parent)) extra.add(m.parent);
          if(m.getChildren){
            for(const ch of m.getChildren()) if(_isMasterCandidate(ch)) extra.add(ch);
          }
        }
        _mastersCache = [...extra];
      }
      _mastersCacheSceneId = scene.uid;
      _mastersCacheVersion++;
    }
  } catch { _mastersCache = []; }
  return _mastersCache;
}

export function getMoleculeMasters(scene, opts={}){
  if(!scene) return [];
  if(!_mastersCache || _mastersCacheSceneId !== scene.uid || !_mastersCache.length){
    refreshMoleculeMasters(scene, { force:true, molState: opts.molState });
  }
  return _mastersCache;
}

export function getAnyMaster(scene, opts={}){
  const m = getMoleculeMasters(scene, opts);
  return m && m.length ? m[0] : null;
}

export function getMoleculeDiag(scene, opts={}){
  const masters = getMoleculeMasters(scene, opts);
  if(!masters.length || typeof BABYLON==='undefined') return 0.0001;
  let min = new BABYLON.Vector3(Number.POSITIVE_INFINITY,Number.POSITIVE_INFINITY,Number.POSITIVE_INFINITY);
  let max = new BABYLON.Vector3(Number.NEGATIVE_INFINITY,Number.NEGATIVE_INFINITY,Number.POSITIVE_INFINITY);
  for(const m of masters){
    try {
      m.refreshBoundingInfo?.();
      const bi = m.getBoundingInfo();
      if(!bi) continue;
      min = BABYLON.Vector3.Minimize(min, bi.boundingBox.minimumWorld);
      max = BABYLON.Vector3.Maximize(max, bi.boundingBox.maximumWorld);
    } catch {}
  }
  return max.subtract(min).length();
}

// Transform a local (molecule-space) vector into world coordinates using
// first master mesh rotation, scaling, and translation plus accumulated
// quaternion if present (optional) - kept simple for now.
export function transformLocalToWorld(scene, local, opts={}){
  // Reverted to manual composition (scale -> rotate -> translate) using first master.
  if(typeof BABYLON==='undefined') return { x: local.x, y: local.y, z: local.z };
  const m = getAnyMaster(scene, opts);
  const v = new BABYLON.Vector3(local.x, local.y, local.z);
  if(m){
    const s = m.scaling || BABYLON.Vector3.One();
    v.x *= s.x; v.y *= s.y; v.z *= s.z;
    if(m.rotationQuaternion) v.applyRotationQuaternionInPlace(m.rotationQuaternion);
    if(m.position) v.addInPlace(m.position);
  }
  return v;
}


export default { getMoleculeMasters, getMoleculeDiag, transformLocalToWorld, getAnyMaster, refreshMoleculeMasters };
