// Utility helpers for VR molecule transforms.
// These functions provide cached access to molecule master meshes and
// world/local transform helpers used by the new VR interaction layer.

let _mastersCache = null;
let _mastersCacheSceneId = null;
let _mastersCacheVersion = 0;


function _isMasterCandidate(m){
  if(!m || !m.name) return false;
  // Exclude canonical highlight meshes only; VR-specific vrBondSel/vrAtomSel removed.
  if(/^highlight_/i.test(m.name)) return false;
  if(/(base_|bond_|atom_|sphere_|cyl|mol|molecule|root)/i.test(m.name)) return true;
  return false;
}

export function refreshMoleculeMasters(scene, { force } = {}){
  if(!scene) return [];
  try {
    if(force || !_mastersCache || _mastersCacheSceneId !== scene.uid){
      _mastersCache = (scene.meshes||[]).filter(_isMasterCandidate);
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

export function getMoleculeMasters(scene){
  if(!scene) return [];
  if(!_mastersCache || _mastersCacheSceneId !== scene.uid || !_mastersCache.length){
    refreshMoleculeMasters(scene, { force:true });
  }
  return _mastersCache;
}

export function getAnyMaster(scene){
  const m = getMoleculeMasters(scene);
  return m && m.length ? m[0] : null;
}

export function getMoleculeDiag(scene){
  const masters = getMoleculeMasters(scene);
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
export function transformLocalToWorld(scene, local){
  // Reverted to manual composition (scale -> rotate -> translate) using first master.
  if(typeof BABYLON==='undefined') return { x: local.x, y: local.y, z: local.z };
  const m = getAnyMaster(scene);
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
