// public/vr/vr-utils.js
// Shared VR helpers for molecule master meshes and transforms

export function getMoleculeMasters(scene) {
  if (scene._vrMasters && scene._vrMasters.length) return scene._vrMasters;
  const masters = scene.meshes.filter(m => m && m.name && (m.name.startsWith('base_') || m.name.startsWith('bond_')));
  scene._vrMasters = masters;
  return masters;
}

export function getAnyMaster(scene) {
  const masters = getMoleculeMasters(scene);
  return masters && masters.length ? masters[0] : null;
}

export function getMoleculeDiag(scene) {
  if (scene._molDiag) return scene._molDiag;
  const masters = getMoleculeMasters(scene);
  if (!masters || masters.length === 0) { scene._molDiag = 1; return 1; }
  let min = new BABYLON.Vector3(Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY);
  let max = new BABYLON.Vector3(Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY);
  for (const m of masters) {
    try {
      m.refreshBoundingInfo && m.refreshBoundingInfo();
      const bi = m.getBoundingInfo();
      const bmin = bi.boundingBox.minimumWorld;
      const bmax = bi.boundingBox.maximumWorld;
      min = BABYLON.Vector3.Minimize(min, bmin);
      max = BABYLON.Vector3.Maximize(max, bmax);
    } catch {}
  }
  const size = max.subtract(min);
  const diag = Math.max(size.x, size.y, size.z) || 1;
  scene._molDiag = diag;
  return diag;
}

export function transformLocalToWorld(scene, pos) {
  const m = getAnyMaster(scene);
  if (!m) return pos.clone();
  const rot = m.rotationQuaternion || BABYLON.Quaternion.Identity();
  const scl = m.scaling || new BABYLON.Vector3(1,1,1);
  const uni = (Math.abs(scl.x - scl.y) < 1e-6 && Math.abs(scl.x - scl.z) < 1e-6) ? scl.x : scl.length()/Math.sqrt(3);
  const p = pos.scale(uni);
  const rotMat = new BABYLON.Matrix();
  if (BABYLON.Matrix.FromQuaternionToRef) {
    BABYLON.Matrix.FromQuaternionToRef(rot, rotMat);
  } else if (rot.toRotationMatrix) {
    rot.toRotationMatrix(rotMat);
  }
  const out = new BABYLON.Vector3();
  BABYLON.Vector3.TransformCoordinatesToRef(p, rotMat, out);
  return out.addInPlace(m.position || BABYLON.Vector3.Zero());
}
