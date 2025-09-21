// Rotate +Y to align with 'dir'
export function quatYto(dir) {
  const up = BABYLON.Vector3.Up();
  const d = dir.normalizeToNew();
  const dot = BABYLON.Vector3.Dot(up, d);
  if (dot > 0.9999) return BABYLON.Quaternion.Identity();
  if (dot < -0.9999) return BABYLON.Quaternion.RotationAxis(BABYLON.Vector3.Right(), Math.PI);
  const axis = BABYLON.Vector3.Cross(up, d).normalize();
  return BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
}

// Set multiple matrices at once
export function setThinInstanceMatrices(mesh, matrices) {
  const count = matrices.length;
  const buffer = new Float32Array(count * 16);
  for (let i = 0; i < count; i++) matrices[i].copyToArray(buffer, i * 16);
  mesh.thinInstanceSetBuffer("matrix", buffer, 16, true);
  mesh.thinInstanceRefreshBoundingInfo(true);
}

// Commit a single thin-instance matrix update and refresh bounds
export function setThinInstanceMatrix(mesh, index, matrix) {
  mesh.thinInstanceSetMatrixAt(index, matrix);
  mesh.thinInstanceBufferUpdated("matrix");
  mesh.thinInstanceRefreshBoundingInfo(true);
}
