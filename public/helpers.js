// Rotate the +Y axis to align with vector `dir`
export function quatYto(dir) {
  const up = BABYLON.Vector3.Up();
  const d = dir.normalizeToNew();
  const dot = BABYLON.Vector3.Dot(up, d);
  if (dot > 0.9999) return BABYLON.Quaternion.Identity();
  if (dot < -0.9999) return BABYLON.Quaternion.RotationAxis(BABYLON.Vector3.Right(), Math.PI);
  const axis = BABYLON.Vector3.Cross(up, d).normalize();
  return BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
}

// Write an array of matrices into a Float32Array and set as the thin-instance buffer.
export function setThinInstanceMatrices(mesh, matrices) {
  const count = matrices.length;
  const buffer = new Float32Array(count * 16);
  for (let i = 0; i < count; i++) {
    matrices[i].copyToArray(buffer, i * 16);
  }
  mesh.thinInstanceSetBuffer("matrix", buffer, 16, true);
  // ensure bounding info includes instances
  mesh.thinInstanceRefreshBoundingInfo(true);
}

// Get/Set a single thin-instance matrix by index (and commit)
export function getThinInstanceMatrix(mesh, index) {
  const m = new BABYLON.Matrix();
  // Babylon keeps the buffer on GPU; we track our own matrices in JS.
  // If you didn't store them, you can reconstruct from position we hold.
  // For safety here, we return Identity (callers compose their own).
  return BABYLON.Matrix.Identity(); 
}

export function setThinInstanceMatrix(mesh, index, matrix) {
  mesh.thinInstanceSetMatrixAt(index, matrix);
  mesh.thinInstanceBufferUpdated("matrix"); // commit to GPU
}
