// public/forces.js
// Draw per-atom forces as thin-instanced arrows with per-instance opacity.
// Toggle between Normalized-length and True-magnitude length.

import { setThinInstanceMatrices, quatYto } from "./helpers.js";

export function createForceRenderer(
  scene,
  atoms,
  {
    color = new BABYLON.Color3(1, 0.2, 0.9),
    // Visual shape
    shaftRadius = 0.04,
    headRadius  = 0.10,
    headLength  = 0.35,
    // Normalized mode: L = lerp(minNorm, maxNorm, mag/maxMag)
    minNormLength = 0.05,
    maxNormLength = 1.2,
    // True mode: L = clamp(scaleTrue * mag, minTrue, maxTrue)
    scaleTrue     = 0.20,   // visual scaling for |F| in "true" mode
    minTrueLength = 0.02,
    maxTrueLength = 2.0,
    // Opacity ramps with mag/maxMag in BOTH modes
    minAlpha      = 0.15,
    maxAlpha      = 1.0,
    // Start mode
    mode = "normalized"     // "normalized" | "true"
  } = {}
) {
  // Build master arrow (unit along +Y)
  const shaft = BABYLON.MeshBuilder.CreateCylinder("forceShaft",
    { height: 1, diameter: shaftRadius * 2, tessellation: 16 }, scene);
  const head = BABYLON.MeshBuilder.CreateCylinder("forceHead",
    { height: headLength, diameterTop: 0, diameterBottom: headRadius * 2, tessellation: 20 }, scene);
  shaft.position.y = 0.5 * (1 - headLength);
  head.position.y  = 1 - headLength / 2;

  const arrow = BABYLON.Mesh.MergeMeshes([shaft, head], true, true, undefined, false, true);
  arrow.name = "forceArrowMaster";
  shaft.dispose(); head.dispose();

  // Material supports per-instance vertex color with alpha
  const mat = new BABYLON.StandardMaterial("forceMat", scene);
  mat.diffuseColor = color.clone();
  mat.emissiveColor = color.scale(0.25);
  mat.specularColor = new BABYLON.Color3(0, 0, 0);
  mat.useVertexColor = true;
  mat.backFaceCulling = false;
  mat.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
  arrow.material = mat;
  arrow.isPickable = false;

  const count = atoms.length;

  // Init empty buffers
  setThinInstanceMatrices(arrow, []);
  arrow.thinInstanceSetBuffer("color", new Float32Array(0), 4, true);

  // Use shared helper for rotation from +Y to direction

  function setForces(forces) {
    // 1) magnitudes + max
    let maxMag = 0;
    const mags = new Array(count);
    for (let i = 0; i < count; i++) {
      const m = forces[i]?.length() || 0;
      mags[i] = m;
      if (m > maxMag) maxMag = m;
    }
    if (maxMag < 1e-9) maxMag = 1e-9;

    // 2) matrices + colors
    const mats = new Array(count);
    const colors = new Float32Array(count * 4);

    for (let i = 0; i < count; i++) {
      const F = forces[i] ?? BABYLON.Vector3.Zero();
      const p = atoms[i].pos;
      const mag = mags[i];

      // length
      let L;
      if (mode === "normalized") {
        const t = Math.min(1, mag / maxMag);
        L = minNormLength + (maxNormLength - minNormLength) * t;
      } else {
        L = Math.max(minTrueLength, Math.min(maxTrueLength, scaleTrue * mag));
      }

      // rotation
      let rot = BABYLON.Quaternion.Identity();
  if (mag >= 1e-6) rot = quatYto(F);

      mats[i] = BABYLON.Matrix.Compose(new BABYLON.Vector3(1, L, 1), rot, p);

      // opacity always by relative magnitude
      const tAlpha = Math.min(1, mag / maxMag);
      const a = minAlpha + (maxAlpha - minAlpha) * tAlpha;
      const o = i * 4;
      colors[o + 0] = color.r;
      colors[o + 1] = color.g;
      colors[o + 2] = color.b;
      colors[o + 3] = a;
    }

    setThinInstanceMatrices(arrow, mats);
    arrow.thinInstanceSetBuffer("color", colors, 4, true);
  }

  function setEnabled(v) { arrow.setEnabled(!!v); }
  function setMode(m) {
    mode = (m === "true") ? "true" : "normalized";
  }
  function setScaleTrue(s) {
    if (typeof s === "number" && isFinite(s) && s > 0) {
      scaleTrue = s;
    }
  }

  setEnabled(false); // start hidden

  return { setForces, setEnabled, setMode, setScaleTrue, arrow };
}
