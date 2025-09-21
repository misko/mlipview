import { quatYto, setThinInstanceMatrices } from "./helpers.js";

export function buildBenzene(
  scene,
  {
    carbonRingRadius = 1.9,
    hydrogenRadius = 3.1,
    bondRadius = 0.08,
    debugAlwaysActive = true,
    addGroundGrid = false
  } = {}
) {
  const n = 6;

  // Materials
  const matC = new BABYLON.StandardMaterial("matC", scene);
  matC.diffuseColor  = new BABYLON.Color3(0.22, 0.24, 0.28);
  matC.emissiveColor = new BABYLON.Color3(0.06, 0.06, 0.06);

  const matH = new BABYLON.StandardMaterial("matH", scene);
  matH.diffuseColor  = new BABYLON.Color3(0.95, 0.97, 1.0);
  matH.emissiveColor = new BABYLON.Color3(0.08, 0.08, 0.08);

  const matBond = new BABYLON.StandardMaterial("matBond", scene);
  matBond.diffuseColor  = new BABYLON.Color3(0.78, 0.82, 0.90);
  matBond.emissiveColor = new BABYLON.Color3(0.04, 0.04, 0.05);

  // Planar positions (mutable: we'll update when dragging)
  const carbons = [], hydrogens = [];
  for (let i = 0; i < n; i++) {
    const a = (i * 2 * Math.PI) / n;
    carbons.push(new BABYLON.Vector3(carbonRingRadius * Math.cos(a), 0, carbonRingRadius * Math.sin(a)));
    hydrogens.push(new BABYLON.Vector3(hydrogenRadius * Math.cos(a), 0, hydrogenRadius * Math.sin(a)));
  }

  // Masters at identity (keep visible so instances render)
  const baseC = BABYLON.MeshBuilder.CreateSphere("baseC", { diameter: 0.7, segments: 24 }, scene);
  baseC.material = matC;

  const baseH = BABYLON.MeshBuilder.CreateSphere("baseH", { diameter: 0.45, segments: 20 }, scene);
  baseH.material = matH;

  // Instance matrices for atoms
  setThinInstanceMatrices(
    baseC,
    carbons.map(p => BABYLON.Matrix.Compose(BABYLON.Vector3.One(), BABYLON.Quaternion.Identity(), p))
  );
  setThinInstanceMatrices(
    baseH,
    hydrogens.map(p => BABYLON.Matrix.Compose(BABYLON.Vector3.One(), BABYLON.Quaternion.Identity(), p))
  );

  // Bonds
  const bondUnit = BABYLON.MeshBuilder.CreateCylinder("bondUnit", { height: 1, diameter: 1, tessellation: 16 }, scene);
  bondUnit.material = matBond;

  function bondMatrix(aPos, bPos) {
    const mid = aPos.add(bPos).scale(0.5);
    const v = bPos.subtract(aPos);
    const len = v.length();
    const rot = quatYto(v);
    const scale = new BABYLON.Vector3(bondRadius * 2, len, bondRadius * 2);
    return BABYLON.Matrix.Compose(scale, rot, mid);
  }

  const ringPairs = [];
  for (let i = 0; i < n; i++) ringPairs.push(["C", i, "C", (i + 1) % n]);
  const chPairs = [];
  for (let i = 0; i < n; i++) chPairs.push(["C", i, "H", i]);

  const bonds = [...ringPairs, ...chPairs].map(([ta, ia, tb, ib]) => ({ a: { type: ta, index: ia }, b: { type: tb, index: ib } }));

  // Build bond matrices
  const bondMatrices = bonds.map(({ a, b }) => {
    const aPos = (a.type === "C" ? carbons[a.index] : hydrogens[a.index]);
    const bPos = (b.type === "C" ? carbons[b.index] : hydrogens[b.index]);
    return bondMatrix(aPos, bPos);
  });
  setThinInstanceMatrices(bondUnit, bondMatrices);

  // Enable picking on thin instances
  baseC.thinInstanceEnablePicking = true;
  baseH.thinInstanceEnablePicking = true;
  bondUnit.thinInstanceEnablePicking = false; // bonds not draggable

  if (debugAlwaysActive) {
    baseC.alwaysSelectAsActiveMesh = true;
    baseH.alwaysSelectAsActiveMesh = true;
    bondUnit.alwaysSelectAsActiveMesh = true;
  }

  // Expose a tiny API the interaction module can use:
  const atoms = [
    ...carbons.map((p, i) => ({ mesh: baseC, type: "C", index: i, pos: p })),
    ...hydrogens.map((p, i) => ({ mesh: baseH, type: "H", index: i, pos: p }))
  ];

  // Bond updater: call after any atom position changes
  function refreshBonds() {
    bonds.forEach((bond, bi) => {
      const aPos = (bond.a.type === "C" ? carbons[bond.a.index] : hydrogens[bond.a.index]);
      const bPos = (bond.b.type === "C" ? carbons[bond.b.index] : hydrogens[bond.b.index]);
      const m = bondMatrix(aPos, bPos);
      bondUnit.thinInstanceSetMatrixAt(bi, m);
    });
    bondUnit.thinInstanceBufferUpdated("matrix");
  }

  return { baseC, baseH, bondUnit, atoms, bonds, refreshBonds };
}
