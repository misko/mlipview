// molecules/default-molecules.js - Default molecule configurations
export function makeBenzeneAtoms() {
  const atoms = [];
  const n = 6, rC = 1.9, rH = 3.1;
  for (let i = 0; i < n; i++) {
    const a = (i * 2 * Math.PI) / n;
    atoms.push({ element: "C", pos: new BABYLON.Vector3(rC * Math.cos(a), 0, rC * Math.sin(a)) });
  }
  for (let i = 0; i < n; i++) {
    const a = (i * 2 * Math.PI) / n;
    atoms.push({ element: "H", pos: new BABYLON.Vector3(rH * Math.cos(a), 0, rH * Math.sin(a)) });
  }
  return atoms;
}

export const defaultBondStyles = {
  "C-C": { radius: 0.12, color: new BABYLON.Color3(0.35, 0.35, 0.38) },
  "C-H": { radius: 0.07, color: new BABYLON.Color3(0.78, 0.82, 0.90) },
  "C-N": { radius: 0.09, color: new BABYLON.Color3(0.55, 0.63, 0.86) },
  "C-O": { radius: 0.09, color: new BABYLON.Color3(0.85, 0.45, 0.45) },
  "C-S": { radius: 0.10, color: new BABYLON.Color3(0.85, 0.78, 0.35) },
};
