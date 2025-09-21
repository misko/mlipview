// Very small element table; extend as needed.
export const ELEMENTS = {
  H:  { covRad: 0.31, color: new BABYLON.Color3(0.95, 0.97, 1.00), vdw: 1.20, scale: 0.45 },
  C:  { covRad: 0.76, color: new BABYLON.Color3(0.22, 0.24, 0.28), vdw: 1.70, scale: 0.70 },
  N:  { covRad: 0.71, color: new BABYLON.Color3(0.38, 0.56, 0.97), vdw: 1.55, scale: 0.65 },
  O:  { covRad: 0.66, color: new BABYLON.Color3(0.95, 0.33, 0.33), vdw: 1.52, scale: 0.65 },
  S:  { covRad: 1.05, color: new BABYLON.Color3(0.98, 0.85, 0.27), vdw: 1.80, scale: 0.80 },
  X:  { covRad: 0.80, color: new BABYLON.Color3(0.70, 0.70, 0.70), vdw: 1.60, scale: 0.70 },
};

export function elInfo(sym) {
  return ELEMENTS[sym] || ELEMENTS.X;
}
