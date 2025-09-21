// public/main.js
import { createBasicScene } from "./scene.js";
import { buildMolecule } from "./molecule.js";
import { enableAtomDragging } from "./interaction.js";
import { loadXYZFromURL } from "./loader_xyz.js";
import { loadBondsSpec } from "./loader_bonds.js";
import { createTorsionController } from "./torsion.js";
import { createMockMLIP } from "./physics_mock.js";
import { createForceRenderer } from "./forces.js";
import { enableBondPicking } from "./bond_pick.js";
import { createStateStore } from "./state.js";

// --- fallback benzene in case ROY.xyz is missing ---
function makeBenzeneAtoms() {
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

async function buildDefault(scene) {
  try {
    const { atoms } = await loadXYZFromURL("./molecules/roy.xyz");
    console.log("[loader] Loaded ROY.xyz with", atoms.length, "atoms");
    return buildMolecule(scene, {
      atoms,
      bondScale: 1.25,
      mode: "ballstick",
      debugAlwaysActive: true,
      bondStyles: {
        "C-C": { radius: 0.12, color: new BABYLON.Color3(0.35, 0.35, 0.38) },
        "C-H": { radius: 0.07, color: new BABYLON.Color3(0.78, 0.82, 0.90) },
        "C-N": { radius: 0.09, color: new BABYLON.Color3(0.55, 0.63, 0.86) },
        "C-O": { radius: 0.09, color: new BABYLON.Color3(0.85, 0.45, 0.45) },
        "C-S": { radius: 0.10, color: new BABYLON.Color3(0.85, 0.78, 0.35) },
      }
    });
  } catch (err) {
    console.warn("[loader] Could not load ROY.xyz, falling back to benzene.", err);
    return buildMolecule(scene, {
      atoms: makeBenzeneAtoms(),
      bondScale: 1.30,
      mode: "ballstick",
      debugAlwaysActive: true,
      bondStyles: {
        "C-C": { radius: 0.12, color: new BABYLON.Color3(0.35, 0.35, 0.38) },
        "C-H": { radius: 0.07, color: new BABYLON.Color3(0.78, 0.82, 0.90) },
      }
    });
  }
}

async function loadRotatableSpec() {
  try {
    const items = await loadBondsSpec("./molecules/ROY.BONDS");
    console.log("[torsion] Loaded rotatable spec:", items);
    return items;
  } catch (e) {
    console.warn("[torsion] No ROY.BONDS found (optional).", e);
    return [];
  }
}

// --- bootstrap ---
const canvas = document.getElementById("renderCanvas");
const engine = new BABYLON.Engine(canvas, true, { antialias: true });
const scene = createBasicScene(engine, canvas, { /* ground: false */ });

const mol = await buildDefault(scene);
const { atoms, refreshBonds } = mol;

// --- Persistent state (records torsions & can reconstruct/export)
const state = createStateStore(mol);
window.stateJson = () => state.getStateJSON();
window.stateExportXYZ = (name) => state.exportXYZ(name);

// Drag interaction
enableAtomDragging(scene, { atoms, refreshBonds, molecule: mol });

// Torsion controller (records into state)
const torsion = createTorsionController(mol, state);
window.rotateBond = (i, j, side = "j", angleDeg = 5, recompute = false) =>
  torsion.rotateAroundBond({ i, j, side, angleDeg, recompute });

// Rotatable bonds spec (optional)
const rotatable = await loadRotatableSpec();

// Bond picking + big controller-friendly UI (default: allow any bond)
enableBondPicking(scene, {
  molecule: mol,
  torsion,
  rotatableSpec: rotatable,
  strict: false
});

// --- Mock MLIP + force visualization ---
const mlip = createMockMLIP(mol); // compute() => { energy, forces }
const forceVis = createForceRenderer(scene, atoms, {
  color: new BABYLON.Color3(1, 0.2, 0.9)   // magenta-ish; length mode toggle below
});

// Energy HUD + Forces/Length toggles
const hud = document.createElement("div");
hud.style.position = "absolute";
hud.style.top = "12px";
hud.style.right = "12px";
hud.style.padding = "8px 10px";
hud.style.background = "rgba(15,18,24,0.7)";
hud.style.border = "1px solid rgba(255,255,255,0.08)";
hud.style.borderRadius = "8px";
hud.style.color = "#f5ffd1";
hud.style.font = "14px system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
hud.innerHTML = `
  <div style="font-size:12px; color:#a4b0c0; margin-bottom:4px">Mock Energy</div>
  <div id="energyVal" style="font-weight:700; font-size:22px; letter-spacing:0.4px; margin-bottom:8px;">0.000</div>
  <div style="display:flex; gap:6px; flex-wrap:wrap;">
    <button id="btnForces">Forces: OFF</button>
    <button id="btnLenMode" title="Toggle between normalized length and true magnitude">Length: Normalized</button>
  </div>
`;
document.body.appendChild(hud);
const energyVal = hud.querySelector("#energyVal");
const btnForces = hud.querySelector("#btnForces");
const btnLenMode = hud.querySelector("#btnLenMode");

let showForces = false;
btnForces.onclick = () => {
  showForces = !showForces;
  forceVis.setEnabled(showForces);
  btnForces.textContent = `Forces: ${showForces ? "ON" : "OFF"}`;
};

// Length mode toggle (normalized <-> true)
let lengthMode = "normalized";
btnLenMode.onclick = () => {
  lengthMode = (lengthMode === "normalized") ? "true" : "normalized";
  forceVis.setMode(lengthMode);
  btnLenMode.textContent = `Length: ${lengthMode === "normalized" ? "Normalized" : "True"}`;
};

// Optional: tune true-length scale at runtime from console
window.setForceScale = (s) => forceVis.setScaleTrue(s);

// Optional: tiny state bar to rebuild & export XYZ
const stateBar = document.createElement("div");
stateBar.style.position = "absolute";
stateBar.style.top = "12px";
stateBar.style.left = "50%";
stateBar.style.transform = "translateX(-50%)";
stateBar.style.background = "rgba(15,18,24,0.7)";
stateBar.style.border = "1px solid rgba(255,255,255,0.08)";
stateBar.style.borderRadius = "8px";
stateBar.style.padding = "6px 8px";
stateBar.style.display = "flex";
stateBar.style.gap = "8px";
stateBar.style.color = "#cfe3ff";
stateBar.style.font = "12px system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
stateBar.innerHTML = `
  <button id="btnRebuild">Rebuild from rotations</button>
  <button id="btnExport">Export XYZ</button>
`;
document.body.appendChild(stateBar);

stateBar.querySelector("#btnRebuild").onclick = () => {
  state.recomputeAndCommit();
};

stateBar.querySelector("#btnExport").onclick = () => {
  const xyz = state.exportXYZ("ROY_from_rotations");
  const blob = new Blob([xyz], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "reconstructed.xyz";
  a.click();
  URL.revokeObjectURL(url);
};

// Render loop: compute E,F; update energy always; update force arrows when enabled
engine.runRenderLoop(() => {
  const { energy, forces } = mlip.compute();

  energyVal.textContent = energy.toFixed(3);

  if (showForces) {
    forceVis.setForces(forces);
  }

  scene.render();
});

addEventListener("resize", () => engine.resize());
