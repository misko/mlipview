// molecules/molecule-loader.js - Molecule loading logic
import { loadXYZFromURL } from "../loader_xyz.js";
import { buildMolecule } from "../molecule.js";
import { makeBenzeneAtoms, defaultBondStyles } from "./default-molecules.js";

export async function buildDefault(scene) {
  try {
    const { atoms } = await loadXYZFromURL("./molecules/roy.xyz");
    console.log("[loader] Loaded ROY.xyz with", atoms.length, "atoms");
    return buildMolecule(scene, {
      atoms,
      bondScale: 1.25,
      mode: "ballstick",
      debugAlwaysActive: true,
      bondStyles: defaultBondStyles
    });
  } catch (err) {
    console.warn("[loader] Could not load ROY.xyz, falling back to benzene.", err);
    return buildMolecule(scene, {
      atoms: makeBenzeneAtoms(),
      bondScale: 1.30,
      mode: "ballstick",
      debugAlwaysActive: true,
      bondStyles: {
        "C-C": defaultBondStyles["C-C"],
        "C-H": defaultBondStyles["C-H"]
      }
    });
  }
}
