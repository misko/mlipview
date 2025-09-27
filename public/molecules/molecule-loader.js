// molecules/molecule-loader.js - Molecule loading logic
import { loadXYZFromURL } from "../loader_xyz.js";
import { buildMolecule } from "../molecule.js";
import { makeBenzeneAtoms, defaultBondStyles } from "./default-molecules.js";

// Explicitly define the canonical default molecule so BOTH desktop and VR
// pathways (which both call buildDefault via setupMolecule) stay aligned.
// If we ever change the default, updating this constant is sufficient.
export const DEFAULT_MOLECULE = 'roy.xyz';

function getRequestedMolecule() {
  try {
    const params = new URLSearchParams(window.location.search);
    let m = params.get('molecule');
    if (!m || !m.trim()) return null;
    m = m.trim().toLowerCase();
    // Allow passing with or without .xyz
    if (!m.endsWith('.xyz')) m += '.xyz';
    // Basic sanitization: no path separators
    if (m.includes('/') || m.includes('..')) return null;
    return m;
  } catch (_) {
    return null;
  }
}

export async function buildDefault(scene) {
  // Fallback chain order:
  //   1. User-requested (?molecule=)
  //   2. DEFAULT_MOLECULE (roy.xyz)
  //   3. benzene.xyz
  //   4. Procedural benzene (as a final safety)
  // This guarantees ROY loads by default in both desktop and VR modes unless
  // the user explicitly requests another molecule.
  const requested = getRequestedMolecule();
  const tried = [];
  const candidates = [];
  if (requested) candidates.push(requested);
  // Ensure DEFAULT_MOLECULE appears once
  if (!candidates.includes(DEFAULT_MOLECULE)) candidates.push(DEFAULT_MOLECULE);
  if (!candidates.includes('benzene.xyz')) candidates.push('benzene.xyz');
  // Attempt sequentially
  for (const file of candidates) {
    try {
      const { atoms } = await loadXYZFromURL(`./molecules/${file}`);
      console.log(`[loader] Loaded ${file} with`, atoms.length, 'atoms');
      // Heuristic: slightly different bond scale for benzene
      const bondScale = file === 'benzene.xyz' ? 1.30 : 1.25;
      return buildMolecule(scene, {
        atoms,
        bondScale,
        mode: 'ballstick',
        debugAlwaysActive: true,
        bondStyles: defaultBondStyles
      });
    } catch (e) {
      tried.push({ file, error: e.message });
      console.warn(`[loader] Failed loading ${file}`, e);
    }
  }
  console.warn('[loader] All file-based molecule loads failed, falling back to procedural benzene', tried);
  return buildMolecule(scene, {
    atoms: makeBenzeneAtoms(),
    bondScale: 1.30,
    mode: 'ballstick',
    debugAlwaysActive: true,
    bondStyles: {
      'C-C': defaultBondStyles['C-C'],
      'C-H': defaultBondStyles['C-H']
    }
  });
}
