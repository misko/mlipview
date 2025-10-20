// Validation constraints for incoming molecules (XYZ or SMILES-derived)
// Allowed periodic table symbols (expandable)
export const ALLOWED_ELEMENTS = new Set([
  'H',
  'He',
  'Li',
  'Be',
  'B',
  'C',
  'N',
  'O',
  'F',
  'Ne',
  'Na',
  'Mg',
  'Al',
  'Si',
  'P',
  'S',
  'Cl',
  'Ar',
  'K',
  'Ca',
  'Sc',
  'Ti',
  'V',
  'Cr',
  'Mn',
  'Fe',
  'Co',
  'Ni',
  'Cu',
  'Zn',
  'Ga',
  'Ge',
  'As',
  'Se',
  'Br',
  'Kr',
  'Rb',
  'Sr',
  'Y',
  'Zr',
  'Nb',
  'Mo',
  'Tc',
  'Ru',
  'Rh',
  'Pd',
  'Ag',
  'Cd',
  'In',
  'Sn',
  'Sb',
  'Te',
  'I',
  'Xe',
  'Cs',
  'Ba',
  'La',
  'Ce',
  'Pr',
  'Nd',
  'Pm',
  'Sm',
  'Eu',
  'Gd',
  'Tb',
  'Dy',
  'Ho',
  'Er',
  'Tm',
  'Yb',
  'Lu',
  'Hf',
  'Ta',
  'W',
  'Re',
  'Os',
  'Ir',
  'Pt',
  'Au',
  'Hg',
  'Tl',
  'Pb',
  'Bi',
]);

export const MAX_ATOMS = 170;

// Validate a parsed XYZ object { elements, positions }
export function validateParsedXYZ(parsed) {
  if (!parsed || !Array.isArray(parsed.elements) || !Array.isArray(parsed.positions)) {
    return { ok: false, error: 'Malformed molecule data' };
  }
  const n = parsed.elements.length;
  if (n !== parsed.positions.length) {
    return { ok: false, error: 'Mismatched elements/positions count' };
  }
  if (n <= 0) return { ok: false, error: 'Empty molecule' };
  if (n > MAX_ATOMS) return { ok: false, error: `Too many atoms (${n} > ${MAX_ATOMS})` };
  for (let i = 0; i < n; i++) {
    const sym = parsed.elements[i];
    if (!ALLOWED_ELEMENTS.has(sym)) {
      return { ok: false, error: `Element not allowed: ${sym}` };
    }
  }
  return { ok: true };
}
