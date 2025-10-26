// Centralized periodic table data and helpers
// Exposes:
// - SYMBOLS: ordered list of element symbols (1..118)
// - SYMBOL_TO_Z / Z_TO_SYMBOL mappings
// - getElementInfo(symbol): returns { Z, symbol, name?, covRad, vdw, color: BABYLON.Color3, scale }
// - ELEMENTS_BY_SYMBOL: precomputed map symbol -> info
//
// Notes:
// - For visual-only needs, covRad/vdw/scale are approximations for many elements.
// - We provide specific overrides for commonly used elements to preserve existing look.
// - Elements not explicitly overridden receive reasonable defaults so rendering/bonds still work.

export const OMOL25_ELEMENTS = [
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
];

const SYMBOLS = [
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
  'Po',
  'At',
  'Rn',
  'Fr',
  'Ra',
  'Ac',
  'Th',
  'Pa',
  'U',
  'Np',
  'Pu',
  'Am',
  'Cm',
  'Bk',
  'Cf',
  'Es',
  'Fm',
  'Md',
  'No',
  'Lr',
  'Rf',
  'Db',
  'Sg',
  'Bh',
  'Hs',
  'Mt',
  'Ds',
  'Rg',
  'Cn',
  'Fl',
  'Lv',
  'Ts',
  'Og',
];

export const SYMBOL_TO_Z = Object.fromEntries(SYMBOLS.map((s, i) => [s, i + 1]));
export const Z_TO_SYMBOL = Object.fromEntries(SYMBOLS.map((s, i) => [i + 1, s]));

function defaultColorForZ(z) {
  // Deterministic pastel color from Z; avoids hard dependencies
  const h = (z * 23) % 360; // spread hues
  const s = 0.45,
    l = 0.65;
  function hslToRgb(h, s, l) {
    const c = (1 - Math.abs(2 * l - 1)) * s;
    const hp = h / 60;
    const x = c * (1 - Math.abs((hp % 2) - 1));
    let r = 0,
      g = 0,
      b = 0;
    if (hp >= 0 && hp < 1) {
      r = c;
      g = x;
    } else if (hp < 2) {
      r = x;
      g = c;
    } else if (hp < 3) {
      g = c;
      b = x;
    } else if (hp < 4) {
      g = x;
      b = c;
    } else if (hp < 5) {
      r = x;
      b = c;
    } else {
      r = c;
      b = x;
    }
    const m = l - c / 2;
    return [r + m, g + m, b + m];
  }
  const [r, g, b] = hslToRgb(h, s, l);
  try {
    return new BABYLON.Color3(r, g, b);
  } catch {
    return {
      r,
      g,
      b,
      toHexString() {
        return '#888';
      },
    };
  }
}

// Base defaults for all elements; tuned to reasonable visualization values
function baseDefaults(z, symbol) {
  const covRad = 0.5 + 0.9 * (Math.atan((z - 10) / 25) / Math.PI + 0.5); // 0.5..1.4 approx
  const vdw = 1.2 + 0.8 * (Math.atan((z - 10) / 25) / Math.PI + 0.5); // 1.2..2.0 approx
  const color = defaultColorForZ(z);
  const scale = Math.max(0.4, Math.min(0.85, vdw * 0.45));
  return { covRad: Number(covRad.toFixed(2)), vdw: Number(vdw.toFixed(2)), color, scale };
}

// Specific overrides to preserve prior appearance/values where they existed
const OVERRIDES = {
  H: { covRad: 0.31, vdw: 1.2, color: new BABYLON.Color3(0.95, 0.97, 1.0), scale: 0.45 },
  C: { covRad: 0.76, vdw: 1.7, color: new BABYLON.Color3(0.22, 0.24, 0.28), scale: 0.7 },
  N: { covRad: 0.71, vdw: 1.55, color: new BABYLON.Color3(0.38, 0.56, 0.97), scale: 0.65 },
  O: { covRad: 0.66, vdw: 1.52, color: new BABYLON.Color3(0.95, 0.33, 0.33), scale: 0.65 },
  S: { covRad: 1.05, vdw: 1.8, color: new BABYLON.Color3(0.98, 0.85, 0.27), scale: 0.8 },
  F: { covRad: 0.57, vdw: 1.47, color: new BABYLON.Color3(0.56, 0.93, 0.56), scale: 0.6 },
  Cl: { covRad: 1.02, vdw: 1.75, color: new BABYLON.Color3(0.2, 0.95, 0.2), scale: 0.72 },
  P: { covRad: 1.07, vdw: 1.8, color: new BABYLON.Color3(1.0, 0.5, 0.0), scale: 0.78 },
  Si: { covRad: 1.11, vdw: 2.1, color: new BABYLON.Color3(0.94, 0.78, 0.62), scale: 0.8 },
  Br: { covRad: 1.2, vdw: 1.85, color: new BABYLON.Color3(0.65, 0.16, 0.16), scale: 0.78 },
  I: { covRad: 1.39, vdw: 1.98, color: new BABYLON.Color3(0.58, 0.0, 0.58), scale: 0.8 },
  Sn: { covRad: 1.39, vdw: 2.17, color: new BABYLON.Color3(0.66, 0.66, 0.72), scale: 0.82 },
  Fe: { covRad: 1.26, vdw: 1.94, color: new BABYLON.Color3(0.8, 0.36, 0.36), scale: 0.78 },
  Cu: { covRad: 1.32, vdw: 1.96, color: new BABYLON.Color3(0.72, 0.45, 0.2), scale: 0.8 },
  Zn: { covRad: 1.22, vdw: 2.01, color: new BABYLON.Color3(0.49, 0.5, 0.69), scale: 0.8 },
  Ag: { covRad: 1.45, vdw: 2.1, color: new BABYLON.Color3(0.75, 0.75, 0.8), scale: 0.84 },
  Au: { covRad: 1.36, vdw: 2.16, color: new BABYLON.Color3(0.85, 0.68, 0.13), scale: 0.84 },
  X: { covRad: 0.8, vdw: 1.6, color: new BABYLON.Color3(0.7, 0.7, 0.7), scale: 0.7 },
};

export function getElementInfo(symbol) {
  const sym = SYMBOLS.includes(symbol) ? symbol : symbol === 'X' ? 'X' : null;
  if (!sym) return { Z: 0, symbol: 'X', ...OVERRIDES.X };
  const Z = SYMBOL_TO_Z[sym];
  const base = baseDefaults(Z, sym);
  const ovr = OVERRIDES[sym] || {};
  return {
    Z,
    symbol: sym,
    covRad: ovr.covRad ?? base.covRad,
    vdw: ovr.vdw ?? base.vdw,
    color: ovr.color ?? base.color,
    scale: ovr.scale ?? base.scale,
  };
}

export const ELEMENTS_BY_SYMBOL = Object.fromEntries(
  SYMBOLS.map((sym) => [sym, getElementInfo(sym)])
);

export default { SYMBOLS, SYMBOL_TO_Z, Z_TO_SYMBOL, getElementInfo, ELEMENTS_BY_SYMBOL };
