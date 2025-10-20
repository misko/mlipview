// Utilities to convert SMILES -> XYZ text via a public resolver (client-side only).
// We use PubChem PUG REST API for a simple mock-able path.
// Note: For offline/airgapped or to avoid external calls, tests will mock fetch responses.

function buildPubChemUrl(smiles) {
  const enc = encodeURIComponent(smiles.trim());
  // Return 3D SDF and convert to XYZ locally? PubChem can output SDF; however,
  // to keep it simple and mockable, we request SDF and use a tiny SDF->XYZ fallback.
  // PubChem SDF URL: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/<SMILES>/SDF?record_type=3d
  return `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${enc}/SDF?record_type=3d`;
}

// Parse minimal subset of SDF mol block to XYZ. This is intentionally tiny and only
// supports what we need for small organics; tests will mock a minimal SDF.
function sdfToXYZ(sdfText) {
  const lines = sdfText.split(/\r?\n/);
  // Find counts line: usually the 4th line in V2000 mol block (after 3 header lines)
  let countsIdx = -1;
  for (let i = 0; i < Math.min(20, lines.length); i++) {
    if (/V2000/.test(lines[i]) && i >= 3) {
      countsIdx = i;
      break;
    }
  }
  if (countsIdx < 0) throw new Error('SDF parse failed: counts line');
  const cnt = lines[countsIdx];
  // Positions 0-2: atom count, bond count in fixed width, but allow whitespace split fallback
  let natoms = parseInt(cnt.slice(0, 3).trim(), 10);
  if (!Number.isFinite(natoms)) {
    const parts = cnt.trim().split(/\s+/);
    natoms = parseInt(parts[0], 10);
  }
  if (!Number.isFinite(natoms) || natoms <= 0) throw new Error('SDF parse failed: atom count');
  const atomStart = countsIdx + 1;
  const atomLines = lines.slice(atomStart, atomStart + natoms);
  const atoms = [];
  for (const l of atomLines) {
    // Columns: x y z symbol ... (fixed-width but use whitespace split)
    const parts = l.trim().split(/\s+/);
    if (parts.length < 4) throw new Error('SDF parse failed: atom line');
    const x = parseFloat(parts[0]);
    const y = parseFloat(parts[1]);
    const z = parseFloat(parts[2]);
    const sym = parts[3];
    if (![x, y, z].every(Number.isFinite)) throw new Error('SDF parse failed: coords');
    atoms.push({ sym, x, y, z });
  }
  // Build XYZ text
  const out = [];
  out.push(String(atoms.length));
  out.push('Generated from SMILES via PubChem');
  for (const a of atoms) out.push(`${a.sym} ${a.x} ${a.y} ${a.z}`);
  return out.join('\n');
}

export async function smilesToXYZ(smiles) {
  const url = buildPubChemUrl(smiles);
  const res = await fetch(url);
  if (!res.ok) throw new Error(`SMILES fetch failed (${res.status})`);
  const sdf = await res.text();
  return sdfToXYZ(sdf);
}

export function isLikelySmiles(str) {
  if (!str || typeof str !== 'string') return false;
  const s = str.trim();
  if (!s) return false;
  // Very loose whitelist of characters commonly in SMILES; block spaces
  return /^[A-Za-z0-9@+\-=#()\[\]\/\\]+$/.test(s);
}

export { sdfToXYZ };
