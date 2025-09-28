// Minimal XYZ parser + loader
// Format:
//   <N>
//   <comment line>
//   <element> <x> <y> <z>
//   ... (N lines)

export function parseXYZ(text) {
  const rawLines = text.split(/\r?\n/); // do not trim globally; preserve indexing for messages
  if (rawLines.length < 2) throw new Error("Invalid XYZ: missing header lines");
  const n = parseInt(rawLines[0].trim(), 10);
  if (!Number.isFinite(n) || n <= 0) {
    throw new Error("Invalid XYZ: first line must be positive integer atom count");
  }
  // Second line is nominal comment; subsequent comment (#...) or blank lines are skipped.
  const atoms = [];
  let i = 2; // start scanning potential atom lines
  while (i < rawLines.length && atoms.length < n) {
    const line = rawLines[i].trim();
    i++;
    if (!line || line.startsWith('#')) continue; // skip blank/comment
    const parts = line.split(/\s+/);
    if (parts.length < 4) {
      throw new Error(`Invalid XYZ atom line near original line ${i}`);
    }
    const [el, xs, ys, zs] = parts;
    const x = parseFloat(xs), y = parseFloat(ys), z = parseFloat(zs);
    if (![x, y, z].every(Number.isFinite)) {
      throw new Error(`Invalid XYZ coords near original line ${i}`);
    }
    atoms.push({ element: el, pos: new BABYLON.Vector3(x, y, z) });
  }
  if (atoms.length !== n) {
    throw new Error(`Invalid XYZ: expected ${n} atoms, parsed ${atoms.length}`);
  }
  return { atoms };
}

export async function loadXYZFromURL(url) {
  const res = await fetch(url, { cache: "no-cache" });
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status} ${res.statusText}`);
  const text = await res.text();
  return parseXYZ(text);
}
