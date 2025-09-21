// Minimal XYZ parser + loader
// Format:
//   <N>
//   <comment line>
//   <element> <x> <y> <z>
//   ... (N lines)

export function parseXYZ(text) {
  const lines = text.trim().split(/\r?\n/);
  if (lines.length < 3) throw new Error("Invalid XYZ: too few lines");
  const n = parseInt(lines[0].trim(), 10);
  if (!Number.isFinite(n) || n <= 0) throw new Error("Invalid XYZ: first line must be a positive integer atom count");

  if (lines.length < 2 + n) throw new Error(`Invalid XYZ: expected ${n} atom lines, got ${lines.length - 2}`);

  const atoms = [];
  for (let i = 2; i < 2 + n; i++) {
    const parts = lines[i].trim().split(/\s+/);
    if (parts.length < 4) throw new Error(`Invalid XYZ atom line at ${i+1}`);
    const el = parts[0];
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);
    if (![x,y,z].every(Number.isFinite)) throw new Error(`Invalid XYZ coords at line ${i+1}`);
    atoms.push({
      element: el,
      pos: new BABYLON.Vector3(x, y, z)
    });
  }
  return { atoms };
}

export async function loadXYZFromURL(url) {
  const res = await fetch(url, { cache: "no-cache" });
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status} ${res.statusText}`);
  const text = await res.text();
  return parseXYZ(text);
}
