// Quick sanity check replicating parseXYZ logic under Node for roy.xyz
const fs = require('fs');
const path = require('path');

function parseXYZ(text) {
  const lines = text.trim().split(/\r?\n/);
  if (lines.length < 3) throw new Error('Invalid XYZ: too few lines');
  const n = parseInt(lines[0].trim(), 10);
  if (!Number.isFinite(n) || n <= 0) throw new Error('Invalid XYZ: first line');
  if (lines.length < 2 + n)
    throw new Error(`Invalid XYZ: expected ${n} atom lines, got ${lines.length - 2}`);
  const atoms = [];
  for (let i = 2; i < 2 + n; i++) {
    const parts = lines[i].trim().split(/\s+/);
    if (parts.length < 4) throw new Error(`Invalid line ${i + 1}`);
    const [el, xs, ys, zs] = parts;
    const x = parseFloat(xs),
      y = parseFloat(ys),
      z = parseFloat(zs);
    if (![x, y, z].every(Number.isFinite)) throw new Error(`Bad coords line ${i + 1}`);
    atoms.push({ el, x, y, z });
  }
  return atoms;
}

const file = path.join(process.cwd(), 'public', 'molecules', 'roy.xyz');
const txt = fs.readFileSync(file, 'utf8');
const atoms = parseXYZ(txt);
console.log(
  'Parsed ROY atoms:',
  atoms.length,
  'first:',
  atoms[0],
  'last:',
  atoms[atoms.length - 1]
);
