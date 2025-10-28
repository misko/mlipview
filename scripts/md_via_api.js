#!/usr/bin/env node
// Perform a /md call via the FastAPI backend and print JSON summary.
// Usage: node scripts/md_via_api.js --calc uma --steps 500 --temp 298 --coords "[[0,0,0],[0.96,0,0],[-0.24,0.93,0]]" --atomic-numbers "[8,1,1]"

function parseArg(name, def) {
  const i = process.argv.indexOf('--' + name);
  return i === -1 ? def : process.argv[i + 1];
}
const baseUrl = parseArg('url', 'http://127.0.0.1:8000');
const calc = parseArg('calc', 'uma');
const steps = parseInt(parseArg('steps', '500'), 10);
const temp = parseFloat(parseArg('temp', '298'));
const anStr = parseArg('atomic-numbers', '[8,1,1]');
const coordsStr = parseArg('coords', '[[0,0,0],[0.96,0,0],[-0.24,0.93,0]]');
let atomic_numbers, coordinates;
try {
  atomic_numbers = JSON.parse(anStr);
} catch {
  console.error('Bad atomic-numbers JSON');
  process.exit(1);
}
try {
  coordinates = JSON.parse(coordsStr);
} catch {
  console.error('Bad coords JSON');
  process.exit(1);
}

async function main() {
  const body = {
    atomic_numbers,
    coordinates,
    steps,
    temperature: temp,
    calculator: calc,
    timestep_fs: 1.0,
  };
  const resp = await fetch(baseUrl + '/md', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  if (!resp.ok) {
    const txt = await resp.text();
    console.error('HTTP error', resp.status, txt);
    process.exit(2);
  }
  const data = await resp.json();
  // Basic explode checks (large distance or NaN positions)
  let maxPair = 0,
    exploded = false;
  const pos = data.positions || [];
  for (let i = 0; i < pos.length; i++) {
    const a = pos[i];
    if (!Array.isArray(a) || a.length !== 3 || a.some((v) => !isFinite(v))) {
      exploded = true;
      break;
    }
    for (let j = i + 1; j < pos.length; j++) {
      const b = pos[j];
      const dx = a[0] - b[0],
        dy = a[1] - b[1],
        dz = a[2] - b[2];
      const d = Math.hypot(dx, dy, dz);
      if (d > maxPair) maxPair = d;
    }
  }
  if (maxPair > 10 || exploded) {
    data.warning = 'possible_explosion';
  }
  console.log(JSON.stringify(data, null, 2));
}

main().catch((e) => {
  console.error(e);
  process.exit(3);
});
