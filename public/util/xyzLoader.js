// Minimal XYZ parser.
// Format:
//   Line 1: integer atom count N
//   Line 2: comment (may contain metadata). Examples:
//            "Some comment Lattice=Ax,Ay,Az;Bx,By,Bz;Cx,Cy,Cz tag1=foo temperature=300"
//            Free-form text allowed; key=value pairs extracted; Lattice parsed into cell vectors.
//   Next N lines: symbol x y z (whitespace separated)
// Returns { elements, positions, comment, tags, cell, temperature }

export function parseXYZ(text) {
  // Preserve blank lines because XYZ allows an empty comment line (line 2)
  const rawLines = text.split(/\r?\n/);
  if (rawLines.length < 2) throw new Error('XYZ: not enough lines');
  const n = parseInt((rawLines[0] || '').trim(), 10);
  if (!Number.isInteger(n) || n <= 0) throw new Error('XYZ: invalid atom count');
  // We require at least N atom lines after the comment line
  if (rawLines.length < 2 + n) throw new Error('XYZ: file shorter than declared atom count');
  const comment = rawLines[1] ?? '';
  const tags = {};
  let temperature = null;
  let cell = null;
  // Extract key=value pairs (split by whitespace) but keep free text. Lattice handled specially.
  // Helper to coerce numeric tokens that may have trailing commas/semicolons
  const toNum = (v) => {
    const t = String(v).replace(/[;,]$/, '');
    const n = parseFloat(t);
    return Number.isFinite(n) ? n : NaN;
  };
  const abcParts = { a: null, b: null, c: null };
  (comment || '').split(/\s+/).forEach((tok) => {
    const eq = tok.indexOf('=');
    if (eq > 0) {
      const k = tok.slice(0, eq);
      const v = tok.slice(eq + 1);
      const key = k.toLowerCase();
      if (key === 'lattice') {
        // Expect Ax,Ay,Az;Bx,By,Bz;Cx,Cy,Cz
        const parts = v.split(';');
        if (parts.length === 3) {
          const vecs = parts.map((p) => p.split(',').map(Number));
          if (vecs.every((vv) => vv.length === 3 && vv.every(Number.isFinite))) {
            cell = {
              a: { x: vecs[0][0], y: vecs[0][1], z: vecs[0][2] },
              b: { x: vecs[1][0], y: vecs[1][1], z: vecs[1][2] },
              c: { x: vecs[2][0], y: vecs[2][1], z: vecs[2][2] },
              enabled: true,
              originOffset: { x: 0, y: 0, z: 0 },
            };
          }
        }
      } else if (key === 'abc') {
        // Allow abc=Ax,Bx,Cx (lengths only)
        const arr = v.split(',').map((s) => toNum(s));
        if (arr.length === 3 && arr.every(Number.isFinite)) {
          // store as tags; we'll resolve with angles if provided
          tags.__abc = arr;
        }
      } else if (key === 'angles' || key === 'alphabetaGamma' /*alt*/) {
        const arr = v.split(',').map((s) => toNum(s));
        if (arr.length === 3 && arr.every(Number.isFinite)) {
          tags.__angles = arr; // alpha,beta,gamma
        }
      } else if (key === 'alpha' || key === 'beta' || key === 'gamma') {
        const num = toNum(v);
        if (Number.isFinite(num)) tags[`__${key}`] = num;
      } else if (key === 'a' || key === 'b' || key === 'c') {
        // Support separated a=, b=, c= definitions often found in 'cell: a=.., b=.., c=..'
        const num = toNum(v);
        if (Number.isFinite(num)) abcParts[key] = num;
      } else if (key === 'temperature' || key === 'temp' || key === 't') {
        const t = parseFloat(v);
        if (Number.isFinite(t)) temperature = t;
      } else {
        tags[k] = v;
      }
    }
  });
  // If individual a,b,c were found, synthesize abc vector
  try {
    if (
      !tags.__abc &&
      Number.isFinite(abcParts.a) &&
      Number.isFinite(abcParts.b) &&
      Number.isFinite(abcParts.c)
    ) {
      tags.__abc = [abcParts.a, abcParts.b, abcParts.c];
    }
  } catch {}
  // If abc + angles provided, synthesize cell parameters
  try {
    if (
      !cell &&
      (tags.__abc ||
        tags.__angles ||
        tags.__alpha != null ||
        tags.__beta != null ||
        tags.__gamma != null)
    ) {
      const abc = tags.__abc || [1, 1, 1];
      const a = abc[0],
        b = abc[1],
        c = abc[2];
      const alpha = tags.__angles ? tags.__angles[0] : (tags.__alpha ?? 90);
      const beta = tags.__angles ? tags.__angles[1] : (tags.__beta ?? 90);
      const gamma = tags.__angles ? tags.__angles[2] : (tags.__gamma ?? 90);
      if ([a, b, c, alpha, beta, gamma].every(Number.isFinite)) {
        // Avoid importing buildCellFromParameters here to keep this module lightweight; inline simple build
        const ca = Math.cos((alpha * Math.PI) / 180),
          cb = Math.cos((beta * Math.PI) / 180),
          cg = Math.cos((gamma * Math.PI) / 180);
        const sg = Math.sin((gamma * Math.PI) / 180) || 1e-12;
        const ax = a,
          ay = 0,
          az = 0;
        const bx = b * cg,
          by = b * sg,
          bz = 0;
        const cx = c * cb,
          cy = c * ((ca - cb * cg) / sg),
          cz = c * Math.sqrt(Math.max(0, 1 - cb * cb - ((ca - cb * cg) / sg) ** 2));
        cell = {
          a: { x: ax, y: ay, z: az },
          b: { x: bx, y: by, z: bz },
          c: { x: cx, y: cy, z: cz },
          enabled: true,
          originOffset: { x: 0, y: 0, z: 0 },
        };
      }
    }
  } catch {}
  const elements = [];
  const positions = [];
  for (let i = 0; i < n; i++) {
    const line = rawLines[2 + i];
    if (typeof line !== 'string') throw new Error('XYZ: file shorter than declared atom count');
    const trimmed = line.trim();
    const parts = trimmed.split(/\s+/);
    if (parts.length < 4) throw new Error(`XYZ: malformed atom line ${i}`);
    const [sym, xs, ys, zs] = parts;
    const x = parseFloat(xs),
      y = parseFloat(ys),
      z = parseFloat(zs);
    if (![x, y, z].every((v) => Number.isFinite(v)))
      throw new Error(`XYZ: invalid coords line ${i}`);
    elements.push(sym);
    positions.push({ x, y, z });
  }
  return { elements, positions, comment, tags, cell, temperature };
}

// Adapter: push parsed molecule into an existing moleculeState (clears existing contents)
export function applyXYZToState(molState, parsed) {
  molState.elements = parsed.elements.slice();
  molState.positions = parsed.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
  molState.bonds = []; // topology unknown; caller should recompute if desired
  if (parsed.cell) {
    molState.cell = parsed.cell;
    molState.markCellChanged?.();
  }
  if (parsed.temperature != null) {
    // Store as default target temperature (global + state mirror)
    try {
      if (typeof window !== 'undefined')
        window.__MLIP_TARGET_TEMPERATURE = Number(parsed.temperature);
    } catch {}
    try {
      molState.dynamics = molState.dynamics || {};
      molState.dynamics.targetTemperature = Number(parsed.temperature);
    } catch {}
  }
  molState.markPositionsChanged();
  molState.markBondsChanged(); // indicates topology reset
  return molState;
}
