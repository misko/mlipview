// Minimal XYZ parser.
// Format:
//   Line 1: integer atom count N
//   Line 2: comment (may contain metadata). Examples:
//            "Some comment Lattice=Ax,Ay,Az;Bx,By,Bz;Cx,Cy,Cz tag1=foo temperature=300"
//            Free-form text allowed; key=value pairs extracted; Lattice parsed into cell vectors.
//   Next N lines: symbol x y z (whitespace separated)
// Returns { elements, positions, comment, tags, cell }

export function parseXYZ(text) {
  const lines = text.split(/\r?\n/).filter(l=>l.trim().length>0);
  if (lines.length < 2) throw new Error('XYZ: not enough lines');
  const n = parseInt(lines[0].trim(),10);
  if (!Number.isInteger(n) || n <= 0) throw new Error('XYZ: invalid atom count');
  if (lines.length < 2 + n) throw new Error('XYZ: file shorter than declared atom count');
  const comment = lines[1];
  const tags = {};
  let cell = null;
  // Extract key=value pairs (split by whitespace) but keep free text. Lattice handled specially.
  comment.split(/\s+/).forEach(tok => {
    const eq = tok.indexOf('=');
    if (eq>0) {
      const k = tok.slice(0,eq);
      const v = tok.slice(eq+1);
      if (k.toLowerCase()==='lattice') {
        // Expect Ax,Ay,Az;Bx,By,Bz;Cx,Cy,Cz
        const parts = v.split(';');
        if (parts.length===3) {
          const vecs = parts.map(p=>p.split(',').map(Number));
          if (vecs.every(vv=>vv.length===3 && vv.every(Number.isFinite))) {
            cell = { a:{x:vecs[0][0],y:vecs[0][1],z:vecs[0][2]}, b:{x:vecs[1][0],y:vecs[1][1],z:vecs[1][2]}, c:{x:vecs[2][0],y:vecs[2][1],z:vecs[2][2]}, enabled:true, originOffset:{x:0,y:0,z:0} };
          }
        }
      } else {
        tags[k] = v;
      }
    }
  });
  const elements = [];
  const positions = [];
  for (let i=0;i<n;i++) {
    const parts = lines[2+i].trim().split(/\s+/);
    if (parts.length < 4) throw new Error(`XYZ: malformed atom line ${i}`);
    const [sym, xs, ys, zs] = parts;
    const x = parseFloat(xs), y=parseFloat(ys), z=parseFloat(zs);
    if (![x,y,z].every(v=>Number.isFinite(v))) throw new Error(`XYZ: invalid coords line ${i}`);
    elements.push(sym);
    positions.push({ x, y, z });
  }
  return { elements, positions, comment, tags, cell };
}

// Adapter: push parsed molecule into an existing moleculeState (clears existing contents)
export function applyXYZToState(molState, parsed) {
  molState.elements = parsed.elements.slice();
  molState.positions = parsed.positions.map(p=>({x:p.x,y:p.y,z:p.z}));
  molState.bonds = []; // topology unknown; caller should recompute if desired
  if (parsed.cell) {
    molState.cell = parsed.cell; molState.markCellChanged?.();
  }
  molState.markPositionsChanged();
  molState.markBondsChanged(); // indicates topology reset
  return molState;
}
