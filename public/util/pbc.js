// PBC utilities: build simple cells and wrap atom positions into the primary cell.
// Contract:
// - computeOrthoCellFromPositions(positions, padA): returns { a,b,c, originOffset, enabled:true, synthetic:true }
//   where a,b,c are vectors spanning an axis-aligned box around all positions with a 2*padA margin (1Å default each side).
// - wrapPositionsInPlace(positions, cell): mutates positions so each lies inside the primary cell
//   defined by origin O and lattice matrix M=[a b c]. Works for general (non-orthogonal) cells.

function _mat3FromCell(cell){
  const a = cell?.a||{x:0,y:0,z:0};
  const b = cell?.b||{x:0,y:0,z:0};
  const c = cell?.c||{x:0,y:0,z:0};
  // Column-major 3x3 as flat array m = [a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z]
  return [a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z];
}

function _cross(u,v){ return { x: u.y*v.z - u.z*v.y, y: u.z*v.x - u.x*v.z, z: u.x*v.y - u.y*v.x }; }
function _dot(u,v){ return u.x*v.x + u.y*v.y + u.z*v.z; }
function _scale(u,s){ return { x:u.x*s, y:u.y*s, z:u.z*s }; }

function _fract(x){
  // Return x wrapped into [0,1)
  const f = x - Math.floor(x);
  // Guard against -0 due to FP
  return f === 1 ? 0 : f;
}

export function wrapPositionsInPlace(positions, cell){
  if (!cell || !cell.a || !cell.b || !cell.c || !positions || !positions.length) return false;
  const O = cell.originOffset || {x:0,y:0,z:0};
  // Use reciprocal vectors computed via vector triple products to avoid numerical issues and
  // to ensure correct mapping for non-orthogonal cells. Given lattice vectors a,b,c, the reciprocal
  // basis (a*, b*, c*) satisfies a*·a = 1, a*·b = 0, a*·c = 0 (and cyclic permutations).
  const a = cell.a, b = cell.b, c = cell.c;
  const v = _dot(a, _cross(b,c));
  if (!Number.isFinite(v) || Math.abs(v) < 1e-12) return false;
  const aStar = _scale(_cross(b,c), 1.0/v);
  const bStar = _scale(_cross(c,a), 1.0/v);
  const cStar = _scale(_cross(a,b), 1.0/v);
  const ax = a.x, ay = a.y, az = a.z;
  const bx = b.x, by = b.y, bz = b.z;
  const cx = c.x, cy = c.y, cz = c.z;
  for (let i=0;i<positions.length;i++){
    const p = positions[i]; if(!p) continue;
    const d = { x: p.x - O.x, y: p.y - O.y, z: p.z - O.z };
    // Fractional components along a,b,c
    let fu = _dot(aStar, d), fv = _dot(bStar, d), fw = _dot(cStar, d);
    fu = _fract(fu); fv = _fract(fv); fw = _fract(fw);
    const nx = O.x + ax*fu + bx*fv + cx*fw;
    const ny = O.y + ay*fu + by*fv + cy*fw;
    const nz = O.z + az*fu + bz*fv + cz*fw;
    p.x = nx; p.y = ny; p.z = nz;
  }
  return true;
}

export function computeOrthoCellFromPositions(positions, padA=1.0){
  const pad = Number.isFinite(padA) ? Math.max(0, padA) : 1.0;
  if (!Array.isArray(positions) || positions.length === 0) {
    // unit cube with padding centered at origin
    const L = 2*pad || 1;
    return {
      a:{x:L,y:0,z:0}, b:{x:0,y:L,z:0}, c:{x:0,y:0,z:L},
      originOffset:{x:-pad,y:-pad,z:-pad}, enabled:true, synthetic:true
    };
  }
  let minX=Infinity,minY=Infinity,minZ=Infinity,maxX=-Infinity,maxY=-Infinity,maxZ=-Infinity;
  for (const p of positions){ if(!p) continue; const x=p.x, y=p.y, z=p.z; if(x<minX)minX=x; if(y<minY)minY=y; if(z<minZ)minZ=z; if(x>maxX)maxX=x; if(y>maxY)maxY=y; if(z>maxZ)maxZ=z; }
  // Add 1 Å padding on all sides => total length increases by 2*pad
  const dx = (maxX - minX) + 2*pad;
  const dy = (maxY - minY) + 2*pad;
  const dz = (maxZ - minZ) + 2*pad;
  // Ensure strictly positive lengths (fall back to 2*pad along collapsed axes)
  const Lx = (dx>1e-8? dx : (2*pad||1));
  const Ly = (dy>1e-8? dy : (2*pad||1));
  const Lz = (dz>1e-8? dz : (2*pad||1));
  const O = { x: minX - pad, y: minY - pad, z: minZ - pad };
  return {
    a:{x:Lx,y:0,z:0}, b:{x:0,y:Ly,z:0}, c:{x:0,y:0,z:Lz},
    originOffset: O,
    enabled: true,
    synthetic: true
  };
}

export function cellToMatrixArray(cell){
  if(!cell || !cell.a || !cell.b || !cell.c) return null;
  return [
    [cell.a.x, cell.a.y, cell.a.z],
    [cell.b.x, cell.b.y, cell.b.z],
    [cell.c.x, cell.c.y, cell.c.z]
  ];
}

export default { wrapPositionsInPlace, computeOrthoCellFromPositions, cellToMatrixArray };
