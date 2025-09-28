import { computeBondsNoState } from '../bond_render.js';

export function createBondService(molState) {
  function baseAtomArray() {
    return molState.positions.map((p, i) => ({ element: molState.elements[i], pos:[p.x,p.y,p.z] }));
  }
  function hasCell() {
    const c = molState.cell; return !!(c && c.enabled && c.a && c.b && c.c);
  }
  function computePeriodicBonds() {
    // Non-periodic path: preserve full opacity shaping from computeBondsNoState (legacy smooth transparency)
    if (!hasCell()) return computeBondsNoState(baseAtomArray()).map(b=>({ i:b.i,j:b.j,length:b.length,opacity:b.opacity ?? 1 }));
    const { a, b, c } = molState.cell;
    const shifts = [
      {x:0,y:0,z:0}, a, b, c,
      {x:-a.x,y:-a.y,z:-a.z}, {x:-b.x,y:-b.y,z:-b.z}, {x:-c.x,y:-c.y,z:-c.z}
    ];
    const aug = [];
    for (let si=0; si<shifts.length; si++) {
      const S = shifts[si];
      for (let i=0;i<molState.positions.length;i++) {
        const p = molState.positions[i];
        aug.push({ element: molState.elements[i], pos:[p.x+S.x,p.y+S.y,p.z+S.z], baseIndex:i, shiftIndex:si });
      }
    }
    const augBonds = computeBondsNoState(aug.map(a=>({ element:a.element, pos:a.pos })));
    const av=a,bv=b,cv=c;
    const M=[av.x,bv.x,cv.x,av.y,bv.y,cv.y,av.z,bv.z,cv.z];
    const det = (
      M[0]*(M[4]*M[8]-M[5]*M[7]) -
      M[1]*(M[3]*M[8]-M[5]*M[6]) +
      M[2]*(M[3]*M[7]-M[4]*M[6])
    );
    let inv=null;
    if (Math.abs(det)>1e-14) {
      const invDet=1/det;
      inv=[
        (M[4]*M[8]-M[5]*M[7])*invDet,
        (M[2]*M[7]-M[1]*M[8])*invDet,
        (M[1]*M[5]-M[2]*M[4])*invDet,
        (M[5]*M[6]-M[3]*M[8])*invDet,
        (M[0]*M[8]-M[2]*M[6])*invDet,
        (M[2]*M[3]-M[0]*M[5])*invDet,
        (M[3]*M[7]-M[4]*M[6])*invDet,
        (M[1]*M[6]-M[0]*M[7])*invDet,
        (M[0]*M[4]-M[1]*M[3])*invDet
      ];
    }
    function frac(p){ if(!inv) return {u:p.x,v:p.y,w:p.z}; return { u:inv[0]*p.x+inv[1]*p.y+inv[2]*p.z, v:inv[3]*p.x+inv[4]*p.y+inv[5]*p.z, w:inv[6]*p.x+inv[7]*p.y+inv[8]*p.z }; }
    function wrap(f){ return { u:f.u-Math.floor(f.u), v:f.v-Math.floor(f.v), w:f.w-Math.floor(f.w) }; }
    function cart(f){ return { x: av.x*f.u + bv.x*f.v + cv.x*f.w, y: av.y*f.u + bv.y*f.v + cv.y*f.w, z: av.z*f.u + bv.z*f.v + cv.z*f.w }; }
    const seen=new Set();
    const out=[];
    for (const eb of augBonds) {
      const A=aug[eb.i], B=aug[eb.j];
      if (A.shiftIndex===B.shiftIndex && A.shiftIndex!==0) continue;
      const pA={x:A.pos[0],y:A.pos[1],z:A.pos[2]};
      const pB={x:B.pos[0],y:B.pos[1],z:B.pos[2]};
      const fA=wrap(frac(pA)); const fB=wrap(frac(pB));
      const cA=cart(fA); const cB=cart(fB);
      const dx=cA.x-cB.x, dy=cA.y-cB.y, dz=cA.z-cB.z; const dist=Math.sqrt(dx*dx+dy*dy+dz*dz);
      const i=Math.min(A.baseIndex,B.baseIndex); const j=Math.max(A.baseIndex,B.baseIndex);
      let du=fB.u-fA.u,dv=fB.v-fA.v,dw=fB.w-fA.w; du-=Math.round(du); dv-=Math.round(dv); dw-=Math.round(dw);
      const key = i+'_'+j+':'+du.toFixed(4)+','+dv.toFixed(4)+','+dw.toFixed(4);
      if (seen.has(key)) continue; seen.add(key);
      const fAraw=frac(pA), fBraw=frac(pB);
      const ndu=(fBraw.u-fAraw.u)-du, ndv=(fBraw.v-fAraw.v)-dv, ndw=(fBraw.w-fAraw.w)-dw;
      const crossing = (Math.abs(ndu)>1e-6 || Math.abs(ndv)>1e-6 || Math.abs(ndw)>1e-6);
      out.push({ i, j, length:dist, opacity: crossing?0.5:1.0, crossing });
    }
    // Fallback: if periodic expansion produced no bonds, fall back to non-periodic with shaped opacity
    if (!out.length) {
      return computeBondsNoState(baseAtomArray()).map(b=>({ i:b.i,j:b.j,length:b.length,opacity:b.opacity ?? 1 }));
    }
    return out;
  }
  function recomputeAndStore() {
    const bonds = computePeriodicBonds();
    // Store opacity for renderer so it can apply group alpha (harmless to existing logic using only i,j)
    molState.bonds = bonds.map(b=>({ i:b.i, j:b.j, opacity:b.opacity }));
    molState.markBondsChanged();
    return bonds;
  }
  return { computePeriodicBonds, recomputeAndStore };
}
