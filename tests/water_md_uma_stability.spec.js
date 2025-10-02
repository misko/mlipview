/** @jest-environment node */
// 250-step UMA MD stability test: ensure coordinates remain finite and bounded.
async function post(path, body, timeoutMs=5000){
  const ac = new AbortController();
  const t = setTimeout(()=>ac.abort(), timeoutMs);
  try {
    const r = await fetch('http://127.0.0.1:8000'+path, { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body), signal: ac.signal });
    const txt = await r.text();
    let json=null; try { json = JSON.parse(txt); } catch{}
    return { status:r.status, ok:r.ok, json, raw:txt };
  } finally { clearTimeout(t); }
}

async function haveServer(){
  const r = await post('/health',{});
  return r.ok;
}

function maxPairDistance(positions){
  let maxD=0; for(let i=0;i<positions.length;i++) for(let j=i+1;j<positions.length;j++){
    const dx=positions[i][0]-positions[j][0]; const dy=positions[i][1]-positions[j][1]; const dz=positions[i][2]-positions[j][2];
    const d=Math.sqrt(dx*dx+dy*dy+dz*dz); if(d>maxD) maxD=d;
  } return maxD;
}

describe('UMA MD 250-step stability', () => {
  test('UMA 250 steps stable', async () => {
    if(!(await haveServer())){ console.warn('[skip] server not running'); return; }
    let positions=[[0,0,0],[0.95,0,0],[-0.24,0.93,0]];
    for(let step=0; step<250; step++){
      const body = { atomic_numbers:[8,1,1], coordinates:positions, steps:1, temperature:298, timestep_fs:1.0, friction:0.02, calculator:'uma' };
      const r = await post('/md', body);
      expect(r.ok).toBe(true);
      const resp = r.json;
      positions = resp.positions;
      for(const p of positions){ expect(p.length).toBe(3); for(const v of p){ expect(Number.isFinite(v)).toBe(true); } }
      const maxD = maxPairDistance(positions);
      expect(maxD).toBeLessThan(10.0);
    }
  }, 60000);
});
