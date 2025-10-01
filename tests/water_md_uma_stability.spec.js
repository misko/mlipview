/** @jest-environment node */
// 250-step UMA MD stability test: ensure coordinates remain finite and bounded.
import http from 'http';

function post(path, body){
  return new Promise((resolve,reject)=>{
    const data = JSON.stringify(body);
    const req = http.request('http://127.0.0.1:8000'+path,{method:'POST',headers:{'Content-Type':'application/json','Content-Length':Buffer.byteLength(data)}},res=>{
      const chunks=[];res.on('data',c=>chunks.push(c));res.on('end',()=>{
        const txt=Buffer.concat(chunks).toString('utf8');
        try{ resolve({ status:res.statusCode, ok:res.statusCode>=200&&res.statusCode<300, json: JSON.parse(txt), raw:txt }); }
        catch{ resolve({ status:res.statusCode, ok:false, raw:txt }); }
      });});
    req.on('error',reject);req.write(data);req.end();
  });
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
