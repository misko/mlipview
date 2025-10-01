/** @jest-environment node */
// UMA MD test: short MD run (2 steps) to verify endpoint responds using UMA calculator.
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

describe('UMA MD water', () => {
  test('2-step UMA MD returns energies & positions', async () => {
    if(!(await haveServer())){ console.warn('[skip] server not running'); return; }
    const body = { atomic_numbers:[8,1,1], coordinates:[[0,0,0],[0.95,0,0],[-0.24,0.93,0]], steps:2, temperature:298, timestep_fs:1.0, friction:0.02, calculator:'uma' };
    const r = await post('/md', body);
    expect(r.ok).toBe(true);
    const resp = r.json;
    expect(Array.isArray(resp.positions)).toBe(true);
    expect(resp.positions.length).toBe(3);
    expect(typeof resp.final_energy).toBe('number');
    expect(resp.calculator).toBe('uma');
  });
});
