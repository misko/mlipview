/** @jest-environment node */
// UMA MD test: short MD run (2 steps) to verify endpoint responds using UMA calculator.
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
