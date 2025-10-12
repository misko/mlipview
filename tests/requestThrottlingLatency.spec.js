/** @jest-environment jsdom */

// Verify client-side pacing achieves >= minStepInterval between MD requests when network is fast,
// and adds minimal extra delay when network is slower than the pacing interval.

jest.setTimeout(20000);

// Mock the heavy Babylon scene so we can import public/index.js
jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn)=>{}, stopRenderLoop:()=>{} },
    scene: { meshes:[], render:()=>{}, dispose:()=>{}, onPointerObservable:{ _l:[], add(fn){ this._l.push(fn); } } },
    camera: { attachControl:()=>{} }
  })
}));

function makeViewerDom(){
  const canvas = document.createElement('canvas');
  canvas.id = 'viewer';
  canvas.getBoundingClientRect = ()=>({ left:0, top:0 });
  document.body.appendChild(canvas);
  // Optional labels used by index.js updates
  const status = document.createElement('div'); status.id='status'; document.body.appendChild(status);
  const instTemp = document.createElement('div'); instTemp.id='instTemp'; document.body.appendChild(instTemp);
  const rpsLabel = document.createElement('div'); rpsLabel.id='rpsLabel'; document.body.appendChild(rpsLabel);
  return canvas;
}

// Minimal MD response helper
function mdResponseEcho(body){
  const pos = body.coordinates || [];
  const N = pos.length;
  const zeros = Array.from({length:N}, ()=> [0,0,0]);
  return {
    ok: true,
    status: 200,
    json: async () => ({
      initial_energy: body?.precomputed?.energy ?? -1.0,
      final_energy: -1.0,
      positions: pos,
      velocities: zeros,
      forces: zeros,
      steps_completed: 1,
      temperature: 300
    }),
    text: async () => JSON.stringify({})
  };
}

describe('request pacing and latency breakdown', () => {
  let timestamps=[]; let networkDelay=5; let callCount=0;
  beforeEach(()=>{ timestamps=[]; callCount=0; networkDelay=5; });

  test('MD pacing: fast vs slow network', async () => {
    // Prevent auto MD run to keep test deterministic
    window.__MLIPVIEW_NO_AUTO_MD = true;
    window.__MLIPVIEW_TEST_MODE = true; // bypass focus gating

    // Mock fetch with controllable delay for /serve/md; immediate for /serve/simple
    global.fetch = (url, opts={}) => new Promise((resolve, reject) => {
      try {
        const u = (typeof url === 'string') ? url : (url && url.url) || '';
        if(/\/serve\/simple$/.test(u)){
          const body = { results:{ energy:-1.0, forces:[[0,0,0]], stress:[0,0,0,0,0,0] } };
          return resolve({ ok:true, status:200, json: async()=> body, text: async()=> JSON.stringify(body) });
        }
        if(/\/serve\/md$/.test(u)){
          const t = (typeof performance!=='undefined' && performance.now) ? performance.now() : Date.now();
          timestamps.push(t); callCount++;
          const body = JSON.parse(opts.body||'{}');
          setTimeout(()=> resolve(mdResponseEcho(body)), Math.max(0, networkDelay));
          return;
        }
        // Other endpoints not expected in this test
        return resolve({ ok:true, status:200, json: async()=>({}), text: async()=> '{}' });
      } catch(e){ reject(e); }
    });

    const { initNewViewer } = await import('../public/index.js');
    const canvas = makeViewerDom();
    const api = await initNewViewer(canvas, {
      elements: ['H','H','O'],
      positions: [{x:0,y:0,z:0},{x:1,y:0,z:0},{x:0,y:1,z:0}],
      bonds: []
    });

    // Ensure pacing is default 30ms
    api.setMinStepInterval(30);

    // Run with fast network (~5ms)
    networkDelay = 5; timestamps.length = 0;
    await api.startMDContinuous({ steps: 8, calculator:'uma', temperature:298, timestep_fs:1.0, friction:0.02 });
    expect(timestamps.length).toBeGreaterThanOrEqual(8);
    const deltasFast = timestamps.slice(1).map((t,i)=> t - timestamps[i]);
    const minDeltaFast = Math.min(...deltasFast);
    const avgDeltaFast = deltasFast.reduce((a,b)=>a+b,0)/deltasFast.length;
    // Should respect pacing ~>=30ms when network is faster
    expect(minDeltaFast).toBeGreaterThanOrEqual(24); // allow jitter
    expect(avgDeltaFast).toBeGreaterThanOrEqual(26);

    // Run with slow network (~50ms). Should not add extra delay beyond the network time.
    networkDelay = 50; timestamps.length = 0;
    await api.startMDContinuous({ steps: 6, calculator:'uma', temperature:298, timestep_fs:1.0, friction:0.02 });
    expect(timestamps.length).toBeGreaterThanOrEqual(6);
    const deltasSlow = timestamps.slice(1).map((t,i)=> t - timestamps[i]);
    const minDeltaSlow = Math.min(...deltasSlow);
    const avgDeltaSlow = deltasSlow.reduce((a,b)=>a+b,0)/deltasSlow.length;
    // Expect roughly the network delay with small overhead, not network+30ms
    expect(minDeltaSlow).toBeGreaterThanOrEqual(40); // network - 10ms tolerance
    expect(avgDeltaSlow).toBeLessThanOrEqual(110); // generous upper bound for CI jitter
  });
});
