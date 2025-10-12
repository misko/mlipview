/** @jest-environment jsdom */

// Test: ?debugLatency=1 logs timing between requests and per-iteration breakdown

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

function wait(ms){ return new Promise(r=>setTimeout(r, ms)); }

describe('debugLatency logs inter-request timing', () => {
  test('logs sincePrev and iter breakdown', async () => {
    // Enable latency debug via search param and flag on existing jsdom window
    try { window.history.pushState({}, '', '?debugLatency=1'); } catch {}
    window.__MLIP_DEBUG_LATENCY = true;
    window.__MLIPVIEW_TEST_MODE = true;
    window.__MLIPVIEW_NO_AUTO_MD = true;
    // Fast minStepInterval so we get multiple iterations
    window.__MLIP_CONFIG = { minStepIntervalMs: 10, mdFriction: 0.5 };
    document.body.innerHTML = '<canvas id="viewer"></canvas><div id="app"></div>';

    // Capture console lines
    const logs = []; const origLog = console.log; const origWarn = console.warn;
    console.log = (...a)=>{ logs.push(['log', a.map(String).join(' ')]); origLog.apply(console, a); };
    console.warn = (...a)=>{ logs.push(['warn', a.map(String).join(' ') ]); origWarn.apply(console, a); };

    // Minimal fetch stub to simulate endpoints quickly with a small delay to produce non-zero times
    global.fetch = async function(url, opts){
      const isRelax = /\/serve\/relax$|(^|\/)relax$/.test(String(url));
      const isMd = /\/serve\/md$|(^|\/)md$/.test(String(url));
      await wait(5);
      if (isRelax) return { ok:true, status:200, json: async ()=> ({ positions:[[0,0,0]], forces:[[0,0,0]], final_energy: -1.0 }) };
      if (isMd) return { ok:true, status:200, json: async ()=> ({ positions:[[0,0,0]], forces:[[0,0,0]], final_energy: -1.0, velocities:[[0,0,0]], temperature: 300 }) };
      return { ok:true, status:200, json: async ()=> ({ results:{} }) };
    };

    const { initNewViewer, setMinStepInterval } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const viewer = await initNewViewer(canvas, { elements:[{Z:8}], positions:[{x:0,y:0,z:0}], bonds:[] });
    // Run a very short MD sequence to trigger multiple requests
    await viewer.startMDContinuous({ steps: 3, temperature: 300 });

    // Give log pipeline time to flush
    await wait(20);

    // Find latency lines
    const iterBreakdown = logs.filter(l => /\[latency:md\] iter/.test(l[1]));
  // We expect at least 1 iteration breakdown entry
    expect(iterBreakdown.length).toBeGreaterThanOrEqual(1);
    expect(iterBreakdown.some(l => /gapSincePrevSend/.test(l[1]))).toBe(true);

    // restore
    console.log = origLog; console.warn = origWarn;
  });
});
