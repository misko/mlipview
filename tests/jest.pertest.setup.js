// Per-test setup: assert cached Ray Serve GPU health (populated in globalSetup) and log device/model load status.
beforeEach(async () => {
  let health = global.__MLIP_HEALTH_SNAPSHOT;
  if(!health){
    // Fallback: single on-demand fetch (should be rare). Use AbortController to avoid lingering sockets.
    try {
      const ac = new AbortController();
      const t = setTimeout(()=>ac.abort(), 2000);
      try {
        const r = await fetch((global.__MLIP_API_URL||'http://127.0.0.1:8000') + '/serve/health', { signal: ac.signal });
        if(r.ok){ health = await r.json(); global.__MLIP_HEALTH_SNAPSHOT = health; }
      } finally { clearTimeout(t); }
    } catch {}
  }
  if(!health) throw new Error('health snapshot missing');
  if(health.device !== 'cuda' || health.cuda_available !== true){
    throw new Error('GPU health check failed: '+JSON.stringify(health));
  }
  // eslint-disable-next-line no-console
  console.log(`[per-test-health] cuda_available=${health.cuda_available} model_loaded=${health.model_loaded}`);
});

afterEach(()=>{
  try {
    if(global.window && global.window.__MLIPVIEW_DEBUG_API){
      // eslint-disable-next-line no-console
      console.log('[per-test-health] disabling API debug logging');
      global.window.__MLIPVIEW_DEBUG_API = false;
    }
  } catch {}
});

afterAll(()=>{
  try {
    const w = global.window;
    if(w && w.__XR_HUD_POLL_INTERVAL){
      clearInterval(w.__XR_HUD_POLL_INTERVAL);
    }
    const cleanups = (w && w.__MLIPVIEW_CLEANUP) || global.__MLIPVIEW_CLEANUP;
    if(Array.isArray(cleanups)){
      for(const fn of cleanups){ try { fn(); } catch {} }
    }
  } catch {}
});
