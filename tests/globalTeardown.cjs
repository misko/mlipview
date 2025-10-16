module.exports = async () => {
  // Run any registered cleanup callbacks (viewer loops, etc.)
  try {
    const cleanups = (global.window && global.window.__MLIPVIEW_CLEANUP) || global.__MLIPVIEW_CLEANUP;
    if (Array.isArray(cleanups)) {
      for (const fn of cleanups) { try { fn(); } catch {} }
    }
  } catch {}

  // Clear global watchdog if present
  try { if(global.__MLIPVIEW_WATCHDOG){ clearInterval(global.__MLIPVIEW_WATCHDOG); delete global.__MLIPVIEW_WATCHDOG; } } catch {}

  // Close log stream if open (open fs WriteStream prevents process exit)
  try { if(global.__MLIP_LOG_STREAM){
    const s = global.__MLIP_LOG_STREAM;
    if(!s.destroyed) {
      await new Promise(res=>{ try { s.end(()=>res()); setTimeout(res, 500); } catch { res(); } });
    }
  } } catch {}

  function killAndWait(proc, name){
    return new Promise(resolve => {
      if(!proc) return resolve();
      let settled = false;
      const timeout = setTimeout(()=>{
        if(!settled){
          try { proc.kill('SIGKILL'); } catch{}
          settled = true; resolve();
        }
      }, 3000);
      proc.on('close', ()=>{ if(!settled){ settled = true; clearTimeout(timeout); resolve(); } });
      try { proc.kill('SIGTERM'); } catch {}
    });
  }

  await Promise.all([
    killAndWait(global.__MLIP_NODE_SERVER, 'node-server'),
    killAndWait(global.__MLIP_PY_SERVER, 'py-server')
  ]);
};
