const { spawn, execSync } = require('child_process');
const fs = require('fs');
const net = require('net');

function waitForPort(port, { timeout = 20000, logPrefix='[waitForPort]' } = {}) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    let attempts = 0;
    (function tryOnce() {
      attempts++;
      const sock = net.createConnection(port, '127.0.0.1');
      sock.once('connect', () => { sock.destroy(); console.log(`${logPrefix} port ${port} up after ${attempts} attempts in ${Date.now()-start}ms`); resolve(true); });
      sock.once('error', () => {
        sock.destroy();
        if (Date.now() - start > timeout) return reject(new Error('Timeout waiting for port '+port));
        if (attempts % 15 === 0) console.log(`${logPrefix} still waiting on ${port} elapsed=${Date.now()-start}ms attempts=${attempts}`);
        setTimeout(tryOnce, 200);
      });
    })();
  });
}

function killPort(port){
  try {
    const out = execSync(`lsof -ti tcp:${port} || true`).toString().trim();
    if(out){
      out.split(/\n/).filter(Boolean).forEach(pid=>{
        try { process.kill(parseInt(pid), 'SIGKILL'); } catch {}
      });
    }
  } catch {}
}

async function cudaPreflight(pyEnv){
  return new Promise((resolve, reject)=>{
    const code = "import torch, json, os; print(json.dumps({'cuda_available': torch.cuda.is_available(), 'cuda_device_count': torch.cuda.device_count(), 'visible': os.environ.get('CUDA_VISIBLE_DEVICES')}), flush=True)";
    const p = spawn(pyEnv, ['-c', code], { env: { ...process.env } });
    let buf='';
    p.stdout.on('data', d=>buf+=d.toString());
    p.stderr.on('data', d=>buf+=d.toString());
    p.on('close', ()=>{
      try { const j = JSON.parse(buf.trim().split(/\n/).pop()); resolve(j); } catch(e){ reject(new Error('CUDA preflight parse failed: '+buf)); }
    });
  });
}

module.exports = async () => {
  const log = fs.createWriteStream('./test-server.log');
  // Expose log stream for teardown so it can be closed (open fs handle can keep Node alive).
  global.__MLIP_LOG_STREAM = log;
  const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd()+ '/mlipview_venv/bin/python';
  // Install global watchdog (moved from per-test setup to avoid open handle at test end)
  const MAX_WALL_MS = parseInt(process.env.MLIPVIEW_JEST_WALL_TIMEOUT || '90000', 10); // 90s default
  if(!global.__MLIPVIEW_WATCHDOG){
    const start = Date.now();
    global.__MLIPVIEW_WATCHDOG = setInterval(()=>{
      const elapsed = Date.now() - start;
      if(elapsed > MAX_WALL_MS){
        console.error(`[jest-watchdog] Exceeded wall timeout ${MAX_WALL_MS}ms; forcing exit.`);
        try { clearInterval(global.__MLIPVIEW_WATCHDOG); } catch{}
        try { const cleanups = global.__MLIPVIEW_CLEANUP; if(Array.isArray(cleanups)){ for(const fn of cleanups){ try{ fn(); } catch{} } } } catch {}
        process.exit(0);
      }
    }, 3000);
  }
  // Preflight CUDA
  const pre = await cudaPreflight(pyEnv).catch(e=>{ console.error('[globalSetup] CUDA preflight failed', e); return null; });
  if(!pre || !pre.cuda_available){
    console.error('[globalSetup] CUDA not available in preflight', pre);
    throw new Error('GPU-only mode but CUDA not available');
  }
  log.write('[preflight] '+JSON.stringify(pre)+'\n');

  // Kill any lingering processes on required ports
  killPort(8000);
  killPort(4000);

  // 1. Start Ray Serve backend (GPU-only enforced)
  const pyCode = "from fairchem_local_server.serve_app import deploy; deploy(); import time; print('SERVE_READY', flush=True);\nwhile True: time.sleep(60)";
  const pyArgs = ['-c', pyCode];
  const pyEnvVars = { ...process.env, UMA_DEVICE:'cuda' };
  // Do NOT inject empty CUDA_VISIBLE_DEVICES; only preserve if already set.
  if(process.env.CUDA_VISIBLE_DEVICES){
    pyEnvVars.CUDA_VISIBLE_DEVICES = process.env.CUDA_VISIBLE_DEVICES;
  }
  const pyProc = spawn(pyEnv, pyArgs, { stdio: ['ignore', 'pipe', 'pipe'], env: pyEnvVars });
  log.write('[env] starting Ray Serve with CUDA_VISIBLE_DEVICES='+(pyEnvVars.CUDA_VISIBLE_DEVICES||'UNSET')+'\n');
  pyProc.stdout.on('data', d => log.write('[py] '+d));
  pyProc.stderr.on('data', d => log.write('[py][err] '+d));
  global.__MLIP_PY_SERVER = pyProc;

  // 2. Start node static app via npx (node server.js)
  const nodeProc = spawn('npx', ['node', 'server.js'], { stdio: ['ignore', 'pipe', 'pipe'] });
  nodeProc.stdout.on('data', d => log.write('[node] '+d));
  nodeProc.stderr.on('data', d => log.write('[node][err] '+d));
  global.__MLIP_NODE_SERVER = nodeProc;

  // Wait for ports, then poll serve health
  try {
    console.log('[globalSetup] waiting for ports 8000 (serve) & 4000 (node)');
    await Promise.all([
      waitForPort(8000, { logPrefix:'[waitForPort:serve]' }),
      waitForPort(4000, { logPrefix:'[waitForPort:node]' })
    ]);
    const start = Date.now();
    let healthy = false;
    let healthPayload = null;
    let pollCount = 0;
    console.log('[globalSetup] begin health polling /serve/health (timeout 30000ms)');
    while(Date.now()-start < 30000){
      pollCount++;
      try {
        const r = await fetch('http://127.0.0.1:8000/serve/health');
        if(r.ok){ healthy=true; healthPayload = await r.json(); break; }
      } catch {}
      if(pollCount % 10 === 0) console.log(`[globalSetup] health still pending elapsed=${Date.now()-start}ms polls=${pollCount}`);
      await new Promise(r=>setTimeout(r,300));
    }
    if(!healthy){
      console.error('[globalSetup] Ray Serve health check failed after 30s');
      throw new Error('Ray Serve did not become healthy');
    }
    if(!healthPayload || healthPayload.device !== 'cuda' || healthPayload.cuda_available !== true){
      console.error('[globalSetup] Expected cuda device but got', healthPayload && healthPayload.device);
      throw new Error('Serve not running with CUDA device');
    }
    // Cache health snapshot for per-test assertions without repeated network fetches.
    global.__MLIP_HEALTH_SNAPSHOT = healthPayload;
  } catch (e) { console.error('Server startup wait failed', e); throw e; }

  // Provide base URLs for tests
  global.__MLIP_BASE_URL = 'http://localhost:4000';
  global.__MLIP_API_URL = 'http://localhost:8000';
};
