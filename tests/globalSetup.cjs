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
    const code = (
      'import json, os, subprocess\n'
      + 'torch_ok=True\n'
      + 'torch_err=None\n'
      + 'try:\n'
      + '    import torch as _t\n'
      + 'except Exception as e:\n'
      + '    torch_ok=False\n'
      + '    torch_err=str(e)\n'
      + 'cuda_avail=False\n'
      + 'cuda_cnt=0\n'
      + 'torch_ver=None\n'
      + 'cuda_ver=None\n'
      + 'cudnn_ver=None\n'
      + 'dev_names=[]\n'
      + 'if torch_ok:\n'
      + '    try:\n'
      + '        cuda_avail = bool(_t.cuda.is_available())\n'
      + '        cuda_cnt = int(_t.cuda.device_count())\n'
      + '        torch_ver = getattr(_t, "__version__", None)\n'
      + '        cuda_ver = getattr(_t.version, "cuda", None)\n'
      + '        try:\n'
      + '            cudnn_ver = _t.backends.cudnn.version()\n'
      + '        except Exception:\n'
      + '            cudnn_ver = None\n'
      + '        try:\n'
      + '            dev_names = [_t.cuda.get_device_name(i) for i in range(cuda_cnt)]\n'
      + '        except Exception:\n'
      + '            dev_names = []\n'
      + '    except Exception as e:\n'
      + '        torch_err = str(e)\n'
      + 'try:\n'
      + '    nvsmi = subprocess.check_output(["bash","-lc","nvidia-smi -L || true"], text=True)\n'
      + 'except Exception:\n'
      + '    nvsmi = None\n'
      + 'out = {\n'
      + '  "torch_import_ok": torch_ok,\n'
      + '  "torch_error": torch_err,\n'
      + '  "torch_version": torch_ver,\n'
      + '  "torch_cuda_version": cuda_ver,\n'
      + '  "cuda_available": cuda_avail,\n'
      + '  "cuda_device_count": cuda_cnt,\n'
      + '  "cuda_device_names": dev_names,\n'
      + '  "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),\n'
      + '  "NV_GPU": os.environ.get("NV_GPU"),\n'
      + '  "nvidia_smi_L": nvsmi,\n'
      + '}\n'
      + 'print(json.dumps(out), flush=True)\n'
    );
    const p = spawn(pyEnv, ['-c', code], { env: { ...process.env } });
    let buf='';
    p.stdout.on('data', d=>buf+=d.toString());
    p.stderr.on('data', d=>buf+=d.toString());
    p.on('close', ()=>{
      const last = buf.trim().split(/\n/).pop();
      try { const j = JSON.parse(last); resolve(j); } catch(e){ reject(new Error('CUDA preflight parse failed: '+buf)); }
    });
  });
}

async function waitForServeHealth({ timeout = 5000, interval = 300 } = {}) {
  const start = Date.now();
  while (Date.now() - start < timeout) {
    try {
      const res = await fetch('http://127.0.0.1:8000/serve/health');
      if (res.ok) {
        const payload = await res.json();
        if (payload && payload.model_loaded === true) {
          return { ok: true, payload };
        }
      }
    } catch {}
    await new Promise((r) => setTimeout(r, interval));
  }
  return { ok: false };
}

module.exports = async () => {
  const sharedServerMode = process.env.MLIPVIEW_SHARED_SERVER === '1';
  const assignBaseUrls = () => {
    global.__MLIP_BASE_URL = 'http://localhost:4000';
    global.__MLIP_API_URL = 'http://localhost:8000';
  };
  assignBaseUrls();
  if (process.env.MLIPVIEW_FAST_JSDOM === '1') {
    return;
  }
  if (process.env.MLIPVIEW_SKIP_SERVERS === '1') {
    return; // skip heavy setup
  }
  if (sharedServerMode) {
    const health = await waitForServeHealth({ timeout: 5000 });
    if (health.ok) {
      console.log('[globalSetup] Reusing existing shared test server.');
      global.__MLIP_SHARED_KEEP = true;
      return;
    }
  }
  const log = fs.createWriteStream('./test-server.log');
  // Expose log stream for teardown so it can be closed (open fs handle can keep Node alive).
  global.__MLIP_LOG_STREAM = log;
  const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd()+ '/mlipview_venv/bin/python';
  // Install global watchdog (moved from per-test setup to avoid open handle at test end)
  if (!sharedServerMode) {
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
  }
  // Preflight CUDA
  const pre = await cudaPreflight(pyEnv).catch(e=>{ console.error('[globalSetup] CUDA preflight failed', e); return null; });
  console.log('[globalSetup] Python env:', pyEnv);
  console.log('[globalSetup] CUDA preflight:', pre);
  log.write('[preflight] '+JSON.stringify(pre)+'\n');
  if(!pre || !pre.cuda_available){
    console.error('[globalSetup] CUDA not available in preflight', pre);
    if(pre && pre.torch_import_ok && !pre.torch_cuda_version){
      console.error('[globalSetup] DETECTED CPU-only torch build. Install a CUDA-enabled torch (torch.version.cuda is null).');
    }
    if(pre && pre.nvidia_smi_L){
      console.error('[globalSetup] nvidia-smi -L output:\n'+pre.nvidia_smi_L);
    }
    throw new Error('GPU-only mode but CUDA not available');
  }

  // Kill any lingering processes on required ports
  killPort(8000);
  killPort(4000);

  // 1. Start Ray Serve backend (GPU-only enforced) in WS mode
  // Use fairchem_local_server2 WebSocket app with explicit ngpus=1 and http-port=8000
  const pyArgs = ['-m', 'fairchem_local_server2.serve_ws_app', '--ngpus', '1', '--http-port', '8000'];
  const pyEnvVars = { ...process.env, UMA_DEVICE:'cuda' };
  // Do NOT inject empty CUDA_VISIBLE_DEVICES; only preserve if already set.
  if(process.env.CUDA_VISIBLE_DEVICES){
    pyEnvVars.CUDA_VISIBLE_DEVICES = process.env.CUDA_VISIBLE_DEVICES;
  }
  const pyProc = spawn(pyEnv, pyArgs, { stdio: ['ignore', 'pipe', 'pipe'], env: pyEnvVars });
  log.write('[env] starting WS Ray Serve (fairchem_local_server2) with CUDA_VISIBLE_DEVICES='+(pyEnvVars.CUDA_VISIBLE_DEVICES||'UNSET')+'\n');
  console.log('[globalSetup] Spawning WS Ray Serve (fairchem_local_server2) ngpus=1 on :8000 with UMA_DEVICE=cuda and CUDA_VISIBLE_DEVICES='+(pyEnvVars.CUDA_VISIBLE_DEVICES||'UNSET'));
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
        if(r.ok){
          const hp = await r.json();
            if (hp && hp.device === 'cuda' && hp.model_loaded === true) {
              if (hp.cuda_available !== true) {
                if (pollCount % 5 === 0) console.warn('[globalSetup] WARNING: cuda_available reported false in health but proceeding (Option A)');
              }
              healthy = true; healthPayload = hp; break;
            }
            if (pollCount % 5 === 0) console.log('[globalSetup] waiting for model/device readiness', hp);
        }
      } catch (e) {
        if (pollCount % 10 === 0) console.log('[globalSetup] health fetch error', e?.message||String(e));
      }
      if(pollCount % 10 === 0) console.log(`[globalSetup] health still pending elapsed=${Date.now()-start}ms polls=${pollCount}`);
      await new Promise(r=>setTimeout(r,300));
    }
    if(!healthy){
      console.error('[globalSetup] Ray Serve health check failed after 30s last payload=', healthPayload);
      throw new Error('Ray Serve did not become healthy with CUDA');
    }
    // Cache health snapshot for per-test assertions without repeated network fetches.
    global.__MLIP_HEALTH_SNAPSHOT = healthPayload;
  } catch (e) { console.error('Server startup wait failed', e); throw e; }

  // Provide base URLs for tests
  global.__MLIP_BASE_URL = 'http://localhost:4000';
  global.__MLIP_API_URL = 'http://localhost:8000';
  if (sharedServerMode) {
    global.__MLIP_SHARED_KEEP = true;
  }
};
