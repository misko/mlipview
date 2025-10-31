import { spawn, execSync } from 'child_process';
import fs from 'fs';
import net from 'net';

function waitForPort(port, timeout = 20000) {
  const t0 = Date.now();
  return new Promise((resolve, reject) => {
    (function probe() {
      const s = net.createConnection(port, '127.0.0.1');
      s.once('connect', () => {
        s.destroy();
        resolve(true);
      });
      s.once('error', () => {
        s.destroy();
        if (Date.now() - t0 > timeout) reject(new Error(`Timeout waiting for ${port}`));
        else setTimeout(probe, 200);
      });
    })();
  });
}

async function isPortOpen(port) {
  return new Promise((resolve) => {
    const s = net.createConnection(port, '127.0.0.1');
    s.once('connect', () => {
      s.destroy();
      resolve(true);
    });
    s.once('error', () => {
      s.destroy();
      resolve(false);
    });
  });
}

function pidsOnPort(port) {
  try {
    // Linux: lsof fallback to fuser; ignore failures
    const cmd = `bash -lc "(command -v lsof >/dev/null 2>&1 && lsof -t -i :${port}) || (command -v fuser >/dev/null 2>&1 && fuser ${port}/tcp 2>/dev/null) || true"`;
    const out = execSync(cmd, { stdio: ['ignore', 'pipe', 'pipe'] })
      .toString()
      .trim();
    if (!out) return [];
    return out
      .split(/\s+/)
      .map((s) => Number(s))
      .filter(Boolean);
  } catch {
    return [];
  }
}

async function freePort(port, timeout = 10000) {
  const busy = await isPortOpen(port);
  if (!busy) return true;
  const pids = pidsOnPort(port);
  for (const pid of pids) {
    try {
      process.kill(pid, 'SIGTERM');
    } catch { }
  }
  const t0 = Date.now();
  while (Date.now() - t0 < timeout) {
    const still = await isPortOpen(port);
    if (!still) return true;
    await new Promise((r) => setTimeout(r, 200));
  }
  // escalate
  for (const pid of pids) {
    try {
      process.kill(pid, 'SIGKILL');
    } catch { }
  }
  const t1 = Date.now();
  while (Date.now() - t1 < timeout) {
    const still = await isPortOpen(port);
    if (!still) return true;
    await new Promise((r) => setTimeout(r, 200));
  }
  throw new Error(`Could not free port ${port}`);
}

const SKIP_BACKEND = process.env.MLIPVIEW_SKIP_SERVERS === '1';

export default async function globalSetup() {
  process.env.NO_VITE_HTTPS = '1';

  // 1) Start or verify backend (WS/FastAPI) on 8000
  if (SKIP_BACKEND) {
    console.log('[playwright][setup] MLIPVIEW_SKIP_SERVERS=1 — reusing existing backend on :8000');
    // Ensure the external backend is ready before continuing.
    await waitForPort(8000).catch((err) => {
      throw new Error(`Expected backend on port 8000 when MLIPVIEW_SKIP_SERVERS=1, but it was not reachable: ${err?.message || err}`);
    });
  } else {
    const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd() + '/mlipview_venv/bin/python';
    const wsLog = fs.createWriteStream('./test-ws-e2e.log');
    // Ensure port 8000 is free before launching
    await freePort(8000).catch(() => { });
    const port8000Busy = await isPortOpen(8000);
    if (!port8000Busy) {
      const env = { ...process.env };
      // Enable backend WS debug to surface WAITING_FOR_ACK and frame logs during tests
      env.WS_DEBUG = env.WS_DEBUG || '1';
      // Pass through UMA geom debug toggle to backend if set
      if (process.env.UMA_GEOM_DEBUG == null && process.env.BACKEND_DEBUG_GEOM) {
        env.UMA_GEOM_DEBUG = process.env.BACKEND_DEBUG_GEOM;
      }
      const py = spawn(
        pyEnv,
        [
          '-m',
          'fairchem_local_server2.serve_ws_app',
          '--ngpus',
          '1',
          '--ncpus',
          '2',
          '--nhttp',
          '1',
          '--http-port',
          '8000',
        ],
        { env }
      );
      py.stdout.pipe(wsLog);
      py.stderr.pipe(wsLog);
      process.env.__WS_PID__ = String(py.pid);
      await waitForPort(8000);
    } else {
      console.log('[playwright][setup] Port 8000 already in use; assuming compatible backend is running');
      await waitForPort(8000);
    }
  }

  // 2) Build frontend
  await new Promise((resolve, reject) => {
    const b = spawn(process.platform === 'win32' ? 'npm.cmd' : 'npm', ['run', 'build'], {
      stdio: 'inherit',
    });
    b.on('exit', (code) => (code === 0 ? resolve() : reject(new Error('vite build failed'))));
  });

  // 3) Serve dist/ with vite preview on 5174
  await freePort(5174).catch(() => { });
  const port5174Busy = await isPortOpen(5174);
  if (!port5174Busy) {
    const preview = spawn(
      process.platform === 'win32' ? 'npx.cmd' : 'npx',
      ['vite', 'preview', '--strictPort', '--port', '5174'],
      { stdio: 'inherit' }
    );
    process.env.__VITE_PREVIEW_PID__ = String(preview.pid);
    await waitForPort(5174);
  }

  // 4) Let tests know where to go
  process.env.BASE_URL = 'http://127.0.0.1:5174';
}
