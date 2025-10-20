import { spawn } from 'child_process';
import fs from 'fs';
import net from 'net';

function waitForPort(port, timeout = 20000) {
  const t0 = Date.now();
  return new Promise((resolve, reject) => {
    (function probe() {
      const s = net.createConnection(port, '127.0.0.1');
      s.once('connect', () => { s.destroy(); resolve(true); });
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
    s.once('connect', () => { s.destroy(); resolve(true); });
    s.once('error', () => { s.destroy(); resolve(false); });
  });
}

export default async function globalSetup() {
  process.env.NO_VITE_HTTPS = '1';

  // 1) Start backend (WS/FastAPI) on 8000
  const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd() + '/mlipview_venv/bin/python';
  const wsLog = fs.createWriteStream('./test-ws-e2e.log');
  const port8000Busy = await isPortOpen(8000);
  if (!port8000Busy) {
    const env = { ...process.env };
    // Pass through UMA geom debug toggle to backend if set
    if (process.env.UMA_GEOM_DEBUG == null && process.env.BACKEND_DEBUG_GEOM) {
      env.UMA_GEOM_DEBUG = process.env.BACKEND_DEBUG_GEOM;
    }
    const py = spawn(pyEnv,
      ['-m','fairchem_local_server2.serve_ws_app','--ngpus','1','--ncpus','2','--nhttp','1','--http-port','8000'],
      { env }
    );
    py.stdout.pipe(wsLog);
    py.stderr.pipe(wsLog);
    process.env.__WS_PID__ = String(py.pid);
    await waitForPort(8000);
  }

  // 2) Build frontend
  await new Promise((resolve, reject) => {
    const b = spawn(process.platform === 'win32' ? 'npm.cmd' : 'npm', ['run','build'], { stdio: 'inherit' });
    b.on('exit', code => code === 0 ? resolve() : reject(new Error('vite build failed')));
  });

  // 3) Serve dist/ with vite preview on 5174
  const port5174Busy = await isPortOpen(5174);
  if (!port5174Busy) {
    const preview = spawn(process.platform === 'win32' ? 'npx.cmd' : 'npx',
      ['vite','preview','--strictPort','--port','5174'],
      { stdio: 'inherit' }
    );
    process.env.__VITE_PREVIEW_PID__ = String(preview.pid);
    await waitForPort(5174);
  }

  // 4) Let tests know where to go
  process.env.BASE_URL = 'http://127.0.0.1:5174';
}
