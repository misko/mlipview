import { execSync } from 'child_process';

async function isPortOpen(port) {
  return new Promise((resolve) => {
    const net = require('net');
    const s = net.createConnection(port, '127.0.0.1');
    s.once('connect', () => { s.destroy(); resolve(true); });
    s.once('error', () => { s.destroy(); resolve(false); });
  });
}

function killPort(port){
  try {
    const cmd = `bash -lc "(command -v lsof >/dev/null 2>&1 && lsof -t -i :${port}) || (command -v fuser >/dev/null 2>&1 && fuser ${port}/tcp 2>/dev/null) || true"`;
    const out = execSync(cmd, { stdio: ['ignore','pipe','pipe'] }).toString().trim();
    if (!out) return;
    const pids = out.split(/\s+/).map(s=>Number(s)).filter(Boolean);
    for (const pid of pids){ try { process.kill(pid, 'SIGTERM'); } catch{} }
  } catch {}
}

export default async function globalTeardown() {
  for (const key of ['__WS_PID__', '__VITE_PREVIEW_PID__']) {
    const pid = Number(process.env[key]);
    if (pid) {
      try { process.kill(pid, 'SIGTERM'); } catch {}
    }
  }
  // Best-effort: ensure ports are freed
  killPort(8000); killPort(5174);
  // small wait to let processes exit
  await new Promise(r=>setTimeout(r, 500));
}
