import { execSync } from 'child_process';

const SKIP_BACKEND = process.env.MLIPVIEW_SKIP_SERVERS === '1';

async function isPortOpen(port) {
  return new Promise((resolve) => {
    const net = require('net');
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

function killPort(port) {
  try {
    const cmd = `bash -lc "(command -v lsof >/dev/null 2>&1 && lsof -t -i :${port}) || (command -v fuser >/dev/null 2>&1 && fuser ${port}/tcp 2>/dev/null) || true"`;
    const out = execSync(cmd, { stdio: ['ignore', 'pipe', 'pipe'] })
      .toString()
      .trim();
    if (!out) return;
    const pids = out
      .split(/\s+/)
      .map((s) => Number(s))
      .filter(Boolean);
    for (const pid of pids) {
      try {
        process.kill(pid, 'SIGTERM');
      } catch {}
    }
  } catch {}
}

export default async function globalTeardown() {
  if (!SKIP_BACKEND) {
    const wsPid = Number(process.env.__WS_PID__);
    if (wsPid) {
      try {
        process.kill(wsPid, 'SIGTERM');
      } catch {}
    }
    killPort(8000);
  } else {
    console.log('[playwright][teardown] MLIPVIEW_SKIP_SERVERS=1 — leaving backend on :8000 running');
  }

  const previewPid = Number(process.env.__VITE_PREVIEW_PID__);
  if (previewPid) {
    try {
      process.kill(previewPid, 'SIGTERM');
    } catch {}
  }
  // Best-effort: ensure preview port is freed (frontend always managed by tests)
  killPort(5174);

  // Best-effort cleanup for any other tracked processes
  for (const key of ['__VITE_PREVIEW_PID__']) {
    const pid = Number(process.env[key]);
    if (pid) {
      try {
        process.kill(pid, 'SIGTERM');
      } catch {}
    }
  }

  // small wait to let processes exit
  await new Promise((r) => setTimeout(r, 500));
}
