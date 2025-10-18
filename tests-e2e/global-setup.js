import { startServers } from '../server.js';
import { spawn } from 'child_process';
import fs from 'fs';
import net from 'net';
import http from 'http';

function waitForPort(port, { timeout = 20000 } = {}) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    (function tryOnce(){
      const s = net.createConnection(port, '127.0.0.1');
      s.once('connect', ()=>{ s.destroy(); resolve(true); });
      s.once('error', ()=>{ s.destroy(); if(Date.now()-start>timeout) return reject(new Error('Timeout waiting for '+port)); setTimeout(tryOnce,200); });
    })();
  });
}

async function isPortOpen(port){
  try { await waitForPort(port, { timeout: 200, }); return true; } catch { return false; }
}

async function waitForHealth(url, { timeout = 20000, interval = 200 } = {}) {
  const t0 = Date.now();
  while (Date.now() - t0 < timeout) {
    try {
      // Node 18+ has global fetch; fall back to http if missing
      if (typeof fetch === 'function') {
        const res = await fetch(url, { cache: 'no-store' });
        if (res.ok) {
          const j = await res.json().catch(()=>({}));
          if (j && (j.ok === true || j.status === 'ok')) return true;
        }
      } else {
        await new Promise((resolve, reject) => {
          const req = http.get(url, (res) => {
            if (res.statusCode && res.statusCode >= 200 && res.statusCode < 300) {
              let data = '';
              res.setEncoding('utf8');
              res.on('data', (c) => (data += c));
              res.on('end', () => {
                try {
                  const j = JSON.parse(data);
                  if (j && (j.ok === true || j.status === 'ok')) return resolve(true);
                } catch {}
                reject(new Error('bad body'));
              });
            } else {
              reject(new Error('status ' + res.statusCode));
            }
          });
          req.on('error', reject);
        });
        return true;
      }
    } catch {}
    await new Promise((r) => setTimeout(r, interval));
  }
  throw new Error('Timeout waiting for health at ' + url);
}

export default async function globalSetup() {
  process.env.NO_HTTPS = '1';
  // Node server: start only if 4000 not already bound
  const have4000 = await isPortOpen(4000).catch(()=>false);
  if (!have4000) {
    startServers();
    await waitForPort(4000);
  }
  // WS server (Ray Serve via our entrypoint) on 8000 for e2e tests
  const have8000 = await isPortOpen(8000).catch(()=>false);
  if (!have8000) {
    const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd() + '/mlipview_venv/bin/python';
    const log2 = fs.createWriteStream('./test-ws-e2e.log');
    const py2 = spawn(pyEnv, ['-m','fairchem_local_server2.serve_ws_app','--ngpus','0','--ncpus','2','--nhttp','1','--http-port','8000']);
    py2.stdout.on('data', d=> log2.write('[ws] '+d));
    py2.stderr.on('data', d=> log2.write('[ws][err] '+d));
    fs.writeFileSync('.e2e-ws.pid', String(py2.pid));
    await waitForPort(8000);
  }
  // Ensure health endpoint is OK before continuing
  await waitForHealth('http://127.0.0.1:8000/serve/health');
}
