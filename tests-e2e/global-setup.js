import { startServers } from '../server.js';
import { spawn } from 'child_process';
import fs from 'fs';
import net from 'net';

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

export default async function globalSetup() {
  process.env.NO_HTTPS = '1';
  // Node server: start only if 4000 not already bound
  const have4000 = await isPortOpen(4000).catch(()=>false);
  if (!have4000) {
    startServers();
    await waitForPort(4000);
  }
  // Python backend: start only if 8000 not already bound
  const have8000 = await isPortOpen(8000).catch(()=>false);
  if (!have8000) {
    const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd() + '/mlipview_venv/bin/python';
    const log = fs.createWriteStream('./test-server-e2e.log');
    const py = spawn(pyEnv, ['-m','uvicorn','fairchem_local_server.server:app','--host','0.0.0.0','--port','8000']);
    py.stdout.on('data', d=> log.write('[py] '+d));
    py.stderr.on('data', d=> log.write('[py][err] '+d));
    fs.writeFileSync('.e2e-py.pid', String(py.pid));
    await waitForPort(8000);
  }
}
