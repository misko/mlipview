const { spawn } = require('child_process');
const fs = require('fs');
const net = require('net');

function waitForPort(port, { timeout = 15000 } = {}) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    (function tryOnce() {
      const sock = net.createConnection(port, '127.0.0.1');
      sock.once('connect', () => { sock.destroy(); resolve(true); });
      sock.once('error', () => {
        sock.destroy();
        if (Date.now() - start > timeout) return reject(new Error('Timeout waiting for port '+port));
        setTimeout(tryOnce, 200);
      });
    })();
  });
}

module.exports = async () => {
  const log = fs.createWriteStream('./test-server.log');
  // 1. Start backend Python API (uvicorn)
  const pyEnv = process.env.MLIPVIEW_PYTHON || process.cwd()+ '/mlipview_venv/bin/python';
  const pyArgs = ['-m', 'uvicorn', 'fairchem_local_server.server:app', '--host', '0.0.0.0', '--port', '8000'];
  const pyProc = spawn(pyEnv, pyArgs, { stdio: ['ignore', 'pipe', 'pipe'] });
  pyProc.stdout.on('data', d => log.write('[py] '+d));
  pyProc.stderr.on('data', d => log.write('[py][err] '+d));
  global.__MLIP_PY_SERVER = pyProc;

  // 2. Start node static app via npx (node server.js)
  const nodeProc = spawn('npx', ['node', 'server.js'], { stdio: ['ignore', 'pipe', 'pipe'] });
  nodeProc.stdout.on('data', d => log.write('[node] '+d));
  nodeProc.stderr.on('data', d => log.write('[node][err] '+d));
  global.__MLIP_NODE_SERVER = nodeProc;

  // Wait for both ports: 8000 (python) and 4000 (node)
  try {
    await Promise.all([waitForPort(8000), waitForPort(4000)]);
  } catch (e) {
    console.error('Server startup wait failed', e);
  }

  // Provide base URLs for tests
  global.__MLIP_BASE_URL = 'http://localhost:4000';
  global.__MLIP_API_URL = 'http://localhost:8000';
};
