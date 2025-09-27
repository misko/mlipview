const { spawn } = require('child_process');
const http = require('http');
const { createFairChemForceField } = require('../public/forcefield/fairchem.js');

const PORT = 8110;
const ENDPOINT = `http://localhost:${PORT}/simple_calculate`;

function waitForServer(url, timeoutMs = 120000) {
  const start = Date.now();
  return new Promise((resolve, reject) => {
    const attempt = () => {
      const req = http.request(url, { method: 'POST' }, res => {
        // We expect 422/400 if body missing – that still means server is up
        if (res.statusCode && res.statusCode < 600) {
          res.resume();
          resolve();
        } else {
          res.resume();
          retry();
        }
      });
      req.on('error', retry);
      req.end();
    };
    const retry = () => {
      if (Date.now() - start > timeoutMs) return reject(new Error('Server did not become ready'));
      setTimeout(attempt, 750);
    };
    attempt();
  });
}

describe('fairchem adapter integration (local server)', () => {
  let proc;
  beforeAll(async () => {
    // spawn uvicorn using venv python if available else system python
    const py = process.env.FAIRCHEM_PYTHON || 'python3';
    proc = spawn(py, ['-m', 'uvicorn', 'fairchem_local_server.server:app', '--port', String(PORT)], {
      stdio: 'inherit'
    });
    await waitForServer(ENDPOINT);
  }, 130000);

  afterAll(() => {
    if (proc && !proc.killed) proc.kill();
  });

  it('computes water energy & forces', async () => {
    const ff = createFairChemForceField({ endpoint: ENDPOINT });
    const Z = [8,1,1];
    const xyz = [
      [0,0,0.119262],
      [0,0.763239,-0.477047],
      [0,-0.763239,-0.477047]
    ];
    const { energy, forces } = await ff.computeRaw({ Z, xyz });
    expect(Number.isFinite(energy)).toBe(true);
    expect(forces).toHaveLength(Z.length);
    for (const f of forces) {
      expect(f).toHaveLength(3);
      f.forEach(v => expect(Number.isFinite(v)).toBe(true));
    }
  }, 60000);
});
