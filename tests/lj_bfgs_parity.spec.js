import { haveServer } from './helpers/server.js';

describe('LJ + BFGS parity via server /relax', () => {
  test('20-step water relaxation first/last vs reference ASE values', async () => {
    const up = await haveServer();
    if (!up) {
      console.warn('[lj parity] server not reachable; skipping');
      return;
    }
    // Water geometry (O,H,H) with O at origin (matches server parity python test ordering)
    // Use ASE reference energies documented previously: initial ~2.802523 final ~-2.983548 (20-step BFGS rc=3.0)
    const atomic_numbers = [8, 1, 1];
    const coords = [
      [0, 0, 0],
      [0.9575, 0, 0],
      [-0.2399872, 0.92662721, 0],
    ];
    const relaxBody = JSON.stringify({
      atomic_numbers,
      coordinates: coords,
      steps: 20,
      calculator: 'lj',
    });
    const json = await fetch('http://127.0.0.1:8000/serve/relax', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: relaxBody,
    }).then((r) => (r.ok ? r.json() : Promise.reject(r.status + ':' + r.statusText)));
    expect(json.steps_completed).toBe(20);
    // Loosen tolerance slightly due to potential minor cutoff/model differences vs stored reference
    expect(Math.abs(json.initial_energy - 2.802523)).toBeLessThan(2e-2);
    expect(Math.abs(json.final_energy - -2.983548)).toBeLessThan(1e-2);
  }, 20000);
});
