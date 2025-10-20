// Tests precomputed energy/forces injection for /serve/md and /serve/relax
// Ensures initial_energy echoes provided precomputed.energy and that
// precomputed_applied contains expected keys.

const API = 'http://127.0.0.1:8000/serve';

async function post(path, body) {
  const resp = await fetch(API + path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  if (!resp.ok) throw new Error(path + ' failed ' + resp.status);
  return resp.json();
}

const water = {
  atomic_numbers: [8, 1, 1],
  coordinates: [
    [0, 0, 0],
    [0.9575, 0, 0],
    [-0.2399872, 0.92662721, 0],
  ],
};

describe('precomputed results seeding', () => {
  test('md: precomputed.energy overrides initial_energy', async () => {
    const base = await post('/md', { ...water, steps: 1, calculator: 'lj' });
    const delta = 0.123456;
    const injected = base.initial_energy + delta;
    const seeded = await post('/md', {
      ...water,
      steps: 1,
      calculator: 'lj',
      precomputed: { energy: injected },
    });
    expect(seeded.initial_energy).toBeCloseTo(injected, 10);
    expect(seeded.precomputed_applied).toContain('energy');
  }, 15000);

  test('relax: precomputed.energy overrides initial_energy', async () => {
    const base = await post('/relax', { ...water, steps: 1, calculator: 'lj' });
    const delta = 0.111111;
    const injected = base.initial_energy + delta;
    const seeded = await post('/relax', {
      ...water,
      steps: 1,
      calculator: 'lj',
      precomputed: { energy: injected },
    });
    expect(seeded.initial_energy).toBeCloseTo(injected, 10);
    expect(seeded.precomputed_applied).toContain('energy');
  }, 15000);

  test('md: precomputed forces shape accepted', async () => {
    const zeroForces = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    const res = await post('/md', {
      ...water,
      steps: 1,
      calculator: 'lj',
      precomputed: { forces: zeroForces },
    });
    // We cannot guarantee they are used (no energy provided) but API should echo application
    expect(res.precomputed_applied).toContain('forces');
    expect(res.forces.length).toBe(3);
  }, 15000);
});
