const { createFairChemForceField } = require('../public/forcefield/fairchem.js');

// Mock fetch for unit test (phase 1)
global.fetch = async (url, opts) => {
  const body = JSON.parse(opts.body);
  return {
    ok: true,
    status: 200,
    statusText: 'OK',
    async json() {
      return {
        results: {
          energy: -10.123,
          forces: body.atomic_numbers.map(() => [0.1, -0.2, 0.3])
        }
      };
    },
    async text() { return ''; }
  };
};

// CommonJS test file
describe('fairchem adapter (unit mock)', () => {
  it('returns energy and forces with proper lengths', async () => {
    const ff = createFairChemForceField({ endpoint: 'http://mock/simple_calculate' });
    const Z = [8,1,1];
    const xyz = [
      [0,0,0.119262],
      [0,0.763239,-0.477047],
      [0,-0.763239,-0.477047]
    ];
    const { energy, forces } = await ff.computeRaw({ Z, xyz });
    expect(typeof energy).toBe('number');
    expect(forces).toHaveLength(Z.length);
    forces.forEach(f => expect(f).toHaveLength(3));
  });
});
