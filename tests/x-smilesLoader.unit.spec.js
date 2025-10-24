/** @jest-environment node */

import { smilesToXYZ, isLikelySmiles } from '../public/util/smilesLoader.js';

describe('smilesLoader basic', () => {
  const origFetch = global.fetch;

  afterEach(() => {
    global.fetch = origFetch;
  });

  test('smilesToXYZ fetches and converts simple molecule', async () => {
    const sdf = `
water
  PubChem

  3  2  0  0  0  0              0 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$
`.trim();
    global.fetch = jest.fn(async () => ({
      ok: true,
      status: 200,
      text: async () => sdf,
    }));
    const xyz = await smilesToXYZ('O');
    expect(xyz.split('\n')[0]).toBe('3');
    expect(xyz.includes('O 0 0 0')).toBe(true);
    expect(global.fetch).toHaveBeenCalledTimes(1);
  });

  test('isLikelySmiles filters simple invalid strings', () => {
    expect(isLikelySmiles('CCO')).toBe(true);
    expect(isLikelySmiles('C C O')).toBe(false);
    expect(isLikelySmiles('')).toBe(false);
  });
});
