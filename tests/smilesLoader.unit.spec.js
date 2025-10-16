describe('smilesLoader basic', () => {
  test('isLikelySmiles accepts simple strings', async () => {
    const { isLikelySmiles } = await import('../public/util/smilesLoader.js');
    expect(isLikelySmiles('CCO')).toBe(true);
    expect(isLikelySmiles('C1=CC=CC=C1')).toBe(true);
    expect(isLikelySmiles('')).toBe(false);
    expect(isLikelySmiles('H2O')).toBe(true);
    expect(isLikelySmiles('with space')).toBe(false);
  });

  test('sdfToXYZ parses minimal SDF to XYZ', async () => {
    const { sdfToXYZ } = await import('../public/util/smilesLoader.js');
    const sdf = `\n  Mrv2014 07092006002D\n\n  3  2  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2400    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$`;
    const xyz = sdfToXYZ(sdf);
    expect(xyz.split(/\n/)[0]).toBe('3');
    expect(/O 0/.test(xyz)).toBe(true);
  });
});
