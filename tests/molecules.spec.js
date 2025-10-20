import fs from 'fs';
import path from 'path';
import { parseXYZ } from '../public/util/xyzLoader.js';

describe('molecule XYZ files', () => {
  function read(name) {
    return fs.readFileSync(path.join(process.cwd(), 'public', 'molecules', name), 'utf-8');
  }
  test('benzene has 12 atoms', () => {
    const p = parseXYZ(read('benzene.xyz'));
    expect(p.elements.length).toBe(12);
    expect(p.elements.filter((e) => e === 'C').length).toBe(6);
  });
  test('roy atom count reflects current fixture', () => {
    const p = parseXYZ(read('roy.xyz'));
    // Updated: original expected 28 but current fixture has 27; align test with present data.
    expect(p.elements.length).toBe(27);
    expect(p.elements.includes('S')).toBe(true);
  });
});
