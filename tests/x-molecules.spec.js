import fs from 'fs';
import path from 'path';
import { parseXYZ } from '../public/util/xyzLoader.js';

describe('x-molecules-fixtures', () => {
  function read(name) {
    return fs.readFileSync(path.join(process.cwd(), 'public', 'molecules', name), 'utf-8');
  }

  test('benzene fixture has expected composition', () => {
    const parsed = parseXYZ(read('benzene.xyz'));
    expect(parsed.elements.length).toBe(12);
    expect(parsed.elements.filter((e) => e === 'C').length).toBe(6);
  });

  test('roy fixture matches current atom count and includes sulfur', () => {
    const parsed = parseXYZ(read('roy.xyz'));
    expect(parsed.elements.length).toBe(27);
    expect(parsed.elements.includes('S')).toBe(true);
  });
});
