import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
import { parseXYZ } from '../public/util/xyzLoader.js';

describe('molecule XYZ files', () => {
  function read(name){ return fs.readFileSync(path.join(__dirname,'../public/molecules', name),'utf-8'); }
  test('benzene has 12 atoms', () => {
    const p = parseXYZ(read('benzene.xyz'));
    expect(p.elements.length).toBe(12);
    expect(p.elements.filter(e=>e==='C').length).toBe(6);
  });
  test('roy has 28 atoms', () => {
    const p = parseXYZ(read('roy.xyz'));
    expect(p.elements.length).toBe(28);
    expect(p.elements.includes('S')).toBe(true);
  });
});
