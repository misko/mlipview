import fs from 'fs';
import path from 'path';

describe('x-vr module availability', () => {
  test('main-vr.js exists with initVRApp export', () => {
    const target = path.join(process.cwd(), 'public', 'vr', 'main-vr.js');
    expect(fs.existsSync(target)).toBe(true);
    const source = fs.readFileSync(target, 'utf8');
    expect(source.includes('export async function initVRApp')).toBe(true);
  });
});
