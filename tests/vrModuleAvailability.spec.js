// Updated regression: ensure the VR entry module exists and defines initVRApp symbol (no legacy re-export).
import fs from 'fs';
import path from 'path';

describe('VR module availability', () => {
  test('public/vr/main-vr.js present and contains initVRApp export token', () => {
    const target = path.join(process.cwd(), 'public', 'vr', 'main-vr.js');
    expect(fs.existsSync(target)).toBe(true);
    const text = fs.readFileSync(target, 'utf8');
    expect(text.includes('export async function initVRApp')).toBe(true);
    // Guard: ensure no forbidden legacy re-export remains
    expect(text.includes('public_to_be_deleted')).toBe(false);
  });
});

