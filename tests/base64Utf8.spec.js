/** @jest-environment jsdom */

import { base64EncodeUtf8, base64DecodeUtf8 } from '../public/ui/moleculeSelect.js';
import { parseXYZ } from '../public/util/xyzLoader.js';

describe('UTF-8 base64 helpers and XYZ parsing', () => {
  test('UTF-8 comment survives encode/decode', () => {
    const comment = '水 H2O – тест 漢字';
    const xyz = `3\n${comment}\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n`;
    const b64 = base64EncodeUtf8(xyz);
    const roundTrip = base64DecodeUtf8(b64);
    expect(roundTrip).toBe(xyz);
    const parsed = parseXYZ(roundTrip);
    expect(parsed.comment).toBe(comment);
    expect(parsed.elements.length).toBe(3);
  });
});
