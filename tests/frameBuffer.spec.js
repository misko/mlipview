import { createFrameBuffer } from '../public/core/frameBuffer.js';

describe('frameBuffer', () => {
  it('stores frames up to capacity and exposes offsets/signatures', () => {
    const buf = createFrameBuffer({ capacity: 3 });
    const makeFrame = (offset) => ({
      seq: Math.abs(offset),
      positions: Array.from({ length: 2 }, (_, i) => [offset + i, offset + i + 0.5, offset - i]),
      energy: offset * 1.23,
      kind: 'md',
    });

    buf.record('md', makeFrame(-1));
    buf.record('md', makeFrame(-2));
    let stats = buf.stats();
    expect(stats.size).toBe(2);
    expect(buf.listOffsets()).toEqual([-1, -2]);

    const latest = buf.getByOffset(-1);
    expect(latest.energy).toBeCloseTo(-2 * 1.23);
    expect(Array.isArray(latest.positions)).toBe(true);
    expect(latest.positions.length).toBe(2);

    const sig = buf.getSignature(-1);
    buf.record('md', makeFrame(-3));
    buf.record('md', makeFrame(-4)); // overwrite oldest
    stats = buf.stats();
    expect(stats.size).toBe(3);
    expect(buf.listOffsets()).toEqual([-1, -2, -3]);
    expect(buf.getSignature(-1)).not.toBeNull();
    expect(buf.getSignature(-4)).toBeNull();
    expect(sig).not.toBeNull();
  });
});
