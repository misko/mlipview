import { createFrameBuffer } from '../public/core/frameBuffer.js';

describe('frameBuffer export/import', () => {
  it('round-trips frames with metadata intact', () => {
    const buffer = createFrameBuffer({ capacity: 3 });
    const frameA = {
      positions: [[0, 0, 0], [1, 0, 0]],
      velocities: [[0, 0, 0], [0, 0, 0]],
      forces: [[0.1, 0, 0], [0, 0.1, 0]],
      energy: -10,
      userInteractionCount: 2,
      simStep: 1,
      seq: 11,
    };
    const frameB = {
      positions: [[0, 1, 0], [1, 1, 0]],
      velocities: [[0, 0, 0], [0, 0, 0]],
      forces: [[0, 0.1, 0], [0, 0.2, 0]],
      energy: -9.5,
      userInteractionCount: 3,
      simStep: 2,
      seq: 12,
    };
    buffer.record('md', frameA, { energyIndex: 0, timestamp: 100 });
    buffer.record('md', frameB, { energyIndex: 1, timestamp: 200 });

    const exported = buffer.exportFrames();
    expect(exported).toHaveLength(2);
    expect(exported[0].seq).toBe(11);
    expect(exported[1].energyIndex).toBe(1);

    const clone = createFrameBuffer({ capacity: 3 });
    clone.importFrames(exported);

    expect(clone.size()).toBe(2);
    const latest = clone.getByOffset(-1);
    expect(latest.seq).toBe(12);
    expect(latest.userInteractionCount).toBe(3);
    expect(latest.positions[0][1]).toBeCloseTo(1);
    expect(latest.energy).toBeCloseTo(-9.5);
  });
});
