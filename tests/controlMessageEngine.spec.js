import { createControlMessageEngine } from '../public/core/controlMessageEngine.js';

describe('controlMessageEngine', () => {
  const frameCount = 12;
  const offsetToIndex = (offset) => {
    const off = Math.floor(Number(offset));
    if (!Number.isFinite(off)) return null;
    let idx = frameCount + off;
    if (idx < 0) idx = 0;
    if (idx >= frameCount) idx = frameCount - 1;
    return idx;
  };
  const engine = createControlMessageEngine({
    resolveFrameIndexById: (id) => {
      if (typeof id !== 'string') return null;
      const numeric = Number(id.replace('frame-', ''));
      if (!Number.isFinite(numeric)) return null;
      return Math.max(0, Math.min(frameCount - 1, numeric - 1));
    },
    offsetToIndex,
    getFrameCount: () => frameCount,
  });

  engine.setMessages([
    {
      id: 'speed-focus',
      priority: 1,
      range: { start: { offset: -6 }, end: { offset: -2 } },
      actions: [
        { type: 'timeline.playbackSpeed', fps: 8 },
        {
          type: 'visual.opacityFocus',
          focus: { atoms: [0, 1], includeBonds: 'connected' },
          focusOpacity: { atoms: 1, bonds: 0.95 },
          backgroundOpacity: { atoms: 0.2, bonds: 0.1 },
        },
      ],
    },
    {
      id: 'annotate-window',
      priority: 0,
      range: { start: { offset: -4 }, end: { offset: -3 } },
      actions: [
        {
          type: 'overlay.callout',
          text: 'Window of interest\nObserve trajectory',
          textSize: 0.45,
          panelSize: { width: 1.2, height: 0.6 },
          anchor: { mode: 'atom', atom: 0, offset: [0, 0.8, 0] },
        },
      ],
    },
    {
      id: 'slow-down',
      priority: 2,
      range: { start: { offset: -3 }, end: { offset: -1 } },
      actions: [{ type: 'timeline.playbackSpeed', fps: 4 }],
    },
  ]);

  engine.updateFrameCount(frameCount);

  it('selects highest-priority actions per type for active frame', () => {
    const frameIndex = offsetToIndex(-3); // overlaps all three messages
    const result = engine.evaluate(frameIndex);
    expect(result.speed).toBeTruthy();
    expect(result.speed.fps).toBe(4); // slow-down overrides speed-focus due to higher priority
    expect(result.callout).toBeTruthy();
    expect(result.callout.text).toContain('Window of interest');
    expect(result.opacity).toBeTruthy();
    expect(result.opacity.focus.atoms).toEqual([0, 1]);
  });

  it('returns empty actions when no message applies', () => {
    const frameIndex = offsetToIndex(-11);
    const result = engine.evaluate(frameIndex);
    expect(result.speed).toBeNull();
    expect(result.callout).toBeNull();
    expect(result.opacity).toBeNull();
  });
});
