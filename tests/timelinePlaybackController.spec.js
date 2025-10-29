import { createTimelinePlaybackController } from '../public/core/timelinePlaybackController.js';

describe('timelinePlaybackController', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  });

  afterEach(() => {
    jest.clearAllTimers();
    jest.useRealTimers();
  });

  it('advances according to fps and applies overrides', () => {
    const invocations = [];
    const controller = createTimelinePlaybackController({
      defaultFps: 10,
      onStep: ({ fps }) => {
        invocations.push(fps);
        return invocations.length < 3;
      },
      schedule: (fn, ms) => setTimeout(fn, ms),
      cancel: (id) => clearTimeout(id),
    });

    controller.start();
    jest.advanceTimersByTime(100); // 1000 / 10 fps
    expect(invocations).toHaveLength(1);
    expect(invocations[0]).toBe(10);

    controller.setSpeedOverride({ fps: 5, sourceId: 'test' });
    jest.advanceTimersByTime(200); // 1000 / 5 fps
    expect(invocations).toHaveLength(2);
    expect(invocations[1]).toBe(5);

    controller.stop();
    jest.advanceTimersByTime(500);
    expect(invocations).toHaveLength(2);
  });

  it('captures and applies playback snapshot data', () => {
    const controller = createTimelinePlaybackController();
    controller.setBaseConfig({
      defaultFps: 22,
      autoPlay: true,
      loop: true,
      loopRange: { start: { offset: -10 }, end: { offset: -2 } },
      startFrame: { frameId: 'frame-0005' },
    });
    const snapshot = controller.getSnapshot();
    expect(snapshot.defaultFps).toBe(22);
    expect(snapshot.autoPlay).toBe(true);
    expect(snapshot.loop).toBe(true);
    expect(snapshot.loopRange.start.offset).toBe(-10);

    controller.applySnapshot({ defaultFps: 12, autoPlay: false });
    const updated = controller.getSnapshot();
    expect(updated.defaultFps).toBe(12);
    expect(updated.autoPlay).toBe(false);
    expect(updated.loop).toBe(true); // unchanged
  });
});
