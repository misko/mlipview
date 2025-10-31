import { test, expect } from './fixtures.js';
import {
  hoverTimelinePanel,
  waitForTimelineBuffer,
  waitForTimelineState,
  selectTimelineOffset,
  computeTimelineOffset,
} from './utils/timeline.js';

test('timeline control policy enforcement', async ({ page, loadViewerPage, uiMode }) => {
  test.setTimeout(180_000);
  console.log(`[timeline-controls] starting test (uiMode=${uiMode})`);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
  console.log('[timeline-controls] viewer page loaded');

  await page.evaluate(() => {
    console.log('[timeline-controls][browser] starting MD continuous run');
    window.viewerApi.startMDContinuous({ steps: 400, temperature: 1500 });
  });

  await waitForTimelineBuffer(page, 60, { timeout: 45_000, label: 'buffer>=60' });

  await hoverTimelinePanel(page, 'initialHover');

  await page.evaluate(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    console.log('[timeline-controls][browser] status after hover', status);
  });

  const playBtn = page.locator('[data-testid="timeline-play"]');
  const pauseBtn = page.locator('[data-testid="timeline-pause"]');
  const liveBtn = page.locator('[data-testid="timeline-live"]');

  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeDisabled();
  await expect(liveBtn).toBeDisabled();

  const targetOffset =
    (await computeTimelineOffset(page, { offset: -12 })) ??
    (await computeTimelineOffset(page, { mode: 'middle' }));
  expect(targetOffset).not.toBeNull();
  console.log(`[timeline-controls] computed target offset ${targetOffset}`);

  await selectTimelineOffset(page, targetOffset, { label: 'initial-select' });

  await waitForTimelineState(
    page,
    (expected) => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.offset === expected && status.mode !== 'live' && !status.playing;
    },
    targetOffset,
    { label: 'waitPausedInitial' }
  );

  await page.evaluate(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    const snapshot = {
      status,
      loop: timeline?.getLoopConfig?.(),
      fps:
        timeline?.getPlaybackConfig?.()?.defaultFps ??
        timeline?.getPlaybackSnapshot?.()?.defaultFps ??
        null,
    };
    console.log('[timeline-controls][browser] after select', snapshot);
  });

  await expect(playBtn).toBeEnabled();
  await expect(pauseBtn).toBeDisabled();
  await expect(liveBtn).toBeEnabled();

  await playBtn.click();
  console.log('[timeline-controls] play clicked');
  await page.waitForTimeout(250);
  await page.evaluate(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    console.log('[timeline-controls][browser] status after play', status);
  });

  await waitForTimelineState(
    page,
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.mode !== 'live' && !!status.playing;
    },
    null,
    { label: 'waitPlayingPostPlay' }
  );

  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeEnabled();

  await pauseBtn.click();
  console.log('[timeline-controls] pause clicked');

  await waitForTimelineState(
    page,
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.mode !== 'live' && !status.playing;
    },
    null,
    { label: 'waitPausedPostPause' }
  );

  await playBtn.click();
  console.log('[timeline-controls] play clicked again');
  await page.waitForTimeout(250);
  await page.evaluate(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    console.log('[timeline-controls][browser] status after resume', status);
  });

  await waitForTimelineState(
    page,
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.mode !== 'live' && !!status.playing;
    },
    null,
    { label: 'waitPlayingPostResume' }
  );

  await waitForTimelineState(
    page,
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.mode === 'live' && !status.active;
    },
    null,
    { timeout: 30_000, label: 'waitReturnToLive' }
  );

  await hoverTimelinePanel(page, 'postLiveHover');
  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeDisabled();

  if (uiMode === 'react') {
    const advanced = page.locator('[data-testid="timeline-advanced"]');
    await expect(advanced).toBeVisible();

    const loopToggle = advanced.locator('[data-testid="timeline-loop-toggle"]');
    await loopToggle.check();
    await page.waitForFunction(() => {
      const cfg = window.viewerApi.timeline.getLoopConfig?.();
      return !!cfg && cfg.loop === true;
    });

    await hoverTimelinePanel(page);
    await selectTimelineOffset(page, targetOffset);

    const startBtn = advanced.locator('[data-testid="timeline-loop-start"]');
    await startBtn.click();

    const endCandidate =
      (await computeTimelineOffset(page, { offset: Number(targetOffset) - 6 })) ??
      (await computeTimelineOffset(page, { mode: 'middle' }));
    await selectTimelineOffset(page, endCandidate);

    const endBtn = advanced.locator('[data-testid="timeline-loop-end"]');
    await endBtn.click();
    console.log('[timeline-controls] loop end set');

    await advanced.locator('[data-testid="timeline-loop-apply"]').click();
    console.log('[timeline-controls] loop apply clicked');

    const expectedRange = {
      start: Math.min(Number(targetOffset), Number(endCandidate)),
      end: Math.max(Number(targetOffset), Number(endCandidate)),
    };

    await page.waitForFunction(
      (range) => {
        const cfg = window.viewerApi.timeline.getLoopConfig?.();
        if (!cfg?.loopRange) return false;
        const startOffset = cfg.loopRange.start?.offset ?? null;
        const endOffset = cfg.loopRange.end?.offset ?? null;
        return (
          (startOffset === range.start && endOffset === range.end) ||
          (startOffset === range.end && endOffset === range.start)
        );
      },
      expectedRange,
      { timeout: 10_000 }
    );

    const fpsInput = advanced.locator('[data-testid="timeline-fps-input"]');
    await fpsInput.fill('');
    await fpsInput.type('48');
    await advanced.getByRole('button', { name: 'Set' }).click();

    await page.waitForFunction(() => {
      const timeline = window.viewerApi.timeline;
      if (!timeline) return false;
      const snapshot =
        typeof timeline.getPlaybackConfig === 'function'
          ? timeline.getPlaybackConfig()
          : typeof timeline.getPlaybackSnapshot === 'function'
            ? timeline.getPlaybackSnapshot()
            : null;
      return snapshot && Number(snapshot.defaultFps) === 48;
    });

    await advanced.locator('[data-testid="timeline-loop-clear"]').click();
    await page.waitForFunction(() => {
      const cfg = window.viewerApi.timeline.getLoopConfig?.();
      return !cfg?.loopRange;
    });
    console.log('[timeline-controls] loop cleared');

    await loopToggle.uncheck();
    await page.waitForFunction(() => {
      const cfg = window.viewerApi.timeline.getLoopConfig?.();
      return !cfg || cfg.loop !== true;
    });
    console.log('[timeline-controls] loop toggle unchecked');
  }

  await page.evaluate(() => window.viewerApi.stopSimulation());
  console.log('[timeline-controls] simulation stopped, test complete');
});
