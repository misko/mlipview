// Timeline rewind/resume acceptance test
import { test, expect } from './fixtures.js';
import {
  waitForTimelineBuffer,
  hoverTimelinePanel,
  selectTimelineOffset,
} from './utils/timeline.js';

test('timeline rewind preserves frame identity and resumes live stream', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  // Kick off an MD run long enough to populate the frame buffer.
  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 200, temperature: 1200 });
  });

  await waitForTimelineBuffer(page, 60, { timeout: 45_000 });

  const initialSize = await page.evaluate(() => window.viewerApi.timeline.bufferStats().size);

  // Reveal the timeline dock.
  await hoverTimelinePanel(page);

  const { firstOffset, secondOffset } = await page.evaluate(() => {
    const timeline = window.viewerApi?.timeline;
    if (!timeline || typeof timeline.getOffsets !== 'function') {
      return { firstOffset: null, secondOffset: null };
    }
    const offsets = timeline.getOffsets();
    if (!Array.isArray(offsets) || offsets.length === 0) {
      return { firstOffset: null, secondOffset: null };
    }
    const filtered = Array.from(new Set(offsets.filter((value) => value < -1)));
    if (filtered.length <= 1) {
      const fallback = Array.from(new Set(offsets));
      return {
        firstOffset: fallback[Math.min(1, fallback.length - 1)] ?? fallback[0] ?? null,
        secondOffset: fallback[Math.min(2, fallback.length - 1)] ?? fallback[0] ?? null,
      };
    }
    const firstIndex = Math.min(15, filtered.length - 1);
    const first = filtered[firstIndex];
    const second = filtered[Math.min(firstIndex + 1, filtered.length - 1)];
    return { firstOffset: first ?? null, secondOffset: second ?? null };
  });

  expect(firstOffset).not.toBeNull();
  expect(secondOffset).not.toBeNull();

  let primaryOffset = firstOffset;
  let secondaryOffset = secondOffset;
  if (primaryOffset === secondaryOffset) {
    const alternate = await page.evaluate((offset) => {
      const timeline = window.viewerApi?.timeline;
      if (!timeline || typeof timeline.getOffsets !== 'function') return null;
      const offsets = timeline.getOffsets();
      if (!Array.isArray(offsets) || offsets.length === 0) return null;
      const idx = offsets.indexOf(offset);
      if (idx === -1) {
        return offsets.find((value) => value !== offset) ?? null;
      }
      const next = offsets[idx + 1] ?? offsets[idx - 1] ?? null;
      return next && next !== offset ? next : offsets.find((value) => value !== offset) ?? null;
    }, primaryOffset);
    if (alternate != null) secondaryOffset = alternate;
  }
  expect(secondaryOffset).not.toBe(primaryOffset);

  // Select frame -15.
  await selectTimelineOffset(page, primaryOffset);

  await page.waitForFunction(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    return status && status.active && status.offset === primaryOffset;
  }, null, { timeout: 10_000 });

  const sigFirst = await page.evaluate(
    (offset) => window.viewerApi.timeline.getSignature(offset),
    primaryOffset
  );
  expect(sigFirst).toBeTruthy();

  // Move to frame -16.
  await selectTimelineOffset(page, secondaryOffset);

  await page.waitForFunction(
    (offset) => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.offset === offset;
    },
    secondaryOffset,
    { timeout: 10_000 }
  );

  // Return to frame -15 and confirm identical payload.
  await selectTimelineOffset(page, primaryOffset);

  await page.waitForFunction(
    (offset) => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.offset === offset;
    },
    primaryOffset,
    { timeout: 10_000 }
  );
  const sigSecond = await page.evaluate(
    (offset) => window.viewerApi.timeline.getSignature(offset),
    primaryOffset
  );
  expect(sigSecond).toBe(sigFirst);

  // Play forward from the selected frame and ensure we land on the latest frame.
  await page.click('[data-testid="timeline-play"]');
  await page.waitForFunction(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    return status && !status.active && status.mode === 'live' && status.offset === -1;
  }, null, { timeout: 30_000 });

  await page.waitForFunction(() => window.viewerApi.getMetrics().running === 'md', null, { timeout: 20_000 });

  const targetSize = initialSize + 30;
  await page.waitForFunction((expected) => {
    const stats = window.viewerApi.timeline.bufferStats();
    return stats.size >= expected;
  }, targetSize, { timeout: 60_000 });

  await page.evaluate(() => window.viewerApi.stopSimulation());
});
