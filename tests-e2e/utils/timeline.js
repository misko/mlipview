import { expect } from '@playwright/test';

export async function waitForTimelineBuffer(
  page,
  minSize,
  { timeout = 45_000, label = 'waitForTimelineBuffer' } = {}
) {
  const start = Date.now();
  console.log(`[timeline-utils] ${label}: waiting for buffer size >= ${minSize}`);
  await page.waitForFunction(
    (target) => {
      const stats = window.viewerApi?.timeline?.bufferStats?.();
      return stats && stats.size >= target;
    },
    minSize,
    { timeout }
  );
  const elapsed = Date.now() - start;
  const stats = await page.evaluate(() => window.viewerApi?.timeline?.bufferStats?.());
  console.log(
    `[timeline-utils] ${label}: buffer ready in ${elapsed}ms (size=${stats?.size ?? 'n/a'})`
  );
}

export async function hoverTimelinePanel(page, label = 'hoverTimelinePanel') {
  console.log(`[timeline-utils] ${label}: hover timeline hitbox`);
  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  await expect(page.locator('[data-testid="timeline-panel"]')).toBeVisible();
  console.log(`[timeline-utils] ${label}: panel visible`);
}

export async function ensureTimelineSliderReady(
  page,
  { timeout = 20_000, label = 'ensureTimelineSliderReady' } = {}
) {
  const start = Date.now();
  console.log(`[timeline-utils] ${label}: waiting for slider to publish offsets`);
  await page.waitForFunction(
    () => {
      const slider = document.querySelector('[data-testid="timeline-slider"]');
      if (!slider) return false;
      const offsets = window.viewerApi?.timeline?.getOffsets?.();
      return Array.isArray(offsets) && offsets.length > 1;
    },
    null,
    { timeout }
  );
  console.log(
    `[timeline-utils] ${label}: slider ready after ${Date.now() - start}ms`,
    await page.evaluate(() => ({
      min: document.querySelector('[data-testid="timeline-slider"]')?.min ?? null,
      max: document.querySelector('[data-testid="timeline-slider"]')?.max ?? null,
      length: window.viewerApi?.timeline?.getOffsets?.()?.length ?? null,
    }))
  );
}

export async function waitForTimelineState(
  page,
  predicate,
  arg,
  { timeout = 20_000, label = 'waitForTimelineState' } = {}
) {
  const start = Date.now();
  console.log(`[timeline-utils] ${label}: waitForFunction start (timeout=${timeout})`);
  await page.waitForFunction(predicate, arg, { timeout });
  const elapsed = Date.now() - start;
  const state = await page.evaluate(() => window.viewerApi?.timeline?.getState?.());
  console.log(`[timeline-utils] ${label}: condition satisfied in ${elapsed}ms`, state);
}

export async function selectTimelineOffset(
  page,
  targetOffset,
  { timeout = 15_000, label = 'selectTimelineOffset' } = {}
) {
  console.log(`[timeline-utils] ${label}: selecting offset ${targetOffset}`);
  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.waitFor({ state: 'visible' });

  const selectionResult = await page.evaluate(async (offset) => {
    const timeline = window.viewerApi?.timeline;
    if (!timeline || typeof timeline.select !== 'function') return null;
    const result = timeline.select(offset);
    if (typeof timeline.ensurePaused === 'function') {
      try {
        await timeline.ensurePaused();
      } catch (err) {
        console.warn('[timeline-utils] ensurePaused failed', err);
      }
    }
    const status =
      (typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    return { result, status };
  }, targetOffset);
  console.log(`[timeline-utils] ${label}: select() returned`, selectionResult?.result ?? selectionResult);
  console.log(`[timeline-utils] ${label}: immediate status`, selectionResult?.status ?? null);

  await page.waitForFunction(
    (offset) => {
      const timeline = window.viewerApi?.timeline;
      if (!timeline) return false;
      const status =
        (typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.offset === offset;
    },
    targetOffset,
    { timeout }
  );
  console.log(
    `[timeline-utils] ${label}: selection applied`,
    await page.evaluate(() => {
      const timeline = window.viewerApi?.timeline;
      if (!timeline) return null;
      return (
        (typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (typeof timeline.getState === 'function' && timeline.getState()) ||
        null
      );
    })
  );
}

export async function computeTimelineOffset(page, selector) {
  return await page.evaluate((options) => {
    const timeline = window.viewerApi?.timeline;
    if (!timeline || typeof timeline.getOffsets !== 'function') return null;
    const offsets = timeline.getOffsets();
    if (!Array.isArray(offsets) || offsets.length === 0) return null;
    if (options?.mode === 'middle') {
      const filtered = offsets.filter((value) => value < -1);
      if (filtered.length === 0) return offsets[0] ?? -1;
      return filtered[Math.floor(filtered.length / 2)];
    }
    const target = Number(options?.offset ?? -1);
    const match = offsets.find((value) => value <= target);
    return match ?? offsets[offsets.length - 1] ?? null;
  }, selector || null);
}
