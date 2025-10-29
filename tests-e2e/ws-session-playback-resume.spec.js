import { test, expect } from './fixtures.js';

const REQUIRED_LIVE_FRAMES = 18;
const MAX_LIVE_FRAMES = 170;
const TIMELINE_WARMUP_FRAMES = 48;
const SNAPSHOT_META = { kind: 'test', label: 'roy-md' };
const MD_WARMUP_STEPS = 120;

async function limitContinuousMd(page, maxSteps) {
  await page.evaluate((limit) => {
    const api = window.viewerApi;
    if (!api || typeof api.startMDContinuous !== 'function') {
      throw new Error('viewerApi.startMDContinuous unavailable');
    }
    if (api.__testMdLimit === limit) return;
    const original = api.startMDContinuous.bind(api);
    api.startMDContinuous = (opts = {}) => {
      const next = { ...(opts || {}) };
      const requested = Number.isFinite(next.steps) ? Math.floor(next.steps) : limit;
      next.steps = Math.max(1, Math.min(requested, limit));
      return original(next);
    };
    Object.defineProperty(api, '__testMdLimit', { value: limit, configurable: true });
  }, maxSteps);
}

async function warmupTimeline(page, { steps, temperature, minFrames }) {
  await page.evaluate(async ({ steps, temperature }) => {
    await window.viewerApi.startMDContinuous({ steps, temperature });
  }, { steps, temperature });

  const target = Math.max(10, minFrames | 0);
  let lastStats = null;
  for (let attempt = 0; attempt < 60; attempt++) {
    lastStats = await page.evaluate(() => {
      const timeline = window.viewerApi?.timeline;
      if (!timeline) return null;
      const stats = timeline.bufferStats?.();
      const offsets = timeline.getOffsets?.() || [];
      return {
        size: Number.isFinite(stats?.size) ? stats.size : offsets.length,
        offsets: offsets.length,
      };
    });
    if ((lastStats?.size ?? 0) >= target || (lastStats?.offsets ?? 0) >= target) {
      return;
    }
    await page.waitForTimeout(500);
  }
  throw new Error(`timeline warmup stuck: wanted >= ${target} frames, saw ${JSON.stringify(lastStats)}`);
}

async function captureSessionSnapshot(page, meta) {
  const json = await page.evaluate((metaInfo) => {
    const snap = window.viewerApi.session.captureSnapshot(metaInfo);
    return JSON.stringify(snap);
  }, meta);
  return JSON.parse(json);
}

async function ensureSimulationStopped(page) {
  await page.evaluate(() => {
    try { window.viewerApi.stopSimulation(); } catch {}
  });

  await page.waitForFunction(
    () => {
      const metrics = window.viewerApi?.getMetrics?.();
      return !metrics || !metrics.running;
    },
    { timeout: 8_000 },
  );
}

async function restoreSnapshot(page, snapshotJson) {
  await page.evaluate(async (json) => {
    const data = JSON.parse(json);
    await window.viewerApi.session.loadSnapshot(data);
  }, snapshotJson);

  await page.waitForFunction(
    () => {
      const offsets = window.viewerApi?.timeline?.getOffsets?.();
      return Array.isArray(offsets) && offsets.length > 0;
    },
    { timeout: 10_000 },
  );
}

async function choosePlaybackOffset(page, fallbackIndex) {
  return await page.evaluate(({ fallbackIdx }) => {
    const offsets = window.viewerApi?.timeline?.getOffsets?.();
    if (!offsets || !offsets.length) return null;
    if (offsets.includes(-5)) return -5;
    const idx = Math.min(fallbackIdx, offsets.length - 1);
    return offsets[idx];
  }, { fallbackIdx: fallbackIndex });
}

async function installLiveFrameCounter(page, { resumeSeq, maxFrames }) {
  await page.evaluate(({ resumeSeq, maxFrames }) => {
    const ws = window.__fairchem_ws__;
    if (!ws || typeof ws.onResult !== 'function') {
      throw new Error('ws.onResult unavailable');
    }
    if (window.__resumeListener) {
      try { window.__resumeListener(); } catch {}
    }
    window.__resumeFrameCount = 0;
    window.__resumeLastSeq = resumeSeq;
    window.__resumeMaxFrames = maxFrames;
    const stopStreams = () => {
      try { window.viewerApi.stopSimulation(); } catch {}
      try { ws.stopSimulation?.(); } catch {}
    };
    window.__resumeListener = ws.onResult((res) => {
      const seq = typeof res?.seq === 'number' ? res.seq : Number(res?.seq);
      if (!Number.isFinite(seq) || seq <= resumeSeq) return;
      window.__resumeLastSeq = seq;
      const next = (window.__resumeFrameCount || 0) + 1;
      window.__resumeFrameCount = next;
      if (next >= maxFrames) stopStreams();
    });
  }, { resumeSeq, maxFrames });
}

async function waitForTimelineLive(page, timeoutMs) {
  await page.waitForFunction(
    () => {
      const state = window.viewerApi?.timeline?.getState?.();
      return state && state.mode === 'live';
    },
    { timeout: timeoutMs },
  );
}

async function awaitLiveFrameWindow(page, { min, max, timeoutMs }) {
  const safeMin = Math.max(1, min | 0);
  const safeMax = Math.max(safeMin, max | 0);
  const deadline = Date.now() + Math.max(1_000, timeoutMs | 0);
  let lastCount = 0;
  while (Date.now() < deadline) {
    const { count } = await page.evaluate(() => ({
      count: window.__resumeFrameCount || 0,
    }));
    lastCount = count;
    if (count > safeMax) {
      throw new Error(`live frame count exceeded limit: ${count} > ${safeMax}`);
    }
    if (count >= safeMin) {
      return count;
    }
    await page.waitForTimeout(250);
  }
  throw new Error(`live frame count stuck below ${safeMin}; last=${lastCount}`);
}

async function teardownLiveCounter(page) {
  await page.evaluate(() => {
    try { window.viewerApi.stopSimulation(); } catch {}
    if (window.__resumeListener) {
      try { window.__resumeListener(); } catch {}
      window.__resumeListener = null;
    }
  });

  await page.waitForFunction(
    () => {
      const metrics = window.viewerApi?.getMetrics?.();
      return !metrics || !metrics.running;
    },
    { timeout: 8_000 },
  );
}

test.describe('session playback resume', () => {
  test('timeline playback restores and resumes MD streaming', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);

    await loadViewerPage({ query: { mol: 'molecules/roy.xyz' } });

    await limitContinuousMd(page, MAX_LIVE_FRAMES);

    await warmupTimeline(page, {
      steps: MD_WARMUP_STEPS,
      temperature: 650,
      minFrames: TIMELINE_WARMUP_FRAMES,
    });

    const snapshot = await captureSessionSnapshot(page, SNAPSHOT_META);
    const snapshotSeq = snapshot?.websocket?.seq ?? 0;
    const snapshotJson = JSON.stringify(snapshot);

    await ensureSimulationStopped(page);

    await restoreSnapshot(page, snapshotJson);

    const targetOffset = await choosePlaybackOffset(page, 4);
    expect(targetOffset).not.toBeNull();

    await installLiveFrameCounter(page, {
      resumeSeq: snapshotSeq,
      maxFrames: MAX_LIVE_FRAMES,
    });

    await page.evaluate((offset) => {
      window.viewerApi.timeline.play(offset);
    }, targetOffset);

    await waitForTimelineLive(page, 20_000);

    const liveCount = await awaitLiveFrameWindow(page, {
      min: REQUIRED_LIVE_FRAMES,
      max: MAX_LIVE_FRAMES,
      timeoutMs: 30_000,
    });

    expect(liveCount).toBeGreaterThanOrEqual(REQUIRED_LIVE_FRAMES);
    expect(liveCount).toBeLessThanOrEqual(MAX_LIVE_FRAMES);

    await teardownLiveCounter(page);
  });
});
