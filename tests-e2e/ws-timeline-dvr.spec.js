import { test, expect } from './fixtures.js';

async function waitForTimelineFrames(page, min, { timeout = 30_000 } = {}) {
  await page.waitForFunction(
    (required) => {
      const dbg = window.viewerApi?.timeline?.debug?.();
      return dbg && dbg.size >= required;
    },
    min,
    { timeout }
  );
}

async function seedMdFrames(page, { minFrames = 60, temperature = 1200 } = {}) {
  await page.evaluate((temp) => {
    try {
      window.viewerApi.startMDContinuous({ steps: 1200, temperature: temp });
    } catch { }
  }, temperature);
  await waitForTimelineFrames(page, minFrames);
}

test.describe('timeline DVR controls', () => {
  test('scrubbing via slider applies stored frame', async ({ page, loadViewerPage }) => {
    test.setTimeout(70_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 40, temperature: 1500 });
    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targetIndex = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      const target = Math.max(0, dbg.liveIndex - 8);
      return target;
    });
    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);
    await page.waitForTimeout(120);

    const result = await page.evaluate((index) => {
      const dbg = window.viewerApi.timeline.debug();
      const frame = window.viewerApi.timeline.getFrameTriples(index) || [];
      const pos = window.viewerApi.state.positions?.[0];
      return {
        mode: dbg.mode,
        playbackIndex: dbg.playbackIndex,
        liveIndex: dbg.liveIndex,
        pos0: pos ? [pos.x, pos.y, pos.z] : null,
        framePos0: frame[0] || null,
      };
    }, targetIndex);

    expect(result.mode).toBe('paused');
    expect(result.playbackIndex).toBe(targetIndex);
    expect(result.framePos0).toBeTruthy();
    expect(result.pos0).toBeTruthy();
    expect(result.pos0[0]).toBeCloseTo(result.framePos0[0], 8);
    expect(result.pos0[1]).toBeCloseTo(result.framePos0[1], 8);
    expect(result.pos0[2]).toBeCloseTo(result.framePos0[2], 8);
  });

  test('scrubbing multiple frames retains bond rendering', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 70, temperature: 1400 });
    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const indices = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      const live = dbg.liveIndex;
      const out = [];
      for (let i = 1; i <= 10; i++) {
        const idx = live - i;
        if (idx >= 0) out.push(idx);
      }
      return out;
    });

    for (const index of indices) {
      await slider.evaluate((node, value) => {
        node.value = String(value);
        node.dispatchEvent(new Event('input', { bubbles: true }));
        node.dispatchEvent(new Event('change', { bubbles: true }));
      }, index);
      await page.waitForFunction(
        (idx) => window.viewerApi.timeline.debug().playbackIndex === idx,
        index,
        { timeout: 5_000 }
      );
      const bondCount = await page.evaluate(() => window.viewerApi.state.bonds?.length ?? 0);
      expect(bondCount).toBeGreaterThan(0);
    }

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('playback runs at 20fps and updates button states', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 80, temperature: 1500 });
    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const playbackMeta = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      const target = Math.max(0, dbg.liveIndex - 12);
      return { liveIndex: dbg.liveIndex, targetIndex: target };
    });

    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, playbackMeta.targetIndex);
    await page.waitForFunction(
      (idx) => window.viewerApi.timeline.debug().playbackIndex === idx,
      playbackMeta.targetIndex,
      { timeout: 5_000 }
    );

    const frameCount = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return dbg.liveIndex - dbg.playbackIndex;
    });

    await page.evaluate(() => {
      window.__TIMELINE_PLAY_START = performance.now();
    });

    const playButton = page.getByRole('button', { name: 'Play from selected frame' });
    const pauseButton = page.getByRole('button', { name: 'Pause simulation' });
    const liveButton = page.getByRole('button', { name: 'Jump to live frame' });

    await playButton.click();

    await expect(playButton).toBeDisabled();
    await expect(pauseButton).toBeEnabled();
    await expect(liveButton).toBeEnabled();

    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'live',
      null,
      { timeout: 30_000 }
    );
    await page.evaluate(() => {
      window.__TIMELINE_PLAY_END = performance.now();
    });

    const elapsed = await page.evaluate(() => {
      return window.__TIMELINE_PLAY_END - window.__TIMELINE_PLAY_START;
    });

    expect(elapsed).toBeGreaterThanOrEqual(frameCount * 35);
    expect(elapsed).toBeLessThanOrEqual(frameCount * 150);

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('playback keeps slider aligned with applied frames', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 90, temperature: 1600 });
    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targetIndex = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return Math.max(0, dbg.liveIndex - 20);
    });

    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);
    await page.waitForFunction(
      (idx) => window.viewerApi.timeline.debug().playbackIndex === idx,
      targetIndex,
      { timeout: 5_000 }
    );

    const hookReady = await page.evaluate(() => {
      try { window.__TIMELINE_INDEX_OFF__?.(); } catch { window.__TIMELINE_INDEX_OFF__ = null; }
      window.__TIMELINE_INDEX_SAMPLES__ = [];
      const controller = window.viewerApi.timeline.controller;
      if (!controller?.on) return false;
      const off = controller.on('index', ({ index }) => {
        try {
          const sliderEl = document.querySelector('.mlip-timeline-slider');
          const sliderVal = sliderEl ? Number(sliderEl.value) : null;
          window.__TIMELINE_INDEX_SAMPLES__.push({
            index,
            slider: sliderVal,
            ts: performance.now(),
          });
        } catch {
          window.__TIMELINE_INDEX_SAMPLES__.push({ index, slider: null, ts: performance.now() });
        }
      });
      window.__TIMELINE_INDEX_OFF__ = off;
      return true;
    });
    expect(hookReady).toBe(true);

    await page.getByRole('button', { name: 'Play from selected frame' }).click();
    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'live',
      null,
      { timeout: 30_000 }
    );

    const playbackSummary = await page.evaluate(() => {
      const data = {
        samples: window.__TIMELINE_INDEX_SAMPLES__
          ? window.__TIMELINE_INDEX_SAMPLES__.map((entry) => ({
              index: entry.index,
              slider: entry.slider,
              ts: entry.ts,
            }))
          : [],
        liveIndex: window.viewerApi.timeline.debug().liveIndex,
      };
      try { window.__TIMELINE_INDEX_OFF__?.(); } catch { }
      window.__TIMELINE_INDEX_OFF__ = null;
      window.__TIMELINE_INDEX_SAMPLES__ = [];
      return data;
    });

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    const { samples, liveIndex } = playbackSummary;
    expect(samples.length).toBeGreaterThan(3);
    const progressed = samples.filter(
      (sample) => typeof sample.index === 'number' && sample.index > targetIndex
    );
    expect(progressed.length).toBeGreaterThan(0);
    for (const sample of progressed) {
      expect(sample.slider).toBe(sample.index);
    }
    const lastSample = samples[samples.length - 1];
    expect(lastSample.index).toBe(liveIndex);
    expect(lastSample.slider).toBe(liveIndex);
  });

  test('play button while live keeps streaming', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 40, temperature: 1400 });
    await page.hover('.mlip-timeline-rail');

    const initial = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return {
        liveIndex: dbg.liveIndex,
        mode: dbg.mode,
      };
    });

    expect(initial.mode).toBe('live');

    const playButton = page.getByRole('button', { name: 'Play from selected frame' });
    await playButton.click();

    await page.waitForFunction(
      (startLiveIdx) => {
        const dbg = window.viewerApi.timeline.debug();
        return dbg.mode === 'live' && dbg.liveIndex > startLiveIdx + 2;
      },
      initial.liveIndex,
      { timeout: 20_000 }
    );

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('scrubbing increments interaction counter and forces pause', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 45, temperature: 1500 });

    const initialInfo = await page.evaluate(() => ({
      uic: window.viewerApi.getVersionInfo().userInteractionVersion,
    }));

    await page.evaluate(() => {
      window.__TIMELINE_TEST_ACTIONS__ = [];
      window.__WS_TEST_HOOK__ = (msg) => {
        try {
          window.__TIMELINE_TEST_ACTIONS__.push(msg);
        } catch { }
      };
    });

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();
    const targetIndex = await page.evaluate(() => Math.max(0, window.viewerApi.timeline.debug().liveIndex - 6));
    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);
    await page.waitForFunction(() => window.viewerApi.timeline.debug().mode === 'paused', null, { timeout: 10_000 });

    const pauseButton = page.getByRole('button', { name: 'Pause simulation' });
    await expect(pauseButton).toBeDisabled();

    const stats = await page.evaluate((initial) => {
      const actions = window.__TIMELINE_TEST_ACTIONS__ || [];
      const stopActions = actions.filter((msg) => msg?.type === 'STOP_SIMULATION').length;
      const userActions = actions.filter((msg) => msg?.type === 'USER_INTERACTION').length;
      window.__WS_TEST_HOOK__ = null;
      return {
        mode: window.viewerApi.timeline.debug().mode,
        uic: window.viewerApi.getVersionInfo().userInteractionVersion,
        stopActions,
        userActions,
      };
    }, initialInfo);

    expect(stats.mode).toBe('paused');
    expect(stats.uic).toBe(initialInfo.uic + 1);
    expect(stats.stopActions).toBeGreaterThan(0);
    expect(stats.userActions).toBe(0);

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('playback catches up and restarts streaming', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 80, temperature: 1600 });

    await page.evaluate(() => {
      window.__TIMELINE_TEST_ACTIONS__ = [];
      window.__WS_TEST_HOOK__ = (msg) => {
        try {
          window.__TIMELINE_TEST_ACTIONS__.push(msg);
        } catch { }
      };
    });

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targetIndex = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return Math.max(0, dbg.liveIndex - 18);
    });
    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);
    await page.waitForTimeout(120);

    await page.getByRole('button', { name: 'Play from selected frame' }).click();
    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'live',
      null,
      { timeout: 25_000 }
    );

    const summary = await page.evaluate(() => {
      const msgs = window.__TIMELINE_TEST_ACTIONS__ || [];
      window.__WS_TEST_HOOK__ = null;
      const starts = [];
      const fullUpdates = [];
      msgs.forEach((msg, idx) => {
        if (msg?.type === 'START_SIMULATION') {
          starts.push({ idx });
        } else if (msg?.type === 'USER_INTERACTION' && msg?.fullUpdate) {
          fullUpdates.push({ idx });
        }
      });
      return {
        mode: window.viewerApi.timeline.debug().mode,
        startsCount: starts.length,
        fullUpdatesCount: fullUpdates.length,
        hasFullUpdateBeforeStart:
          starts.length > 0
            ? fullUpdates.some((entry) => entry.idx < starts[0].idx)
            : false,
      };
    });

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    expect(summary.mode).toBe('live');
    expect(summary.fullUpdatesCount).toBeGreaterThan(0);
    expect(summary.startsCount).toBeGreaterThan(0);
    expect(summary.hasFullUpdateBeforeStart).toBe(true);
  });

  test('timeline slider displays relative tick labels', async ({ page, loadViewerPage }) => {
    test.setTimeout(70_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 55, temperature: 1500 });
    await page.hover('.mlip-timeline-rail');

    const labels = await page.evaluate(() => {
      const list = document.getElementById('mlip-timeline-ticks');
      if (!list) return null;
      const options = Array.from(list.querySelectorAll('option'));
      const dbg = window.viewerApi.timeline.debug();
      const liveIdx = dbg.liveIndex;
      const liveOpt = options.find((opt) => Number(opt.value) === liveIdx);
      const minusOneOpt = options.find((opt) => Number(opt.value) === liveIdx - 1);
      const oldestOpt = options[0];
      return {
        liveLabel: liveOpt?.label ?? null,
        minusOneLabel: minusOneOpt?.label ?? null,
        oldestLabel: oldestOpt?.label ?? null,
        optionCount: options.length,
      };
    });

    expect(labels?.optionCount).toBeGreaterThan(1);
    expect(labels?.liveLabel).toBe('0 (live)');
    expect(labels?.minusOneLabel).toBe('-1');
    expect(labels?.oldestLabel).toMatch(/-\d+/);

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('pause stops simulation until Live resumes', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 50, temperature: 1400 });

    await page.evaluate(() => {
      window.__TIMELINE_TEST_ACTIONS__ = [];
      window.__WS_TEST_HOOK__ = (msg) => {
        try {
          window.__TIMELINE_TEST_ACTIONS__.push(msg);
        } catch { }
      };
    });

    await page.hover('.mlip-timeline-rail');
    await page.getByRole('button', { name: 'Pause simulation' }).click();
    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'paused',
      null,
      { timeout: 10_000 }
    );

    const pauseState = await page.evaluate(() => ({
      mode: window.viewerApi.timeline.debug().mode,
      running: window.viewerApi.getMetrics()?.running ?? null,
      stopActions: (window.__TIMELINE_TEST_ACTIONS__ || []).filter(
        (msg) => msg?.type === 'STOP_SIMULATION'
      ).length,
    }));

    expect(pauseState.mode).toBe('paused');
    expect(pauseState.running).toBeNull();
    expect(pauseState.stopActions).toBeGreaterThan(0);

    await page.getByRole('button', { name: 'Jump to live frame' }).click();
    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'live',
      null,
      { timeout: 20_000 }
    );

    const resumeState = await page.evaluate(() => {
      const starts = (window.__TIMELINE_TEST_ACTIONS__ || []).filter(
        (msg) => msg?.type === 'START_SIMULATION'
      ).length;
      window.__WS_TEST_HOOK__ = null;
      return {
        mode: window.viewerApi.timeline.debug().mode,
        running: window.viewerApi.getMetrics()?.running,
        startActions: starts,
      };
    });

    expect(resumeState.mode).toBe('live');
    expect(resumeState.running).toBe('md');
    expect(resumeState.startActions).toBeGreaterThan(0);
    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('timeline buffer caps at 500 frames', async ({ page, loadViewerPage }) => {
    test.setTimeout(50_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await waitForTimelineFrames(page, 1, { timeout: 10_000 });

    const capResult = await page.evaluate(() => {
      const basePositions = window.viewerApi.timeline.getFrameTriples(
        window.viewerApi.timeline.debug().liveIndex
      ) || [];
      const sample = basePositions.length ? basePositions : [[0, 0, 0]];
      for (let i = 0; i < 620; i++) {
        const offset = i * 0.001;
        const positions = sample.map(([x, y, z]) => [x + offset, y, z]);
        window.viewerApi.timeline.injectTestFrame({
          kind: 'test',
          positions,
          energy: offset,
          simStep: i,
        });
      }
      const dbg = window.viewerApi.timeline.debug();
      const oldest = window.viewerApi.timeline.getFrameMeta(0);
      const newest = window.viewerApi.timeline.getFrameMeta(dbg.liveIndex);
      return {
        size: dbg.size,
        oldestId: oldest?.id ?? null,
        newestId: newest?.id ?? null,
      };
    });

    expect(capResult.size).toBe(500);
    expect(capResult.newestId).toBeGreaterThan(capResult.oldestId);
    expect(capResult.oldestId).toBeGreaterThan(100);
  });

  test('rail is keyboard accessible', async ({ page, loadViewerPage }) => {
    test.setTimeout(50_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 30, temperature: 1300 });
    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    const rail = page.locator('.mlip-timeline-rail');
    await rail.focus();
    await expect(rail).toHaveAttribute('data-expanded', 'true');

    const slider = page.locator('.mlip-timeline-slider');
    await slider.focus();
    const before = await page.evaluate(() => window.viewerApi.timeline.debug().playbackIndex);
    await page.keyboard.press('ArrowLeft');
    await page.waitForTimeout(80);
    const after = await page.evaluate(() => window.viewerApi.timeline.debug().playbackIndex);
    expect(after).toBeLessThanOrEqual(before);
  });

  test('live button resumes streaming after timeline scrub', async ({ page, loadViewerPage }) => {
    test.setTimeout(90_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 60, temperature: 1500 });

    await page.evaluate(() => {
      window.__TIMELINE_TEST_ACTIONS__ = [];
      window.__WS_TEST_HOOK__ = (msg) => {
        try {
          window.__TIMELINE_TEST_ACTIONS__.push(msg);
        } catch { }
      };
    });

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();
    const targetIndex = await page.evaluate(() => Math.max(0, window.viewerApi.timeline.debug().liveIndex - 5));
    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);
    await page.waitForFunction(
      (idx) => window.viewerApi.timeline.debug().playbackIndex === idx && window.viewerApi.timeline.debug().mode === 'paused',
      targetIndex,
      { timeout: 10_000 }
    );

    await page.getByRole('button', { name: 'Jump to live frame' }).click();
    await page.waitForFunction(
      () => window.viewerApi.timeline.debug().mode === 'live',
      null,
      { timeout: 20_000 }
    );

    const startCount = await page.evaluate(() => {
      const actions = window.__TIMELINE_TEST_ACTIONS__ || [];
      window.__WS_TEST_HOOK__ = null;
      return actions.filter((msg) => msg?.type === 'START_SIMULATION').length;
    });

    expect(startCount).toBeGreaterThan(0);

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });
  });

  test('timeline frame switch latency remains under 1s', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 80, temperature: 1500 });
    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    await page.hover('.mlip-timeline-rail');

    const indices = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      const live = dbg.liveIndex;
      const out = [];
      for (let i = 1; i <= 5; i++) {
        const idx = live - i;
        if (idx >= 0) out.push(idx);
      }
      return out;
    });

    expect(indices.length).toBeGreaterThan(0);

    const durations = [];
    for (const index of indices) {
      const duration = await page.evaluate(async (target) => {
        const slider = document.querySelector('.mlip-timeline-slider');
        if (!slider) throw new Error('timeline slider missing');
        const start = performance.now();
        slider.value = String(target);
        slider.dispatchEvent(new Event('input', { bubbles: true }));
        slider.dispatchEvent(new Event('change', { bubbles: true }));
        await new Promise((resolve, reject) => {
          const deadline = performance.now() + 2000;
          const check = () => {
            const dbg = window.viewerApi?.timeline?.debug?.();
            if (dbg && dbg.playbackIndex === target && dbg.mode === 'paused') {
              resolve();
              return;
            }
            if (performance.now() > deadline) {
              reject(new Error('timeline scrub timed out'));
              return;
            }
            requestAnimationFrame(check);
          };
          check();
        });
        const end = performance.now();
        return end - start;
      }, index);
      durations.push(duration);
    }

    const maxDuration = Math.max(...durations);
    expect(maxDuration).toBeLessThan(1000);
  });

  test('rapid sequential scrubs stay under 200ms', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 80, temperature: 1500 });
    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch { }
    });

    await page.evaluate(() => {
      window.__TIMELINE_APPLY_TRACE = [];
      const orig = console.log.bind(console);
      console.log = (...args) => {
        try {
          if (args[0] === '[timeline][apply-frame]' && args[1] && typeof args[1].durationMs === 'number') {
            window.__TIMELINE_APPLY_TRACE.push(args[1]);
          }
        } catch {}
        orig(...args);
      };
      window.__TIMELINE_LOG_ORIG = orig;
    });

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targets = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      const live = dbg.liveIndex;
      return [live - 4, live - 5, live - 6].filter((v) => v >= 0);
    });

    try {
      for (const target of targets) {
        await slider.evaluate((node, value) => {
          node.value = String(value);
          node.dispatchEvent(new Event('input', { bubbles: true }));
          node.dispatchEvent(new Event('change', { bubbles: true }));
        }, target);
        await page.waitForFunction(
          (idx) => window.viewerApi.timeline.debug().playbackIndex === idx,
          target,
          { timeout: 10_000 }
        );
      }
    } finally {
      await page.evaluate(() => {
        if (window.__TIMELINE_LOG_ORIG) {
          console.log = window.__TIMELINE_LOG_ORIG;
          delete window.__TIMELINE_LOG_ORIG;
        }
      });
    }

    const durations = await page.evaluate((count) => {
      const out = (window.__TIMELINE_APPLY_TRACE || []).map((entry) => entry.durationMs);
      return out.slice(-count);
    }, targets.length);

    expect(durations.length).toBeGreaterThan(0);
    const slowest = Math.max(...durations);
    expect(slowest).toBeLessThan(200);
  });

  test('default lattice stays disabled after timeline playback', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      disableAutoMd: true,
      testMode: true,
    });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    await seedMdFrames(page, { minFrames: 70, temperature: 1400 });

    const liveState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cellEnabled: !!(st.cell && st.cell.enabled),
        mode: window.viewerApi.timeline.debug().mode,
      };
    });
    expect(liveState.mode).toBe('live');
    expect(liveState.showCell).toBe(false);
    expect(liveState.cellEnabled).toBe(false);

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targetIndex = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return Math.max(0, dbg.liveIndex - 12);
    });

    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);

    await page.waitForFunction(() => window.viewerApi.timeline.debug().mode === 'paused', null, {
      timeout: 15_000,
    });

    const pausedState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cellEnabled: !!(st.cell && st.cell.enabled),
      };
    });
    expect(pausedState.showCell).toBe(false);
    expect(pausedState.cellEnabled).toBe(false);

    const playButton = page.getByRole('button', { name: 'Play from selected frame' });
    await playButton.click();

    await page.waitForFunction(() => window.viewerApi.timeline.debug().mode === 'live', null, {
      timeout: 30_000,
    });

    const resumedState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cellEnabled: !!(st.cell && st.cell.enabled),
      };
    });
    expect(resumedState.showCell).toBe(false);
    expect(resumedState.cellEnabled).toBe(false);

    await page.evaluate(() => {
      try { window.viewerApi.stopSimulation(); } catch {}
    });
  });

  test('cell visibility persists through timeline pause and resume', async ({ page, loadViewerPage }) => {
    test.setTimeout(120_000);
    await loadViewerPage({
      query: { mol: 'molecules/water.xyz', autoMD: 0 },
      disableAutoMd: true,
      testMode: true,
    });
    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');

    const initialCell = await page.evaluate(async () => {
      const api = window.viewerApi;
      const cell = [
        [8.5, 0, 0],
        [0, 8.5, 0],
        [0, 0, 8.5],
      ];
      const elements = api.state.elements.map((e) => e.symbol || e.sym || e.S || e);
      const positions = api.state.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
      await api.applyFullSnapshot({ elements, positions, cell });
      api.state.showCell = true;
      api.state.markCellChanged?.();
      return {
        showCell: api.state.showCell,
        cell: [
          [api.state.cell.a.x, api.state.cell.a.y, api.state.cell.a.z],
          [api.state.cell.b.x, api.state.cell.b.y, api.state.cell.b.z],
          [api.state.cell.c.x, api.state.cell.c.y, api.state.cell.c.z],
        ],
      };
    });

    expect(initialCell.showCell).toBe(true);

    await seedMdFrames(page, { minFrames: 70, temperature: 1400 });

    const liveState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cell: st.cell
          ? [
              [st.cell.a.x, st.cell.a.y, st.cell.a.z],
              [st.cell.b.x, st.cell.b.y, st.cell.b.z],
              [st.cell.c.x, st.cell.c.y, st.cell.c.z],
            ]
          : null,
        mode: window.viewerApi.timeline.debug().mode,
      };
    });

    expect(liveState.mode).toBe('live');
    expect(liveState.showCell).toBe(true);
    expect(liveState.cell).toEqual(initialCell.cell);

    await page.hover('.mlip-timeline-rail');
    const slider = page.locator('.mlip-timeline-slider');
    await expect(slider).toBeEnabled();

    const targetIndex = await page.evaluate(() => {
      const dbg = window.viewerApi.timeline.debug();
      return Math.max(0, dbg.liveIndex - 10);
    });

    await slider.evaluate((node, value) => {
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetIndex);

    await page.waitForFunction(() => window.viewerApi.timeline.debug().mode === 'paused', null, {
      timeout: 15_000,
    });

    const pausedState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cell: st.cell
          ? [
              [st.cell.a.x, st.cell.a.y, st.cell.a.z],
              [st.cell.b.x, st.cell.b.y, st.cell.b.z],
              [st.cell.c.x, st.cell.c.y, st.cell.c.z],
            ]
          : null,
      };
    });

    expect(pausedState.showCell).toBe(true);
    expect(pausedState.cell).toEqual(liveState.cell);

    const playButton = page.getByRole('button', { name: 'Play from selected frame' });
    await playButton.click();

    await page.waitForFunction(() => window.viewerApi.timeline.debug().mode === 'live', null, {
      timeout: 30_000,
    });

    const resumedState = await page.evaluate(() => {
      const st = window.viewerApi.state;
      return {
        showCell: !!st.showCell,
        cell: st.cell
          ? [
              [st.cell.a.x, st.cell.a.y, st.cell.a.z],
              [st.cell.b.x, st.cell.b.y, st.cell.b.z],
              [st.cell.c.x, st.cell.c.y, st.cell.c.z],
            ]
          : null,
      };
    });

    expect(resumedState.showCell).toBe(true);
    expect(resumedState.cell).toEqual(liveState.cell);

    await page.evaluate(() => {
      try {
        window.viewerApi.stopSimulation();
      } catch {}
    });
  });
});
