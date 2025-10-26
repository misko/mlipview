import { test, expect } from './fixtures.js';

const HEALTH_URL = 'http://127.0.0.1:8000/serve/health';

async function getAtomCount(page) {
  return await page.evaluate(() => window.viewerApi?.state?.positions?.length || 0);
}

async function resetWsLog(page) {
  await page.evaluate(() => {
    window.__WS_HOOK_LOG = [];
  });
}

async function ensureWsHook(page) {
  await page.evaluate(() => {
    if (!window.__WS_HOOK_READY__) {
      window.__WS_HOOK_LOG = [];
      window.__WS_TEST_HOOK__ = (msg) => {
        try {
          const clone = msg ? JSON.parse(JSON.stringify(msg)) : msg;
          window.__WS_HOOK_LOG.push(clone);
        } catch {
          window.__WS_HOOK_LOG.push(msg);
        }
      };
      window.__WS_HOOK_READY__ = true;
    }
  });
}

async function getLastFullUpdate(page) {
  return await page.evaluate(() => {
    const log = window.__WS_HOOK_LOG || [];
    for (let i = log.length - 1; i >= 0; i--) {
      const msg = log[i];
      if (!msg) continue;
      const type = msg.type || msg.payload?.case || msg.payloadType;
      if (type === 'USER_INTERACTION' || type === 'userInteraction') {
        const payload = msg.payload?.value || {};
        const natoms =
          (typeof msg.natoms === 'number' && Number.isFinite(msg.natoms) && msg.natoms >= 0
            ? msg.natoms
            : undefined) ??
          payload.natoms ??
          payload.nAtoms ??
          payload.natom ??
          (payload.atomicNumbers &&
            typeof payload.atomicNumbers.indices === 'object' &&
            payload.atomicNumbers.indices.length) ??
          (typeof msg.positionsCount === 'number' ? msg.positionsCount : undefined);
        const fullUpdate =
          msg.fullUpdate === true ||
          (payload.fullUpdateFlag && payload.fullUpdateFlag.case === 'fullUpdate' && payload.fullUpdateFlag.value === true) ||
          payload.fullUpdate === true ||
          payload.full_update === true ||
          payload.fullupdate === true;
        if (typeof natoms === 'number' || fullUpdate) {
          return { natoms, fullUpdate };
        }
      }
    }
    return null;
  });
}

async function waitForStructureChanged(page, timeoutMs = 8000) {
  return await page.evaluate(
    ({ timeoutMs }) =>
      new Promise((resolve) => {
        try {
          const ws = window.__fairchem_ws__ || window.__WS_API__;
          if (!ws || typeof ws.onResult !== 'function') return resolve(false);
          const existing = window.__WS_HOOK_LOG;
          if (Array.isArray(existing)) {
            for (let i = existing.length - 1; i >= 0; i--) {
              const msg = existing[i];
              const text =
                msg?.message ||
                msg?.notice?.message ||
                (msg?.payload?.case === 'notice' ? msg.payload.value?.message : undefined);
              if (text === 'STRUCTURE_CHANGED') {
                return resolve(true);
              }
            }
          }
          const timeout = setTimeout(() => {
            off();
            resolve(false);
          }, timeoutMs);
          const off = ws.onResult((res) => {
            try {
              const msg =
                res?.message ||
                res?.notice?.message ||
                (res?.notice && res.notice.message);
              if (msg === 'STRUCTURE_CHANGED') {
                clearTimeout(timeout);
                off();
                resolve(true);
              }
            } catch {}
          });
        } catch {
          resolve(false);
        }
      }),
    { timeoutMs }
  );
}

async function screenFromAtomIndex(page, index) {
  return await page.evaluate(({ index }) => {
    const api = window.viewerApi;
    if (!api || !api.scene || !api.camera) return null;
    const state = api.state;
    if (!state || !Array.isArray(state.positions) || index >= state.positions.length) return null;
    const pos = state.positions[index];
    const engine =
      api.scene.getEngine && typeof api.scene.getEngine === 'function'
        ? api.scene.getEngine()
        : null;
    const canvas =
      engine && typeof engine.getRenderingCanvas === 'function'
        ? engine.getRenderingCanvas()
        : null;
    if (!canvas) return null;
    const rect = canvas.getBoundingClientRect();
    const world = new BABYLON.Vector3(pos.x, pos.y, pos.z);
    const vp = (() => {
      try {
        const viewport =
          api.camera.viewport ||
          (api.scene.activeCamera && api.scene.activeCamera.viewport) ||
          new BABYLON.Viewport(0, 0, 1, 1);
        if (!engine) return viewport;
        if (typeof viewport.toGlobal === 'function') {
          return viewport.toGlobal(engine.getRenderWidth(), engine.getRenderHeight());
        }
        const w = engine.getRenderWidth();
        const h = engine.getRenderHeight();
        return new BABYLON.Viewport(
          viewport.x * w,
          viewport.y * h,
          viewport.width * w,
          viewport.height * h
        );
      } catch {
        return new BABYLON.Viewport(0, 0, canvas.width || 0, canvas.height || 0);
      }
    })();
    const viewProj = api.scene.getTransformMatrix();
    const proj = BABYLON.Vector3.Project(world, BABYLON.Matrix.Identity(), viewProj, vp);
    if (!Number.isFinite(proj.x) || !Number.isFinite(proj.y)) return null;
    return {
      x: rect.left + proj.x,
      y: rect.top + proj.y,
      position: { x: pos.x, y: pos.y, z: pos.z },
    };
  }, { index });
}

async function worldFromClient(page, clientX, clientY) {
  return await page.evaluate(({ clientX, clientY }) => {
    const api = window.viewerApi;
    if (!api || !api.scene || !api.camera) return null;
    const engine =
      api.scene.getEngine && typeof api.scene.getEngine === 'function'
        ? api.scene.getEngine()
        : null;
    const canvas =
      engine && typeof engine.getRenderingCanvas === 'function'
        ? engine.getRenderingCanvas()
        : null;
    if (!canvas) return null;
    const rect = canvas.getBoundingClientRect();
    const px = clientX - rect.left;
    const py = clientY - rect.top;
    const camera = api.camera;
    const scene = api.scene;
    const ray = scene.createPickingRay(px, py, BABYLON.Matrix.Identity(), camera);
    const camPos = camera.position
      ? camera.position
      : camera.getFrontPosition
        ? camera.getFrontPosition(0)
        : { x: 0, y: 0, z: -10 };
    const target = camera.target || { x: 0, y: 0, z: 0 };
    const normal = new BABYLON.Vector3(
      target.x - camPos.x,
      target.y - camPos.y,
      target.z - camPos.z
    );
    if (normal.lengthSquared() < 1e-6) {
      normal.set(0, 1, 0);
    }
    normal.normalize();
    const planePoint = new BABYLON.Vector3(target.x, target.y, target.z);
    const denom = BABYLON.Vector3.Dot(normal, ray.direction);
    let point = planePoint;
    if (Math.abs(denom) >= 1e-6) {
      const t = BABYLON.Vector3.Dot(planePoint.subtract(ray.origin), normal) / denom;
      if (Number.isFinite(t)) {
        point = ray.origin.add(ray.direction.scale(t));
      }
    }
    return {
      inside:
        clientX >= rect.left &&
        clientX <= rect.right &&
        clientY >= rect.top &&
        clientY <= rect.bottom,
      point: { x: point.x, y: point.y, z: point.z },
    };
  }, { clientX, clientY });
}

test.describe('Periodic table atom editing', () => {
  test.beforeEach(async ({ page, loadViewerPage }) => {
    const health = await page.request.get(HEALTH_URL);
    if (!health.ok()) test.skip(true, 'Backend health unavailable');
    const meta = await health.json();
    if (!meta || String(meta.device || '').toLowerCase() !== 'cuda') {
      test.skip(true, 'Backend not reporting CUDA device');
    }
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });
    await ensureWsHook(page);
  });

  test('add/remove atoms via periodic table click, drag, and right-click during idle/MD/relax', async ({
    page,
  }) => {
    test.setTimeout(120_000);

    const canvasBox = await page.locator('#viewer').boundingBox();
    expect(canvasBox).not.toBeNull();
    const canvasCenter = {
      x: canvasBox.x + canvasBox.width / 2,
      y: canvasBox.y + canvasBox.height / 2,
    };

    // Ensure the periodic section is expanded so the grid is interactable
    const selectionHeader = page.locator('#section-selection .panel-header');
    if (await selectionHeader.isVisible()) {
      const expanded = await selectionHeader.getAttribute('aria-expanded');
      if (expanded !== 'true') {
        await selectionHeader.click();
        await page.waitForTimeout(150);
      }
    }

    const initialCount = await getAtomCount(page);

    // --- Idle: click add at origin ---
    await resetWsLog(page);
    await page.locator('#miniPeriodic td[data-symbol="C"]').click();
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.positions?.length || 0) === prev + 1,
      initialCount,
      { timeout: 8000 }
    );
    let count = await getAtomCount(page);
    expect(count).toBe(initialCount + 1);
    let lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(count);

    // Right-click remove newly added atom
    await resetWsLog(page);
    let coords = await screenFromAtomIndex(page, count - 1);
    expect(coords).toBeTruthy();
    await page.mouse.click(coords.x, coords.y, { button: 'right' });
    await page.waitForFunction(
      (target) => (window.viewerApi?.state?.positions?.length || 0) === target,
      initialCount,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(initialCount);

    // --- Idle: drag add onto canvas ---
    await resetWsLog(page);
    const dragCell = await page.locator('#miniPeriodic td[data-symbol="O"]').boundingBox();
    expect(dragCell).not.toBeNull();
    await page.mouse.move(dragCell.x + dragCell.width / 2, dragCell.y + dragCell.height / 2);
    await page.mouse.down();
    await page.mouse.move(canvasCenter.x, canvasCenter.y, { steps: 15 });
    await page.mouse.up();
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.positions?.length || 0) === prev + 1,
      initialCount,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount + 1);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(count);
    const dropWorld = await worldFromClient(page, canvasCenter.x, canvasCenter.y);
    const addedPos = await page.evaluate(() => {
      const st = window.viewerApi?.state;
      if (!st || !st.positions || !st.positions.length) return null;
      const p = st.positions[st.positions.length - 1];
      return { x: p.x, y: p.y, z: p.z };
    });
    expect(dropWorld?.point).toBeTruthy();
    const dist = Math.hypot(
      addedPos.x - dropWorld.point.x,
      addedPos.y - dropWorld.point.y,
      addedPos.z - dropWorld.point.z
    );
    expect(dist).toBeLessThan(3);

    // Remove drag-added atom to reset baseline
    await resetWsLog(page);
    coords = await screenFromAtomIndex(page, count - 1);
    await page.mouse.click(coords.x, coords.y, { button: 'right' });
    await page.waitForFunction(
      (target) => (window.viewerApi?.state?.positions?.length || 0) === target,
      initialCount,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(initialCount);

    // --- MD: add/remove during simulation ---
    await resetWsLog(page);
    await page.evaluate(() =>
      window.viewerApi?.startMDContinuous?.({ steps: 400, temperature: 900 })
    );
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'md',
      null,
      { timeout: 8000 }
    );

    const mdCell = await page.locator('#miniPeriodic td[data-symbol="N"]').boundingBox();
    expect(mdCell).not.toBeNull();
    await resetWsLog(page);
    const mdAddNotice = waitForStructureChanged(page);
    await page.mouse.move(mdCell.x + mdCell.width / 2, mdCell.y + mdCell.height / 2);
    await page.mouse.down();
    await page.mouse.move(canvasCenter.x, canvasCenter.y, { steps: 20 });
    await page.mouse.up();
    expect(await mdAddNotice).toBeTruthy();
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.positions?.length || 0) === prev + 1,
      initialCount,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount + 1);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(count);

    // Restart MD to ensure running before removal
    const prevMdForceVersion = await page.evaluate(() => window.viewerApi?.state?.forceCache?.version || 0);
    await page.evaluate(async () => {
      try {
        await window.viewerApi?.startMDContinuous?.({ steps: 200, temperature: 850 });
      } catch {}
    });
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'md',
      null,
      { timeout: 8000 }
    );
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.forceCache?.version || 0) > prev,
      prevMdForceVersion,
      { timeout: 8000 }
    );

    await resetWsLog(page);
    coords = await screenFromAtomIndex(page, count - 1);
    const mdRemoveNotice = waitForStructureChanged(page);
    await page.mouse.click(coords.x, coords.y, { button: 'right' });
    expect(await mdRemoveNotice).toBeTruthy();
    await page.waitForFunction(
      (target) => (window.viewerApi?.state?.positions?.length || 0) === target,
      initialCount,
      { timeout: 8000 }
    );
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'md',
      null,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(initialCount);
    await page.evaluate(() => window.viewerApi?.stopSimulation?.());
    await page.waitForFunction(
      () => !window.viewerApi?.getMetrics?.().running,
      null,
      { timeout: 8000 }
    );

    // --- Relax: add/remove during relax stream ---
    await page.evaluate(() =>
      window.viewerApi?.startRelaxContinuous?.({ maxSteps: 400, fmax: 0.05 })
    );
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'relax',
      null,
      { timeout: 8000 }
    );

    const relaxCell = await page.locator('#miniPeriodic td[data-symbol="S"]').boundingBox();
    expect(relaxCell).not.toBeNull();
    await resetWsLog(page);
    const relaxAddNotice = waitForStructureChanged(page);
    await page.mouse.move(relaxCell.x + relaxCell.width / 2, relaxCell.y + relaxCell.height / 2);
    await page.mouse.down();
    await page.mouse.move(canvasCenter.x, canvasCenter.y, { steps: 20 });
    await page.mouse.up();
    expect(await relaxAddNotice).toBeTruthy();
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.positions?.length || 0) === prev + 1,
      initialCount,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount + 1);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(count);

    const prevRelaxForceVersion = await page.evaluate(() => window.viewerApi?.state?.forceCache?.version || 0);
    await page.evaluate(() =>
      window.viewerApi?.startRelaxContinuous?.({ maxSteps: 200, fmax: 0.1 })
    );
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'relax',
      null,
      { timeout: 8000 }
    );
    await page.waitForFunction(
      (prev) => (window.viewerApi?.state?.forceCache?.version || 0) > prev,
      prevRelaxForceVersion,
      { timeout: 8000 }
    );

    await resetWsLog(page);
    coords = await screenFromAtomIndex(page, count - 1);
    const relaxRemoveNotice = waitForStructureChanged(page);
    await page.mouse.click(coords.x, coords.y, { button: 'right' });
    expect(await relaxRemoveNotice).toBeTruthy();
    await page.waitForFunction(
      (target) => (window.viewerApi?.state?.positions?.length || 0) === target,
      initialCount,
      { timeout: 8000 }
    );
    await page.waitForFunction(
      () => window.viewerApi?.getMetrics?.().running === 'relax',
      null,
      { timeout: 8000 }
    );
    count = await getAtomCount(page);
    expect(count).toBe(initialCount);
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(initialCount);
    await page.evaluate(() => window.viewerApi?.stopSimulation?.());
    await page.waitForFunction(
      () => !window.viewerApi?.getMetrics?.().running,
      null,
      { timeout: 8000 }
    );

    // Final full_update should reflect initial atom count
    lastUpdate = await getLastFullUpdate(page);
    expect(lastUpdate?.fullUpdate).toBeTruthy();
    expect(lastUpdate?.natoms).toBe(initialCount);
    const finalCount = await getAtomCount(page);
    expect(finalCount).toBe(initialCount);
  });
});
