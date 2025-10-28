import { test, expect } from './fixtures.js';

const BACKEND_URL = 'http://127.0.0.1:8000';

async function loadViewer(page, loadViewerPage, ensureWsReady) {
  await loadViewerPage({
    server: BACKEND_URL,
    query: { mol: 'molecules/water.xyz', autoMD: 0 },
    testMode: true,
    disableAutoMd: true,
  });
  expect(await ensureWsReady()).toBeTruthy();
}

test.describe('Bond rotation isolation', () => {
  test('rotation group remains constant even when bonds change', async ({ page, loadViewerPage, ensureWsReady }) => {
    await loadViewer(page, loadViewerPage, ensureWsReady);

    const initial = await page.evaluate(() => {
      const firstBond = window.viewerApi.state.bonds?.[0];
      if (!firstBond) throw new Error('no bonds available');
      const ref = {
        i: firstBond.i,
        j: firstBond.j,
        key: firstBond.key ?? `${firstBond.i}-${firstBond.j}`,
        index: firstBond.index ?? 0,
      };
      window.viewerApi.debugSelectBond(ref);
      const sel = window.viewerApi.debugGetSelection();
      const group = sel?.data?.rotationGroup;
      if (!group) throw new Error('rotation group missing');
      return { sideAtoms: [...group.sideAtoms], anchor: group.anchor, movingRoot: group.movingRoot };
    });

    expect(initial.sideAtoms.length).toBeGreaterThan(0);

    await page.evaluate(() => {
      const bonds = window.viewerApi.state.bonds;
      bonds.push({ i: 0, j: 1, key: '0-1-extra', index: bonds.length, opacity: 1, order: 1 });
    });

    const afterRotation = await page.evaluate(() => {
      window.viewerApi.manipulation.rotateBond(Math.PI / 6);
      const dbg = window.viewerApi.manipulation._debug.getLastRotation();
      return dbg ? { sideAtoms: [...dbg.sideAtoms], anchor: dbg.anchor, movingRoot: dbg.movingRoot } : null;
    });

    expect(afterRotation).not.toBeNull();
    expect(afterRotation.sideAtoms).toEqual(initial.sideAtoms);
    expect(afterRotation.anchor).toBe(initial.anchor);
    expect(afterRotation.movingRoot).toBe(initial.movingRoot);
  });

  test('frame updates do not overwrite rotating atoms', async ({ page, loadViewerPage, ensureWsReady }) => {
    await loadViewer(page, loadViewerPage, ensureWsReady);

    const rotationData = await page.evaluate(() => {
      window.viewerApi.debugSelectBond({ i: 0, j: 1 });
      const sel = window.viewerApi.debugGetSelection();
      const group = sel?.data?.rotationGroup;
      if (!group) throw new Error('rotation group missing');
      window.viewerApi.manipulation.rotateBond(Math.PI / 8);
      const positions = window.viewerApi.state.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
      return {
        sideAtoms: [...group.sideAtoms],
        positions,
      };
    });

    expect(rotationData.sideAtoms.length).toBeGreaterThan(0);

    const debugUIC = await page.evaluate(({ rotationData }) => {
      const ws = window.__fairchem_ws__;
      ws?.pauseIncoming?.(120);
      const beforeState = ws?.getState?.() || {};
      const positions = rotationData.positions.map((p, idx) => {
        if (rotationData.sideAtoms.includes(idx)) {
          return [p.x + 10, p.y + 10, p.z + 10];
        }
        return [p.x + 1, p.y + 1, p.z + 1];
      });
      const injectUIC = Number.isFinite(beforeState.userInteractionCount) ? beforeState.userInteractionCount : 0;
      ws.injectTestResult({ seq: 999, client_seq: 0, userInteractionCount: injectUIC, simStep: 0, positions });
      const afterState = ws?.getState?.() || {};
      console.log('[bond-rotation:test][uic]', beforeState.userInteractionCount, afterState.userInteractionCount);
      return { before: beforeState.userInteractionCount ?? null, after: afterState.userInteractionCount ?? null };
    }, { rotationData });

    await page.waitForTimeout(50);

    const { positions: finalState, deltas, state } = await page.evaluate(({ rotationData }) => {
      const ws = window.__fairchem_ws__;
      const beforeState = ws?.getState?.() || {};
      const positions = window.viewerApi.state.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
      const deltas = positions.map((after, idx) => {
        const before = rotationData.positions[idx];
        const dx = after.x - before.x;
        const dy = after.y - before.y;
        const dz = after.z - before.z;
        return { idx, dx, dy, dz, dist: Math.hypot(dx, dy, dz), side: rotationData.sideAtoms.includes(idx) };
      });
      const afterState = ws?.getState?.() || {};
      const side = new Set(rotationData.sideAtoms);
      const formatted = delta =>
        `idx=${delta.idx}` +
        ` side=${side.has(delta.idx)}` +
        ` dx=${delta.dx.toFixed(6)}` +
        ` dy=${delta.dy.toFixed(6)}` +
        ` dz=${delta.dz.toFixed(6)}` +
        ` dist=${delta.dist.toFixed(6)}`;
      const grouped = {
        side: deltas.filter((d) => side.has(d.idx)).map(formatted),
        other: deltas.filter((d) => !side.has(d.idx)).map(formatted),
      };
      console.log('[bond-rotation][delta-report]', grouped);
      const top = deltas
        .map((d) => ({ ...d, abs: Math.max(Math.abs(d.dx), Math.abs(d.dy), Math.abs(d.dz)) }))
        .sort((a, b) => b.abs - a.abs)
        .slice(0, 10)
        .map((d) => formatted(d));
      console.log('[bond-rotation][delta-top10]', top);
      return { positions, deltas, state: { before: beforeState, after: afterState }, grouped, top };
    }, { rotationData });
    console.log('[bond-rotation][delta-summary]', { top: (state?.top ?? []), state });

    for (const idx of rotationData.sideAtoms) {
      const before = rotationData.positions[idx];
      const after = finalState[idx];
      expect(after.x).toBeCloseTo(before.x, 6);
      expect(after.y).toBeCloseTo(before.y, 6);
      expect(after.z).toBeCloseTo(before.z, 6);
    }

    const otherIdx = finalState.findIndex((_, idx) => !rotationData.sideAtoms.includes(idx));
    expect(otherIdx).toBeGreaterThanOrEqual(0);
    const beforeOther = rotationData.positions[otherIdx];
    const afterOther = finalState[otherIdx];
    expect(afterOther.x).toBeCloseTo(beforeOther.x + 1, 6);
    expect(afterOther.y).toBeCloseTo(beforeOther.y + 1, 6);
    expect(afterOther.z).toBeCloseTo(beforeOther.z + 1, 6);
  });

  test('roy N-C bond rotation holds position after idle period and includes hydrogens', async ({ page, loadViewerPage, ensureWsReady }) => {
    await loadViewerPage({
      server: BACKEND_URL,
      query: { mol: 'molecules/roy.xyz', autoMD: 0 },
      disableAutoMd: true,
      testMode: true,
    });
    expect(await ensureWsReady()).toBeTruthy();

    const rotationData = await page.evaluate(() => {
      const api = window.viewerApi;
      const { state, manipulation } = api;
      const elements = state.elements || [];
      const bonds = (state.bonds || []).map((b, index) => ({
        i: b.i,
        j: b.j,
        key: b.key ?? `${b.i}-${b.j}`,
        index,
      }));
      const candidates = bonds.filter((b) => {
        const ei = elements[b.i];
        const ej = elements[b.j];
        return (
          (ei === 'N' && ej === 'C') ||
          (ei === 'C' && ej === 'N')
        );
      });
      if (!candidates.length) throw new Error('no N-C bonds found in roy.xyz');

      let selected = null;
      for (const bond of candidates) {
        for (const orientation of [0, 1]) {
          api.debugSelectBond({
            i: bond.i,
            j: bond.j,
            key: bond.key,
            index: bond.index,
            orientation,
          });
          const sel = api.debugGetSelection();
          const group = sel?.data?.rotationGroup;
          if (!group || !Array.isArray(group.sideAtoms) || !group.sideAtoms.length) continue;
          const hasHydrogen = group.sideAtoms.some((idx) => elements[idx] === 'H');
          if (hasHydrogen) {
            selected = { bond, orientation, group };
            break;
          }
        }
        if (selected) break;
      }
      if (!selected) throw new Error('no N-C bond with hydrogen rotation group');

      const sideAtoms = [...selected.group.sideAtoms];
      const allBefore = state.positions.map((p) => [p.x, p.y, p.z]);
      const before = sideAtoms.map((idx) => {
        const p = state.positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z, element: elements[idx] };
      });

      if (!manipulation.rotateBond(Math.PI / 6)) {
        throw new Error('rotateBond returned false');
      }

      const after = sideAtoms.map((idx) => {
        const p = state.positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z, element: elements[idx] };
      });

      return {
        sideAtoms,
        before,
        after,
        anchor: selected.group.anchor,
        movingRoot: selected.group.movingRoot,
        allBefore,
      };
    });

    expect(rotationData.sideAtoms.length).toBeGreaterThan(0);
    const hydrogenIndices = rotationData.after.filter((pt) => pt.element === 'H');
    expect(hydrogenIndices.length).toBeGreaterThan(0);

    const displacement = (a, b) => Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);

    for (const ptAfter of rotationData.after) {
      const before = rotationData.before.find((p) => p.idx === ptAfter.idx);
      expect(before).toBeTruthy();
      if (ptAfter.idx === rotationData.anchor || ptAfter.idx === rotationData.movingRoot) {
        continue;
      }
      const moved = displacement(ptAfter, before);
      expect(moved).toBeGreaterThan(0.02);
    }

    await page.waitForTimeout(2200);

    await page.evaluate(({ positions }) => {
      const ws = window.__fairchem_ws__;
      if (!ws || typeof ws.injectTestResult !== 'function') {
        throw new Error('ws injectTestResult unavailable');
      }
      ws.injectTestResult({
        seq: 999,
        client_seq: 0,
        userInteractionCount: 0,
        simStep: 0,
        positions,
      });
    }, { positions: rotationData.allBefore });

    await page.waitForTimeout(100);

    const afterFrame = await page.evaluate(({ indices }) => {
      const { positions, elements } = window.viewerApi.state;
      return indices.map((idx) => {
        const p = positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z, element: elements[idx] };
      });
    }, { indices: rotationData.sideAtoms });

    for (const ptLater of afterFrame) {
      const after = rotationData.after.find((p) => p.idx === ptLater.idx);
      expect(after).toBeTruthy();
      const drift = displacement(ptLater, after);
      expect(drift).toBeLessThan(0.01);
    }
  });

  test('scroll wheel rotates selected bond', async ({ page, loadViewerPage, ensureWsReady }) => {
    await loadViewerPage({
      server: BACKEND_URL,
      query: { mol: 'molecules/roy.xyz', autoMD: 0 },
      disableAutoMd: true,
      testMode: true,
    });
    expect(await ensureWsReady()).toBeTruthy();

    const setup = await page.evaluate(() => {
      const api = window.viewerApi;
      const { state } = api;
      const elements = state.elements || [];
      const bonds = (state.bonds || []).map((b, index) => ({
        i: b.i,
        j: b.j,
        key: b.key ?? `${b.i}-${b.j}`,
        index,
      }));
      const candidates = bonds.filter((b) => {
        const ei = elements[b.i];
        const ej = elements[b.j];
        return (ei === 'N' && ej === 'C') || (ei === 'C' && ej === 'N');
      });
      if (!candidates.length) throw new Error('no N-C bonds found');
      let selected = null;
      for (const bond of candidates) {
        for (const orientation of [0, 1]) {
          api.debugSelectBond({
            i: bond.i,
            j: bond.j,
            key: bond.key,
            index: bond.index,
            orientation,
          });
          const sel = api.debugGetSelection();
          const group = sel?.data?.rotationGroup;
          if (!group || !Array.isArray(group.sideAtoms) || !group.sideAtoms.length) continue;
          const hasHydrogen = group.sideAtoms.some((idx) => elements[idx] === 'H');
          if (hasHydrogen) {
            selected = {
              bond,
              orientation,
              sideAtoms: [...group.sideAtoms],
              anchor: group.anchor,
              movingRoot: group.movingRoot,
            };
            break;
          }
        }
        if (selected) break;
      }
      if (!selected) throw new Error('no N-C bond with hydrogen rotation group');
      const positions = selected.sideAtoms.map((idx) => {
        const p = state.positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z };
      });
      return { ...selected, positions };
    });

    expect(setup.sideAtoms.length).toBeGreaterThan(0);

    await page.$eval('canvas', (canvas) => {
      const evt = new WheelEvent('wheel', { deltaY: -120, bubbles: true, cancelable: true });
      canvas.dispatchEvent(evt);
    });
    await page.waitForTimeout(80);

    const afterScrollUp = await page.evaluate(({ indices }) => {
      return indices.map((idx) => {
        const p = window.viewerApi.state.positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z };
      });
    }, { indices: setup.sideAtoms });

    let maxDelta = 0;
    for (const pt of afterScrollUp) {
      const before = setup.positions.find((p) => p.idx === pt.idx);
      expect(before).toBeTruthy();
      const delta = Math.hypot(pt.x - before.x, pt.y - before.y, pt.z - before.z);
      if (delta > maxDelta) maxDelta = delta;
    }
    expect(maxDelta).toBeGreaterThan(0.01);

    await page.$eval('canvas', (canvas) => {
      const evt = new WheelEvent('wheel', { deltaY: 120, bubbles: true, cancelable: true });
      canvas.dispatchEvent(evt);
    });
    await page.waitForTimeout(80);

    const afterScrollDown = await page.evaluate(({ indices }) => {
      return indices.map((idx) => {
        const p = window.viewerApi.state.positions[idx];
        return { idx, x: p.x, y: p.y, z: p.z };
      });
    }, { indices: setup.sideAtoms });

    let deltaBetween = 0;
    for (const pt of afterScrollDown) {
      const prev = afterScrollUp.find((p) => p.idx === pt.idx);
      expect(prev).toBeTruthy();
      const delta = Math.hypot(pt.x - prev.x, pt.y - prev.y, pt.z - prev.z);
      if (delta > deltaBetween) deltaBetween = delta;
    }
    expect(deltaBetween).toBeGreaterThan(0.005);
  });
});
