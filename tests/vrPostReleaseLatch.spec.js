/** @jest-environment jsdom */

// This test reproduces a snap-back on joystick drag release: an MD response landing
// immediately after endDrag attempts to re-apply the pre-drag depth. We expect the
// just-dragged atom to be latched and excluded briefly, keeping its user-set position.

function mkResp(positions){ return { positions, forces: positions.map(()=>[0,0,0]), final_energy: -1.0, velocities: positions.map(()=>[0,0,0]), temperature: 300 }; }

const queue = [];

beforeEach(()=>{
  document.body.innerHTML = '<canvas id="viewer" width="300" height="200"></canvas>';
  jest.useFakeTimers();
  jest.setSystemTime(0);
  global.fetch = jest.fn(async (url, opts) => {
    const u = String(url||'');
    if (u.includes('/serve/simple')) {
      return new Response(JSON.stringify({ results: { energy: -2.0, forces: [[0,0,0],[0,0,0]] } }), { status: 200 });
    }
    if (u.includes('/serve/md')) {
      const resp = queue.length ? queue.shift() : mkResp([[0,0,0],[1,0,0]]);
      return new Response(JSON.stringify(resp), { status: 200 });
    }
    return new Response(JSON.stringify({ ok: true }), { status: 200 });
  });
});

afterEach(()=>{
  jest.useRealTimers();
});

function almostEqual(a,b,eps=1e-9){ return Math.abs(a-b) <= eps; }

describe('VR joystick drag release does not snap back', () => {
  test('MD response immediately after endDrag does not override just-dragged atom', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const elements = ['H','H'];
    const positions = [{ x:0,y:0,z:0 }, { x:1,y:0,z:0 }];
    const bonds = [];
    const api = await initNewViewer(canvas, { elements, positions, bonds });

    // Select atom 0 and perform a VR-style drag using manipulation API
    api.selection.clickAtom(0);
    const intersector = () => ({ x: 2.5, y: 0, z: 0 });
    expect(api.manipulation.beginDrag(intersector)).toBe(true);
    expect(api.manipulation.updateDrag(intersector)).toBe(true);
    const dragged = { ...api.state.positions[0] };

    // End drag, then immediately enqueue an MD response that would snap atom 0 back to x=0
    api.manipulation.endDrag();
    const snapBack = mkResp([[0,0,0],[1,0,0]]);
    queue.push(snapBack);
    const res = await api.mdStep({});
    expect(res && res.applied).toBe(true);

    // Atom 0 should remain at dragged.x due to latch, atom 1 can be updated by MD
    const p0 = api.state.positions[0];
    const p1 = api.state.positions[1];
    expect(almostEqual(p0.x, dragged.x)).toBe(true);
    expect(almostEqual(p1.x, 1)).toBe(true);

    // Advance time beyond latch window (~0.8s) and send another MD update; now atom 0 should update
    jest.advanceTimersByTime(1000);
    const forward = mkResp([[10,0,0],[11,0,0]]);
    queue.push(forward);
    await api.mdStep({});
    const q0 = api.state.positions[0];
    const q1 = api.state.positions[1];
    expect(almostEqual(q0.x, 10)).toBe(true);
    expect(almostEqual(q1.x, 11)).toBe(true);
  });
});
