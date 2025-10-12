/** @jest-environment jsdom */

// Verify: while dragging an atom (desktop/manipulation path), MD response does not move that atom.

describe('desktop drag excludes MD updates for held atom', () => {
  beforeEach(() => {
    // Basic DOM canvas
    document.body.innerHTML = '<canvas id="viewer" width="300" height="200"></canvas>';
    // Mock fetch for baseline simple + md endpoints
    global.fetch = jest.fn(async (url, opts) => {
      const u = String(url||'');
      // simple_calculate baseline
      if (u.includes('/serve/simple')) {
        return new Response(JSON.stringify({ results: { energy: -1.23, forces: [[0,0,0],[0,0,0]] } }), { status: 200 });
      }
      // md step: try to move both atoms by +10, +0, +0
      if (u.includes('/serve/md')) {
        const body = opts && opts.body ? JSON.parse(opts.body) : {};
        const N = (body && body.atomic_numbers) ? body.atomic_numbers.length : 2;
        const positions = Array.from({ length: N }, (_, i) => [i*1.0 + 10, 0, 0]);
        const forces = Array.from({ length: N }, () => [0,0,0]);
        const velocities = Array.from({ length: N }, () => [0,0,0]);
        return new Response(JSON.stringify({ positions, forces, final_energy: -1.0, velocities, temperature: 300 }), { status: 200 });
      }
      // relax not used here
      return new Response(JSON.stringify({ ok: true }), { status: 200 });
    });
  });

  test('held index remains at dragged position after mdStep', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const elements = ['H','H'];
    const positions = [{ x:0,y:0,z:0 }, { x:1,y:0,z:0 }];
    const bonds = [];
    const api = await initNewViewer(canvas, { elements, positions, bonds });
    // Select atom 0
    api.selection.clickAtom(0);
    // Begin drag with a trivial intersector returning a constant point near 2,0,0
    const intersector = () => ({ x: 2.0, y: 0, z: 0 });
    const started = api.manipulation.beginDrag(intersector, { planePoint:{x:0,y:0,z:0}, planeNormal:{x:0,y:1,z:0} });
    expect(started).toBe(true);
    // Update drag once -> move atom 0 to approx 2,0,0
    const moved = api.manipulation.updateDrag(intersector);
    expect(moved).toBe(true);
    const held = { ...api.state.positions[0] };
    // Now perform an MD step which attempts to move all atoms by +10 on x
    const res = await api.mdStep({ temperature: 298 });
    expect(res && res.applied).toBe(true);
    // Atom 0 should remain at held position; atom 1 should update to ~11
    const p0 = api.state.positions[0];
    const p1 = api.state.positions[1];
    expect(Math.abs(p0.x - held.x)).toBeLessThan(1e-6);
    expect(Math.abs(p1.x - 11.0)).toBeLessThan(1e-6);
    // End drag cleanup
    api.manipulation.endDrag();
  });
});
