/** @jest-environment jsdom */

// Verify: while dragging via VR path (shim), MD response does not move the held atom.

describe('VR drag excludes MD updates for held atom', () => {
  beforeEach(() => {
    document.body.innerHTML = '<canvas id="viewer" width="300" height="200"></canvas>';
    // Mock fetch responses for baseline and md
    global.fetch = jest.fn(async (url, opts) => {
      const u = String(url || '');
      if (u.includes('/serve/simple')) {
        return new Response(
          JSON.stringify({
            results: {
              energy: -2.34,
              forces: [
                [0, 0, 0],
                [0, 0, 0],
              ],
            },
          }),
          { status: 200 }
        );
      }
      if (u.includes('/serve/md')) {
        const body = opts && opts.body ? JSON.parse(opts.body) : {};
        const N = body && body.atomic_numbers ? body.atomic_numbers.length : 2;
        const positions = Array.from({ length: N }, (_, i) => [i * 1.0 + 10, 0, 0]);
        const forces = Array.from({ length: N }, () => [0, 0, 0]);
        const velocities = Array.from({ length: N }, () => [0, 0, 0]);
        return new Response(
          JSON.stringify({ positions, forces, final_energy: -1.0, velocities, temperature: 300 }),
          { status: 200 }
        );
      }
      return new Response(JSON.stringify({ ok: true }), { status: 200 });
    });
  });

  test('held index remains at dragged position after mdStep (VR shim)', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    const elements = ['H', 'H'];
    const positions = [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
    ];
    const bonds = [];
    const api = await initNewViewer(canvas, { elements, positions, bonds });
    // VR shim is created automatically in jsdom by createVRSupport.init()
    // We simulate VR-driven drag by calling through the same manipulation exposed to VR.
    api.selection.clickAtom(0);
    // Emulate VR intersector returning local position ~3,0,0
    const intersector = () => ({ x: 3.0, y: 0, z: 0 });
    // Begin drag (VR path uses same manipulation API injected into VR layer)
    const began = api.manipulation.beginDrag(intersector, {
      planePoint: { x: 0, y: 0, z: 0 },
      planeNormal: { x: 0, y: 1, z: 0 },
    });
    expect(began).toBe(true);
    const moved = api.manipulation.updateDrag(intersector);
    expect(moved).toBe(true);
    const held = { ...api.state.positions[0] };
    // Trigger an MD step which would shift all atoms by +10 on x
    const res = await api.mdStep({ temperature: 298 });
    expect(res && res.applied).toBe(true);
    // Held atom unchanged; other atom updated
    const p0 = api.state.positions[0];
    const p1 = api.state.positions[1];
    expect(Math.abs(p0.x - held.x)).toBeLessThan(1e-6);
    expect(Math.abs(p1.x - 11.0)).toBeLessThan(1e-6);
    api.manipulation.endDrag();
  });
});
