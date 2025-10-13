import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function dist(a,b){ const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return Math.hypot(dx,dy,dz); }

describe('PBC along-X primary long bond suppression (10x10x10, C @ 0.1 and 9.9 on X)', () => {
  test('long primary bonds across X are invisible (opacity 0) and ghost shows short bond', () => {
    const cell = { a:{x:10,y:0,z:0}, b:{x:0,y:10,z:0}, c:{x:0,y:0,z:10}, originOffset:{x:0,y:0,z:0}, enabled:true };
    const state = createMoleculeState({
      elements:['C','C'],
      positions:[ {x:0.1,y:0.1,z:0.1}, {x:9.9,y:0.1,z:0.1} ],
      bonds:[],
      cell
    });
    state.showCell = true;
    state.showGhostCells = true;
    const bondService = createBondService(state);
    bondService.recomputeAndStore();

    const scene = { onPointerObservable:{ add(){} } };
    const view = createMoleculeView(scene, state);
    view.rebuildBonds();
    view.rebuildGhosts();

    // Primary bonds: if any span long distance (>5 Å here), they must be fully invisible (opacity 0)
    const longThreshold = 5.0;
    let sawLongPrimary = 0; let anyVisibleLong = 0;
    for (const [key, g] of view._internals.bondGroups.entries()) {
      for (const b of g.indices) {
        const pA = state.positions[b.i];
        const pB = state.positions[b.j];
        const L = dist(pA,pB);
        if (L > longThreshold) {
          sawLongPrimary++;
          const alpha = (b.opacity != null) ? b.opacity : 1.0;
          if (alpha > 1e-3) anyVisibleLong++;
        }
      }
    }
    // Expect at least one long primary instance (the crossing), but it must be invisible
    expect(sawLongPrimary).toBeGreaterThanOrEqual(1);
    expect(anyVisibleLong).toBe(0);

    // Ghost bonds: at least one should exist bridging the short image across X
    const ghostCount = Array.from(view._internals.ghostBondGroups.values()).reduce((s,g)=> s + g.mats.length, 0);
    expect(ghostCount).toBeGreaterThan(0);
  });
});
