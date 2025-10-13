import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Mock browser-style test: two carbons near opposite faces of a 10x10x10 cell.
// Expectation: there should be NO long opaque primary bond spanning the box.
// If a long bond is rendered, we log detailed diagnostics and fail.

function dist(a,b){ const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return Math.hypot(dx,dy,dz); }

describe('PBC long primary bond crossing (10x10x10 cell, C @ 0.1 and 9.9)', () => {
  test('no long opaque primary bonds across the cell diagonal', () => {
    // Arrange: explicit cell and atom positions
    const cell = {
      a:{x:10,y:0,z:0}, b:{x:0,y:10,z:0}, c:{x:0,y:0,z:10},
      originOffset:{x:0,y:0,z:0}, enabled:true
    };
    const state = createMoleculeState({
      elements:['C','C'],
      positions:[ {x:0.1,y:0.1,z:0.1}, {x:9.9,y:9.9,z:9.9} ],
      bonds:[],
      cell
    });
    state.showCell = true;
    // Turn on ghost rendering just to mirror UI flow; primary vs ghost separation is what we inspect
    state.showGhostCells = true;

    // Use real bond computation respecting PBC/cell.
    const bondService = createBondService(state);
    bondService.recomputeAndStore();

  // Minimal scene stub for view
    const scene = { onPointerObservable:{ add(){} } };
    const view = createMoleculeView(scene, state);
    view.rebuildBonds();
    view.rebuildGhosts();

    // Inspect real bond instances (primary, non-ghost). A long primary bond across the box is unexpected.
    const longThreshold = 6.0; // > cell/2 along any axis would be unexpected; diagonal ~17.3 here
    const longs = [];
    for (const [key, g] of view._internals.bondGroups.entries()) {
      for (const b of g.indices) {
        const pA = state.positions[b.i];
        const pB = state.positions[b.j];
        const L = dist(pA,pB);
        if (L > longThreshold) {
          longs.push({ key, i:b.i, j:b.j, L, opacity: b.opacity ?? 1 });
          // Console diagnostics: dump atom positions and cell
          // Note: We print even if opacity is 0, to reveal invisible-but-present long instances.
          const c = state.cell;
          console.log('[LONG-BOND][primary] key=', key, 'i=', b.i, 'j=', b.j, 'len=', L.toFixed(3), 'opacity=', (b.opacity ?? 1));
          console.log('[LONG-BOND][atoms] A=', pA, 'B=', pB);
          console.log('[LONG-BOND][cell] a=', c.a, 'b=', c.b, 'c=', c.c, 'originOffset=', c.originOffset);
        }
      }
    }

    // Assert: there should be ZERO long primary bonds.
    expect(longs.length).toBe(0);
  });
});
