import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function countGhostBondInstances(view){
  return Array.from(view._internals.ghostBondGroups.values()).reduce((s,g)=> s + (g.mats?.length||0), 0);
}

describe('PBC ghost bond masters disable when bonds disappear', () => {
  test('breaking a C-C bond clears ghost instances and disables masters (no origin artifact)', () => {
    // 10x10x10 box so neighbor images are far; show PBC + ghosts
    const cell = { a:{x:10,y:0,z:0}, b:{x:0,y:10,z:0}, c:{x:0,y:0,z:10}, originOffset:{x:0,y:0,z:0}, enabled:true };
    const state = createMoleculeState({
      elements:['C','C'],
      positions:[ {x:0,y:0,z:0}, {x:1.42,y:0,z:0} ], // bonded distance
      bonds:[],
      cell
    });
    state.showCell = true; state.showGhostCells = true;
    const bondService = createBondService(state);
    bondService.recomputeAndStore();
    const scene = { onPointerObservable:{ add(){} } };
    const view = createMoleculeView(scene, state);
    view.rebuildBonds();
    view.rebuildGhosts();

    // Expect some ghost bonds to be instantiated for same-shift images
    const initialGhosts = countGhostBondInstances(view);
    expect(initialGhosts).toBeGreaterThan(0);

    // Stretch atoms beyond bonding cutoff to break the bond
    state.positions[1] = { x: 6.0, y: 0, z: 0 };
    state.markPositionsChanged();
    bondService.recomputeAndStore();
    view.rebuildBonds();
    view.rebuildGhosts();

    // All ghost bond instance buffers should be empty and masters disabled (not visible)
    let totalGhost = 0; let anyMasterVisible = false;
    for (const g of view._internals.ghostBondGroups.values()) {
      totalGhost += g.mats.length;
      const m = g.master; anyMasterVisible = anyMasterVisible || !!(m && m.isVisible);
    }
    expect(totalGhost).toBe(0);
    expect(anyMasterVisible).toBe(false);
  });
});
