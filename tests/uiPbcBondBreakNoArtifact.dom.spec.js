/** @jest-environment jsdom */

import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Minimal manipulation stub to stretch bond like user dragging would do
function createManipulation(state){
  return {
    rotateBond({ i, j, theta }){
      // Just move j along +x proportional to theta to simulate stretch until bond breaks
      const p = state.positions[j];
      state.positions[j] = { x: p.x + Math.abs(theta)*2.5, y: p.y, z: p.z };
      state.markPositionsChanged();
      return true;
    }
  };
}

describe('UI: PBC enabled, breaking bond does not leave origin artifact', () => {
  test('no stray ghost cylinder/sphere at origin after bond disappears', () => {
    // Arrange: C-C in a PBC box
    const cell = { a:{x:10,y:0,z:0}, b:{x:0,y:10,z:0}, c:{x:0,y:0,z:10}, originOffset:{x:0,y:0,z:0}, enabled:true };
    const state = createMoleculeState({ elements:['C','C'], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], cell });
    state.showCell = true; state.showGhostCells = true;
    const bonds = createBondService(state); bonds.recomputeAndStore();
    const scene = { onPointerObservable:{ add(){} } };
    const view = createMoleculeView(scene, state);
    view.rebuildBonds(); view.rebuildGhosts();

    const manipulation = createManipulation(state);
    // Act: rotateBond with a large theta to break the bond (simulated stretch)
    manipulation.rotateBond({ i:0, j:1, theta: 2.0 });
    bonds.recomputeAndStore();
    view.rebuildBonds(); view.rebuildGhosts();

    // Assert: no ghost instances and masters not visible; avoids a cylinder/sphere at origin
    let ghostInst = 0; let anyVisible=false;
    for (const g of view._internals.ghostBondGroups.values()) { ghostInst += g.mats.length; anyVisible = anyVisible || g.master.isVisible; }
    for (const g of view._internals.ghostAtomGroups.values()) { /* atom ghosts may remain */ if (g.mats.length===0) anyVisible = anyVisible || g.master.isVisible; }
    expect(ghostInst).toBe(0);
    expect(anyVisible).toBe(false);
  });
});
