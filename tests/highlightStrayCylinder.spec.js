// Regression test: ensure highlight bond cylinder remains hidden when no bond is selected.
// Uses global BABYLON stub from tests/jest.setup.js (no extra mocking here).

// Minimal molecule state stub with event bus
function createBus(){const map=new Map();return {on:(e,f)=>{(map.get(e)||map.set(e,[]).get(e)).push(f);}, emit:(e)=>{(map.get(e)||[]).forEach(f=>f());}};}

function createMolState(){
  return { elements:['C','C'], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], bonds:[{i:0,j:1,opacity:1}], selection:null, showCell:false, showGhostCells:false, cell:{enabled:false}, bus:createBus() };
}

// Import the view factory from the new viewer path (no path/fs needed)

describe('highlight stray bond cylinder regression', () => {
  test('highlight bond mesh disabled when no bond selected', async () => {
    const molState = createMolState();
    // dynamic import using ESM-compatible path
  const mod = await import('../public/render/moleculeView.js');
    const { createMoleculeView } = mod;
    const scene = {};
    const view = createMoleculeView(scene, molState);
    const highlight = view._internals.highlight;
    expect(highlight.bond.isVisible).toBe(false);
    molState.selection = { kind:'atom', data:{ index:0 } };
    molState.bus.emit('selectionChanged');
    expect(highlight.bond.isVisible).toBe(false);
    molState.selection = { kind:'bond', data:{ i:0, j:1 } };
    molState.bus.emit('selectionChanged');
    expect(highlight.bond.isVisible).toBe(true);
    molState.selection = null;
    molState.bus.emit('selectionChanged');
    expect(highlight.bond.isVisible).toBe(false);
  });
});
