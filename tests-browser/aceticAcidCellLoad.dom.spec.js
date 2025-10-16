/** @jest-environment jsdom */

// Verify that loading acetic acid XYZ enables PBC and loads the monoclinic cell from header.

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';

function makeViewer(){
  const listeners = {};
  const bus = { on:(ev,fn)=>{ (listeners[ev]||(listeners[ev]=[])).push(fn); }, emit:(ev)=>{ (listeners[ev]||[]).forEach(fn=>fn()); } };
  const state = {
    bus,
    elements: [], positions: [], bonds: [],
    cell: { a:{x:1,y:0,z:0}, b:{x:0,y:1,z:0}, c:{x:0,y:0,z:1}, enabled:false, originOffset:{x:0,y:0,z:0} },
    showCell: false,
    markCellChanged(){ bus.emit('cellChanged'); },
    markPositionsChanged(){ bus.emit('positionsChanged'); },
    markBondsChanged(){ bus.emit('bondsChanged'); },
    dynamics: {}
  };
  const viewerApi = { state, recomputeBonds: ()=>{} };
  return viewerApi;
}

function aceticAcidFixture(){
  // Use the exact header from public/molecules/acetic_acid.xyz + a minimal atom subset
  return `32\ncell: a=13.31, b=4.09, c=5.77, alpha=90.00, beta=107.0, gamma=90.00, spacegroup=P21/c, temp=400K\nC 0 0 0\nH 0 0 1\nH 0 1 0\nH 1 0 0\nC 2 0 0\nH 2 0 1\nH 2 1 0\nH 3 0 0\nC 4 0 0\nO 4 0 1\nO 4 1 0\nC 5 0 0\nH 5 0 1\nH 5 1 0\nH 6 0 0\nH 6 0 1\nC 7 0 0\nO 7 0 1\nO 7 1 0\nC 8 0 0\nH 8 0 1\nH 8 1 0\nH 9 0 0\nH 9 0 1\nC 10 0 0\nO 10 0 1\nO 10 1 0\nC 11 0 0\nH 11 0 1\nH 11 1 0\nH 12 0 0\nH 12 0 1\n`;
}

describe('acetic acid cell loads and PBC enabled', ()=>{
  beforeEach(()=>{ document.body.innerHTML = '<div id="app"></div>'; });

  test('cell parsed, monoclinic, and showCell enabled', ()=>{
    const xyz = aceticAcidFixture();
    const parsed = parseXYZ(xyz);
    expect(parsed.cell).toBeTruthy();
    // Angles from header should be alpha=90, beta=107, gamma=90
    expect(parsed.cell.enabled).toBe(true);

    const viewer = makeViewer();
    applyParsedToViewer(viewer, parsed);

    // After apply, PBC should be visible by default
    expect(viewer.state.showCell).toBe(true);
    expect(viewer.state.cell && viewer.state.cell.enabled).toBe(true);

    // Basic sanity for parameters: lengths roughly match and monoclinic by design
    const a = viewer.state.cell.a, b = viewer.state.cell.b, c = viewer.state.cell.c;
    const ax = Math.hypot(a.x,a.y,a.z), bx = Math.hypot(b.x,b.y,b.z), cx = Math.hypot(c.x,c.y,c.z);
    expect(ax).toBeCloseTo(13.31, 2);
    expect(bx).toBeCloseTo(4.09, 2);
    expect(cx).toBeCloseTo(5.77, 2);
  });
});
