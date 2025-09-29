// Reproduce red (oversized) sphere artifact on molecule switch.
// Hypothesis: highlight atom sphere remains enabled after selection indices become invalid
// when switching from a larger molecule (ROY) to a smaller one (benzene), leaving it at origin
// or with stale scale.

// We'll simulate by: load ROY, select an atom, then switch to Benzene and assert highlight hidden.

import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

async function loadXYZText(name){
  const fs = await import('fs');
  const path = await import('path');
  const p = path.join(process.cwd(),'public','molecules',name);
  if (!fs.existsSync(p)) return null; // allow skip
  return fs.readFileSync(p,'utf8');
}

function buildStateFromXYZ(xyzText){
  const parsed = parseXYZ(xyzText);
  const state = { elements:[], positions:[], bonds:[], selection:null, showCell:false, showGhostCells:false, cell:{enabled:false}, bus:createBus() };
  // Inject minimal event emitters used by xyzLoader adapter
  state.markPositionsChanged = ()=> state.bus.emit('positionsChanged');
  state.markBondsChanged = ()=> state.bus.emit('bondsChanged');
  state.markCellChanged = ()=> state.bus.emit('cellChanged');
  applyXYZToState(state, parsed);
  return { state, parsed };
}

function createBus(){const map=new Map();return {on:(e,f)=>{(map.get(e)||map.set(e,[]).get(e)).push(f);}, emit:(e)=>{(map.get(e)||[]).forEach(fn=>fn());}};}

describe('red sphere artifact regression', () => {
  test('atom highlight disabled after molecule switch', async () => {
    const roy = await loadXYZText('roy.xyz');
    const benz = await loadXYZText('benzene.xyz');
    if (!roy || !benz) {
      console.warn('[redSphereSwitch.spec] Skipped due to missing fixture files');
      return;
    }
    const { state: royState } = buildStateFromXYZ(roy);
    const view = createMoleculeView({}, royState);
    // Select an atom in ROY (e.g., index 0)
    royState.selection = { kind:'atom', data:{ index:0 } };
    royState.bus.emit('selectionChanged');
    const hi = view._internals.highlight;
    expect(hi.atom.isVisible).toBe(true);
    // Switch to Benzene: mutate state in-place simulating loader
    const parsedB = parseXYZ(benz);
    // Provide required mark* functions if lost
    royState.markPositionsChanged = ()=> royState.bus.emit('positionsChanged');
    royState.markBondsChanged = ()=> royState.bus.emit('bondsChanged');
    royState.markCellChanged = ()=> royState.bus.emit('cellChanged');
    applyXYZToState(royState, parsedB); // overrides & emits events
    // Simulate recompute (empty -> should hide existing highlight)
    royState.bus.emit('bondsChanged');
    // After switch highlight atom should be hidden (no lingering large sphere)
    expect(hi.atom.isVisible).toBe(false);
  });
});
