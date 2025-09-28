import { createEventBus } from './eventBus.js';

export function createMoleculeState({ elements = [], positions = [], bonds = [], cell = null } = {}) {
  const bus = createEventBus();
  const state = {
    bus,
    elements: elements.slice(),
    positions: positions.map(p => ({ x:p.x, y:p.y, z:p.z })),
    bonds: bonds.map(b => ({ i:b.i, j:b.j })),
    cell: cell || { 
      a:{x:1,y:0,z:0},
      b:{x:0,y:1,z:0},
      c:{x:0,y:0,z:1},
      enabled:false,
      originOffset:{x:0,y:0,z:0}
    },
    // Ghost atoms: mapping from base atom index to array of ghost image objects { shift:[nx,ny,nz], pos:{x,y,z} }
    // For now just a placeholder container to avoid future shape changes in state.
    ghostImages: [],
    selection: { kind:null, data:null },
    dynamics: { velocities: [], forces: [], mass: [], temperature: 0 },
    versions: { positions:0, bonds:0, cell:0, selection:0, topology:0, dynamics:0 },
    markPositionsChanged() { this.versions.positions++; bus.emit('positionsChanged', state); },
    markBondsChanged() { this.versions.bonds++; bus.emit('bondsChanged', state); },
    markCellChanged() { this.versions.cell++; bus.emit('cellChanged', state); },
    markSelectionChanged() { this.versions.selection++; bus.emit('selectionChanged', state.selection); },
    markDynamicsChanged() { this.versions.dynamics++; bus.emit('dynamicsChanged', state); }
  };
  return state;
}
