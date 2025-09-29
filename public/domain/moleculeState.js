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
  showCell: false,
  showGhostCells: false,
    selection: { kind:null, data:null },
    dynamics: { velocities: [], forces: [], mass: [], temperature: 0 },
    versions: { positions:0, bonds:0, cell:0, selection:0, topology:0, dynamics:0 },
    markPositionsChanged() { this.versions.positions++; bus.emit('positionsChanged', state); },
    markBondsChanged() { this.versions.bonds++; bus.emit('bondsChanged', state); },
  markCellChanged() { this.versions.cell++; bus.emit('cellChanged', state); },
  toggleCellVisibility() { this.showCell = !this.showCell; this.markCellChanged(); },
  toggleGhostCells() { this.showGhostCells = !this.showGhostCells; this.markCellChanged(); },
  // Enhanced toggle: if user turns cell on but underlying cell is disabled (no lattice metadata loaded),
  // we synthesize a bounding box cell from current atom positions so something visible appears.
  // This avoids the confusing "Cell ON" label with no visual change.
  toggleCellVisibilityEnhanced() {
    this.showCell = !this.showCell;
    if (this.showCell && (!this.cell || !this.cell.enabled)) {
      if (!this.cell) this.cell = { a:{x:1,y:0,z:0}, b:{x:0,y:1,z:0}, c:{x:0,y:0,z:1}, enabled:false, originOffset:{x:0,y:0,z:0} };
      // Compute bounding box of current positions
      if (this.positions.length) {
        let minX=Infinity,minY=Infinity,minZ=Infinity,maxX=-Infinity,maxY=-Infinity,maxZ=-Infinity;
        for (const p of this.positions) { if (p.x<minX) minX=p.x; if (p.y<minY) minY=p.y; if (p.z<minZ) minZ=p.z; if (p.x>maxX) maxX=p.x; if (p.y>maxY) maxY=p.y; if (p.z>maxZ) maxZ=p.z; }
        // Ensure non‑zero extents with a small minimum so the box is visible for single atoms
        const eps = 1; // fallback scale when dimension collapses
        const dx = (maxX-minX)||eps;
        const dy = (maxY-minY)||eps;
        const dz = (maxZ-minZ)||eps;
        this.cell.a = { x:dx, y:0,  z:0 };
        this.cell.b = { x:0,  y:dy, z:0 };
        this.cell.c = { x:0,  y:0,  z:dz };
        this.cell.originOffset = { x:minX, y:minY, z:minZ };
  this.cell.enabled = true;
  this.cell.synthetic = true; // mark as synthetic so bond service can avoid periodic expansion
      } else {
        // No atoms yet; fall back to unit cube at origin
        this.cell.a = {x:1,y:0,z:0};
        this.cell.b = {x:0,y:1,z:0};
        this.cell.c = {x:0,y:0,z:1};
        this.cell.originOffset = {x:0,y:0,z:0};
  this.cell.enabled = true;
  this.cell.synthetic = true;
      }
    }
    this.markCellChanged();
  },
    markSelectionChanged() { this.versions.selection++; bus.emit('selectionChanged', state.selection); },
    markDynamicsChanged() { this.versions.dynamics++; bus.emit('dynamicsChanged', state); }
  };
  return state;
}
