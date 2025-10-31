import { createEventBus } from './eventBus.js';
import { computeOrthoCellFromPositions, wrapPositionsInPlace } from '../util/pbc.js';
import { __count } from '../util/funcCount.js';

export function createMoleculeState({
  elements = [],
  positions = [],
  bonds = [],
  cell = null,
} = {}) {
  __count('moleculeState#createMoleculeState');
  const bus = createEventBus();
  const state = {
    bus,
    elements: elements.slice(),
    positions: positions.map((p) => ({ x: p.x, y: p.y, z: p.z })),
    bonds: bonds.map((b) => ({
      i: b.i,
      j: b.j,
      opacity: typeof b.opacity === 'number' ? b.opacity : 1,
      crossing: !!b.crossing,
      imageDelta: Array.isArray(b.imageDelta) ? b.imageDelta.slice(0, 3) : [0, 0, 0],
      cellOffsetA: Array.isArray(b.cellOffsetA) ? b.cellOffsetA.slice(0, 3) : [0, 0, 0],
      cellOffsetB: Array.isArray(b.cellOffsetB) ? b.cellOffsetB.slice(0, 3) : [0, 0, 0],
    })),
    cell: cell || {
      a: { x: 1, y: 0, z: 0 },
      b: { x: 0, y: 1, z: 0 },
      c: { x: 0, y: 0, z: 1 },
      enabled: false,
      originOffset: { x: 0, y: 0, z: 0 },
    },
    // Ghost atoms: mapping from base atom index to array of ghost image objects { shift:[nx,ny,nz], pos:{x,y,z} }
    // For now just a placeholder container to avoid future shape changes in state.
    ghostImages: [],
    ghostBondMeta: [],
    showCell: false,
    showGhostCells: false,
    // Forces visualization toggle: OFF by default, but enable automatically in test mode for visualization tests
    showForces:
      typeof window !== 'undefined' && window.__MLIPVIEW_TEST_MODE === true ? true : false,
    selection: { kind: null, data: null },
    dynamics: { velocities: [], forces: [], mass: [], temperature: 0 },
    versions: { positions: 0, bonds: 0, cell: 0, selection: 0, topology: 0, dynamics: 0 },
    markPositionsChanged() {
      __count('moleculeState#markPositionsChanged');
      // If PBC visual is on and cell is enabled, wrap atom coordinates into primary cell
      try {
        if (this.showCell && this.cell && this.cell.enabled) {
          wrapPositionsInPlace(this.positions, this.cell);
        }
      } catch {}
      this.versions.positions++;
      bus.emit('positionsChanged', state);
    },
    markBondsChanged() {
      __count('moleculeState#markBondsChanged');
      this.versions.bonds++;
      bus.emit('bondsChanged', state);
    },
    markCellChanged() {
      __count('moleculeState#markCellChanged');
      this.versions.cell++;
      bus.emit('cellChanged', state);
    },
    toggleCellVisibility() {
      __count('moleculeState#toggleCellVisibility');
      this.showCell = !this.showCell;
      this.markCellChanged();
    },
    toggleGhostCells() {
      __count('moleculeState#toggleGhostCells');
      this.showGhostCells = !this.showGhostCells;
      this.markCellChanged();
    },
    toggleForceVectorsVisibility() {
      __count('moleculeState#toggleForceVectorsVisibility');
      this.showForces = !this.showForces;
      bus.emit('forcesChanged', state);
    },
    // Enhanced toggle: if user turns cell on but underlying cell is disabled (no lattice metadata loaded),
    // we synthesize a bounding box cell from current atom positions so something visible appears.
    // This avoids the confusing "Cell ON" label with no visual change.
    toggleCellVisibilityEnhanced() {
      __count('moleculeState#toggleCellVisibilityEnhanced');
      this.showCell = !this.showCell;
      if (this.showCell && (!this.cell || !this.cell.enabled)) {
        // Ensure cell object exists
        if (!this.cell)
          this.cell = {
            a: { x: 1, y: 0, z: 0 },
            b: { x: 0, y: 1, z: 0 },
            c: { x: 0, y: 0, z: 1 },
            enabled: false,
            originOffset: { x: 0, y: 0, z: 0 },
          };
        // Compute orthorhombic cell that encloses all atoms with 1 Å padding on all sides
        const computed = computeOrthoCellFromPositions(this.positions, 1.0);
        this.cell.a = computed.a;
        this.cell.b = computed.b;
        this.cell.c = computed.c;
        this.cell.originOffset = computed.originOffset;
        this.cell.enabled = true;
        this.cell.synthetic = true; // mark as synthetic for downstream consumers
      }
      this.markCellChanged();
    },
    markSelectionChanged() {
      __count('moleculeState#markSelectionChanged');
      this.versions.selection++;
      bus.emit('selectionChanged', state.selection);
    },
    markDynamicsChanged() {
      __count('moleculeState#markDynamicsChanged');
      this.versions.dynamics++;
      bus.emit('dynamicsChanged', state);
    },
  };
  return state;
}
