import { ensureBabylonStub } from './utils/babylonStub.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function createBus() {
  const listeners = Object.create(null);
  return {
    on(event, handler) {
      (listeners[event] = listeners[event] || []).push(handler);
    },
    emit(event) {
      const items = listeners[event];
      if (!items) return;
      for (const handler of items) {
        try {
          handler();
        } catch {}
      }
    },
  };
}

function createScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    onPointerObservable: { add() {} },
  };
}

describe('moleculeView mesh modes', () => {
  beforeEach(() => {
    ensureBabylonStub();
    global.window = {
      __MLIPVIEW_TEST_MODE: true,
      location: { search: '' },
    };
  });

  afterEach(() => {
    delete global.window;
  });

  test('bonds with low opacity spawn on the translucent master by default', () => {
    const state = {
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1, opacity: 0.2 }],
      bus: createBus(),
      showGhostCells: false,
      showCell: false,
      cell: { enabled: false },
    };
    state.markPositionsChanged = () => state.bus.emit('positionsChanged');
    state.markCellChanged = () => state.bus.emit('cellChanged');

    const scene = createScene();
    const view = createMoleculeView(scene, state);

    const modes = view.getBondMeshModes();
    expect(modes.default).toEqual(['soft']);
    expect(modes.current).toEqual(['soft']);

    const bondGroup = Array.from(view._internals.bondGroups.values())[0];
    expect(bondGroup.instances.soft.length).toBe(1);
    expect(bondGroup.instances.solid.length).toBe(0);
  });

  test('applyMeshModes migrates atoms and bonds between solid and soft masters', () => {
    const state = {
      elements: ['C', 'O'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.2, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1, opacity: 1 }],
      bus: createBus(),
      showGhostCells: false,
      showCell: false,
      cell: { enabled: false },
    };
    state.markPositionsChanged = () => state.bus.emit('positionsChanged');
    state.markCellChanged = () => state.bus.emit('cellChanged');

    const scene = createScene();
    const view = createMoleculeView(scene, state);

    const baselineAtoms = view.getAtomMeshModes();
    const baselineBonds = view.getBondMeshModes();
    expect(baselineAtoms.current.every((mode) => mode === 'solid')).toBe(true);
    expect(baselineBonds.current.every((mode) => mode === 'solid')).toBe(true);

    view.applyMeshModes({
      atomDefault: ['solid', 'solid'],
      atomCurrent: ['soft', 'solid'],
      bondDefault: ['solid'],
      bondCurrent: ['soft'],
    });

    const currentAtoms = view.getAtomMeshModes();
    const currentBonds = view.getBondMeshModes();
    expect(currentAtoms.current).toEqual(['soft', 'solid']);
    expect(currentBonds.current).toEqual(['soft']);

    const atomGroups = Array.from(view._internals.atomGroups.values());
    const carbonGroup = atomGroups.find((g) => g.key === 'C');
    const oxygenGroup = atomGroups.find((g) => g.key === 'O');
    expect(carbonGroup.instances.soft.length).toBe(1);
    expect(carbonGroup.instances.solid.length).toBe(0);
    expect(oxygenGroup.instances.soft.length).toBe(0);
    expect(oxygenGroup.instances.solid.length).toBe(1);

    const bondGroup = Array.from(view._internals.bondGroups.values())[0];
    expect(bondGroup.instances.soft.length).toBe(1);
    expect(bondGroup.instances.solid.length).toBe(0);
  });
});
