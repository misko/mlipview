import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function vec(x, y, z) {
  return { x, y, z };
}

function makeCell(enabled) {
  return {
    a: { x: 10, y: 0, z: 0 },
    b: { x: 0, y: 10, z: 0 },
    c: { x: 0, y: 0, z: 10 },
    originOffset: { x: 0, y: 0, z: 0 },
    enabled,
  };
}

function ghostInstanceCount(view) {
  return Array.from(view._internals.ghostBondGroups.values()).reduce(
    (sum, grp) => sum + (grp.mats?.length || 0),
    0
  );
}

describe('bond service periodic ghost handling', () => {
  function recompute(state) {
    const service = createBondService(state);
    service.recomputeAndStore();
    return state;
  }

  test('close atoms gain six ghost bonds when periodic is enabled', () => {
    const baseState = createMoleculeState({
      elements: ['C', 'C'],
      positions: [vec(0, 0, 0), vec(1.4, 0, 0)],
      cell: makeCell(false),
    });

    recompute(baseState);
    expect(baseState.bonds).toHaveLength(1);
    expect(baseState.ghostBondMeta).toHaveLength(0);

    baseState.cell.enabled = true;
    recompute(baseState);
    expect(baseState.bonds).toHaveLength(1);
    expect(baseState.ghostBondMeta).toHaveLength(6);
    const shifts = baseState.ghostBondMeta
      .map((meta) => `${meta.shiftA.join(',')}|${meta.shiftB.join(',')}`)
      .sort();
    expect(shifts).toEqual([
      '-1,0,0|-1,0,0',
      '0,-1,0|0,-1,0',
      '0,0,-1|0,0,-1',
      '0,0,1|0,0,1',
      '0,1,0|0,1,0',
      '1,0,0|1,0,0',
    ]);
  });

  test('distant atoms only render ghost bonds across the boundary (near ±half-box)', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [vec(-4.5, 0, 0), vec(4.5, 0, 0)],
      cell: {
        a: { x: 10, y: 0, z: 0 },
        b: { x: 0, y: 10, z: 0 },
        c: { x: 0, y: 0, z: 10 },
        originOffset: { x: -5, y: -5, z: -5 },
        enabled: false,
      },
    });

    recompute(state);
    expect(state.bonds).toHaveLength(0);
    expect(state.ghostBondMeta).toHaveLength(0);

    state.cell.enabled = true;
    recompute(state);
    expect(state.bonds).toHaveLength(0);
    expect(state.ghostBondMeta).toHaveLength(2);
    const shifts = state.ghostBondMeta
      .map((meta) => `${meta.shiftA.join(',')}|${meta.shiftB.join(',')}`)
      .sort();
    expect(shifts).toEqual(['-1,0,0|-1,0,0', '1,0,0|1,0,0']);
  });

  test('atoms near cell edges produce two ghost bonds when periodic is enabled', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [vec(9.5, 0, 0), vec(-9.5, 0, 0)],
      cell: {
        a: { x: 10, y: 0, z: 0 },
        b: { x: 0, y: 10, z: 0 },
        c: { x: 0, y: 0, z: 10 },
        originOffset: { x: 0, y: 0, z: 0 },
        enabled: false,
      },
    });

    recompute(state);
    expect(state.bonds).toHaveLength(0);
    expect(state.ghostBondMeta).toHaveLength(0);

    state.cell.enabled = true;
    recompute(state);
    expect(state.bonds).toHaveLength(0);
    expect(state.ghostBondMeta).toHaveLength(2);
    const shifts = state.ghostBondMeta
      .map((meta) => `${meta.shiftA.join(',')}|${meta.shiftB.join(',')}`)
      .sort();
    expect(shifts).toEqual(['-1,0,0|-1,0,0', '1,0,0|1,0,0']);
  });
});

describe('molecule view ghost refresh', () => {
  test('rebuildBonds propagates ghost instances without manual rebuild', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [vec(0, 0, 0), vec(1.4, 0, 0)],
      cell: makeCell(true),
    });
    state.showCell = true;
    state.showGhostCells = true;
    const view = createMoleculeView({ onPointerObservable: { add() {} } }, state);
    expect(ghostInstanceCount(view)).toBe(0);

    const service = createBondService(state);
    service.recomputeAndStore();

    expect(state.ghostBondMeta).toHaveLength(6);
    expect(ghostInstanceCount(view)).toBe(6);
  });
});
