import { createPickingService } from '../public/core/pickingService.js';

describe('x-picking-service', () => {
  test('pickAtPointer returns null when scene pick misses', () => {
    const scene = {
      pointerX: 0,
      pointerY: 0,
      pick: () => ({ hit: false }),
      onPointerObservable: { add() {} },
    };
    const view = { resolveAtomPick: () => null, resolveBondPick: () => null };
    const selection = { clickAtom: jest.fn(), clickBond: jest.fn(), clear: jest.fn() };
    const picking = createPickingService(scene, view, selection);
    expect(picking.pickAtPointer()).toBeNull();
  });
});
