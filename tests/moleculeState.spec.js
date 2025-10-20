import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('moleculeState', () => {
  test('versions increment on change', () => {
    const st = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
    });
    const v0 = st.versions.positions;
    st.positions[0].x = 0.2;
    st.markPositionsChanged();
    expect(st.versions.positions).toBe(v0 + 1);
  });
});
