// Test that cell centering now applies an originOffset without moving atoms (no jump)
import { buildMolecule } from '../public/molecule.js';

// BABYLON mock provided globally via jest.setup.js
const BABYLON = global.BABYLON;

describe('cell centering originOffset', () => {
  test('originOffset applied without moving atoms', () => {
    const atoms = [
      { element: 'C', pos: new BABYLON.Vector3(0,0,0) },
      { element: 'C', pos: new BABYLON.Vector3(2,0,0) },
      { element: 'C', pos: new BABYLON.Vector3(0,3,0) }
    ];
    // Save original positions references to compare after setCellVectors
    const originals = atoms.map(a => a.pos.clone());

    const mol = buildMolecule({}, { atoms, center: false });
    // Provide minimal groups & bondGroups expectations already set by buildMolecule
    // Define a simple orthorhombic cell bigger than the atom extents so offset non-zero
    mol.setCellVectors(new BABYLON.Vector3(10,0,0), new BABYLON.Vector3(0,10,0), new BABYLON.Vector3(0,0,10), { recenter: true });

    // Atoms should remain unchanged
    atoms.forEach((a,i) => {
      expect(a.pos.x).toBeCloseTo(originals[i].x, 10);
      expect(a.pos.y).toBeCloseTo(originals[i].y, 10);
      expect(a.pos.z).toBeCloseTo(originals[i].z, 10);
    });

    // originOffset should be non-zero if centroid not already at 0.5,0.5,0.5 in fractional space
    const off = mol.__cellState.originOffset;
    expect(off).toBeDefined();
    expect(Math.abs(off.x)+Math.abs(off.y)+Math.abs(off.z)).toBeGreaterThan(0);

    // Toggling cell should not mutate atom positions
    mol.toggleCell(true);
    atoms.forEach((a,i) => {
      expect(a.pos.x).toBeCloseTo(originals[i].x, 10);
    });
  });
});
