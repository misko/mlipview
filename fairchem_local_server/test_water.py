import json

import pytest
from ase.build import molecule
from fastapi.testclient import TestClient

from fairchem_local_server.server import app, encode

client = TestClient(app)


@pytest.mark.timeout(180)
def test_water_energy_forces_stress():
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})
    encoded = encode(atoms)
    if isinstance(encoded, str):
        atoms_json = json.loads(encoded)
    else:
        atoms_json = encoded
    # Provide a pseudo cell so stress is meaningful.
    # Large orthorhombic box to avoid strong PBC interactions.
    cell = [
        [15.0, 0.0, 0.0],
        [0.0, 15.0, 0.0],
        [0.0, 0.0, 15.0],
    ]
    payload = {
        "atoms_json": atoms_json,
        "properties": ["energy", "forces", "stress"],
        "info": atoms.info,
        # For /calculate path the server reconstructs Atoms; cell needs to
        # be baked into atoms_json via encode (future: extend API). Keeping
        # placeholder here for clarity if API expands to accept cell.
        "cell": cell,
    }
    resp = client.post("/calculate", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "results" in data
    results = data["results"]
    assert "energy" in results
    assert isinstance(results["energy"], (int, float))
    assert "forces" in results
    forces = results["forces"]
    assert isinstance(forces, list) and len(forces) == len(atoms)
    # each force vector should have length 3
    for f in forces:
        assert isinstance(f, list) and len(f) == 3
    # Stress may be None; if present validate length 6 (Voigt order)
    if "stress" in results and results["stress"] is not None:
        stress = results["stress"]
        assert isinstance(stress, list) and len(stress) == 6
        # Basic numeric sanity (no NaN)
        for s in stress:
            assert isinstance(s, (int, float)) and s == s
