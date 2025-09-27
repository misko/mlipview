import json

import pytest
from ase.build import molecule
from fastapi.testclient import TestClient

from fairchem_local_server.server import app, encode

client = TestClient(app)


@pytest.mark.timeout(120)
def test_water_energy_forces():
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})
    encoded = encode(atoms)
    if isinstance(encoded, str):
        atoms_json = json.loads(encoded)
    else:
        atoms_json = encoded
    payload = {
        "atoms_json": atoms_json,
        "properties": ["energy", "forces"],
        "info": atoms.info,
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
