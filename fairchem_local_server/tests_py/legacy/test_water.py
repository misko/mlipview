"""Consolidated water single-point test (merged legacy duplicates)."""

import pytest
from ase.build import molecule
from fastapi.testclient import TestClient

from fairchem_local_server.server import app, encode

client = TestClient(app)


@pytest.mark.timeout(180)
def test_water_energy_forces_stress():
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})
    # encode side-effect (ensures no regression in json serialization path)
    encode(atoms)
    cell = [
        [15.0, 0.0, 0.0],
        [0.0, 15.0, 0.0],
        [0.0, 0.0, 15.0],
    ]
    payload = {
        "atomic_numbers": atoms.get_atomic_numbers().tolist(),
        "coordinates": atoms.get_positions().tolist(),
        "properties": ["energy", "forces", "stress"],
        "charge": atoms.info.get("charge", 0),
        "spin_multiplicity": atoms.info.get("spin", 1),
        "cell": cell,
    }
    resp = client.post("/simple_calculate", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    results = data.get("results", {})
    energy = results.get("energy")
    assert isinstance(energy, (int, float))
    forces = results.get("forces")
    assert isinstance(forces, list) and len(forces) == len(atoms)
    for f in forces:
        assert isinstance(f, list) and len(f) == 3
    stress = results.get("stress")
    if stress is not None:
        assert isinstance(stress, list) and len(stress) == 6
        for s in stress:
            assert isinstance(s, (int, float)) and s == s
