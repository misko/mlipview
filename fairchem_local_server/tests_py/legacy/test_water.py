"""Consolidated water single-point test (HTTP app; httpx client)."""

import httpx
import pytest
from ase.build import molecule

# keep encode side-effect from server module
from fairchem_local_server.server import encode


@pytest.mark.timeout(180)
@pytest.mark.asyncio
async def test_water_energy_forces_stress(ray_serve_app):
    """
    Calls /serve/simple and validates energy/forces, with optional stress.
    Uses the ray_serve_app fixture to ensure the server is up and to get the base URL.
    """
    atoms = molecule("H2O")
    # Explicit charge/spin to match UMA defaults (and avoid warnings)
    atoms.info.update({"charge": 0, "spin": 1})

    # encode side-effect (ensure json serialization path doesn’t regress)
    encode(atoms)

    # generous cubic box—stress only meaningful if backend provides it
    cell = [
        [15.0, 0.0, 0.0],
        [0.0, 15.0, 0.0],
        [0.0, 0.0, 15.0],
    ]

    payload = {
        "atomic_numbers": [int(z) for z in atoms.get_atomic_numbers()],
        "coordinates": atoms.get_positions().tolist(),
        "properties": ["energy", "forces", "stress"],
        "charge": int(atoms.info.get("charge", 0)),
        "spin_multiplicity": int(atoms.info.get("spin", 1)),
        "cell": cell,
        # leave pbc False; some calculators won’t return stress otherwise
        "pbc": False,
        # default to UMA; set to "lj" if you want the Lennard-Jones baseline instead
        "calculator": "uma",
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.post(f"{ray_serve_app}/serve/simple", json=payload)
    assert resp.status_code == 200, resp.text

    data = resp.json()
    results = data.get("results", {})
    assert isinstance(results, dict), "Response 'results' must be a dict"

    # energy
    energy = results.get("energy")
    assert isinstance(energy, (int, float)), "Energy must be numeric"

    # forces
    forces = results.get("forces")
    assert isinstance(forces, list) and len(forces) == len(
        atoms
    ), "Forces shape mismatch"
    for f in forces:
        assert isinstance(f, list) and len(f) == 3, "Each force must be a 3-vector"
        for x in f:
            assert isinstance(x, (int, float)), "Force components must be numeric"

    # stress (optional): if provided, validate it
    stress = results.get("stress")
    if stress is not None:
        assert (
            isinstance(stress, list) and len(stress) == 6
        ), "Voigt stress must have 6 components"
        for s in stress:
            assert (
                isinstance(s, (int, float)) and s == s
            ), "Stress components must be finite numbers"
