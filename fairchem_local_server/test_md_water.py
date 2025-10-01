import math

import pytest
from fastapi.testclient import TestClient

from fairchem_local_server.server import app

client = TestClient(app)


def _max_pair_distance(coords):
    m = 0.0
    for i in range(len(coords)):
        xi, yi, zi = coords[i]
        for j in range(i + 1, len(coords)):
            xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d > m:
                m = d
    return m


@pytest.mark.timeout(300)
def test_md_water_uma_stability():
    """Run a 500-step 298K MD on water with UMA; ensure system doesn't explode.

    Criteria:
    - Energies finite (not NaN/Inf)
    - Max pair distance stays below 10 Å (far above typical water bond but below explosion)
    - Returned steps_completed == requested
    - Final temperature within a reasonable band of target (50% - 200%) acknowledging crude rescaling.
    """
    atomic_numbers = [8, 1, 1]
    coordinates = [
        [0.0, 0.0, 0.0],  # O at origin
        [0.9575, 0.0, 0.0],  # H
        [-0.2399872, 0.92662721, 0.0],  # H
    ]
    payload = {
        "atomic_numbers": atomic_numbers,
        "coordinates": coordinates,
        "steps": 500,
        "temperature": 298.0,
        "timestep_fs": 1.0,
        "calculator": "uma",
    }
    resp = client.post("/md", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["steps_completed"] == 500
    assert math.isfinite(data["initial_energy"]) and math.isfinite(data["final_energy"])
    # distance check
    max_dist = _max_pair_distance(data["positions"])
    assert max_dist < 10.0, f"Unphysical expansion (max pair distance {max_dist:.2f} Å)"
    # temperature band (coarse thermostat)
    T_final = data["temperature"]
    assert 0.5 * payload["temperature"] < T_final < 2.0 * payload["temperature"], f"Final T {T_final:.1f}K out of band"
