import math

import httpx
import pytest


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
@pytest.mark.asyncio
async def test_md_water_uma_stability(ray_serve_app):
    """Run a 500-step 298K MD on water with UMA; ensure system doesn't explode.

    Criteria:
    - Energies finite (not NaN/Inf)
    - Max pair distance stays below 10 Å
    - Returned steps_completed == requested
    - Final temperature within 50%–200% of target
    """
    atomic_numbers = [8, 1, 1]
    coordinates = [
        [0.0, 0.0, 0.0],  # O at origin
        [0.9575, 0.0, 0.0],  # H
        [-0.2399872, 0.92662721, 0.0],  # H
    ]
    payload = {
        "atomic_numbers": [int(z) for z in atomic_numbers],
        "coordinates": coordinates,
        "steps": 500,
        "temperature": 298.0,
        "timestep_fs": 1.0,
        "calculator": "uma",
    }

    async with httpx.AsyncClient(timeout=60.0) as client:
        r = await client.post(f"{ray_serve_app}/serve/md", json=payload)

    assert r.status_code == 200, r.text
    data = r.json()

    # steps
    assert data["steps_completed"] == payload["steps"]

    # energies finite
    assert math.isfinite(data["initial_energy"]) and math.isfinite(data["final_energy"])

    # distance check
    max_dist = _max_pair_distance(data["positions"])
    assert max_dist < 10.0, f"Unphysical expansion (max pair distance {max_dist:.2f} Å)"

    # temperature band (coarse thermostat)
    T_final = data["temperature"]
    tgt = payload["temperature"]
    assert (
        (0.5 * tgt) < T_final < (2.0 * tgt + 0.1)
    ), f"Final T {T_final:.1f}K out of band"
