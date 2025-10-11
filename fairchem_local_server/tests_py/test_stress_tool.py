import asyncio

import httpx
import pytest


@pytest.mark.timeout(120)
@pytest.mark.asyncio
async def test_stress_tool_sends_requests(ray_serve_app):
    """Ensure 10 requests for simple/relax/md succeed.

    We don't validate physics here, just 200 responses and JSON shape sanity.
    """
    base = ray_serve_app.rstrip("/")

    # Build canonical payloads matching Pydantic models
    simple = {
        "atomic_numbers": [8, 1, 1],
        "coordinates": [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239987, 0.926627, 0.0],
        ],
        "properties": ["energy", "forces"],
        "charge": 0,
        "spin_multiplicity": 1,
        "cell": [
            [15.0, 0.0, 0.0],
            [0.0, 15.0, 0.0],
            [0.0, 0.0, 15.0],
        ],
        "calculator": "uma",
    }

    relax = {
        "atomic_numbers": simple["atomic_numbers"],
        "coordinates": simple["coordinates"],
        "steps": 1,
        "calculator": "uma",
        "cell": simple["cell"],
        "pbc": [False, False, False],
        "charge": 0,
        "spin_multiplicity": 1,
        "return_trace": False,
    }

    md = {
        "atomic_numbers": simple["atomic_numbers"],
        "coordinates": simple["coordinates"],
        "steps": 1,
        "temperature": 298.0,
        "timestep_fs": 1.0,
        "friction": 0.02,
        "calculator": "uma",
        "cell": simple["cell"],
        "pbc": [False, False, False],
        "charge": 0,
        "spin_multiplicity": 1,
        "return_trajectory": False,
    }

    async with httpx.AsyncClient(timeout=60.0) as client:
        # Simple x10
        tasks = [client.post(base + "/serve/simple", json=simple) for _ in range(10)]
        rs = await asyncio.gather(*tasks)
        for r in rs:
            assert r.status_code == 200, r.text
            j = r.json()
            assert "results" in j

        # Relax x10
        tasks = [client.post(base + "/serve/relax", json=relax) for _ in range(10)]
        rs = await asyncio.gather(*tasks)
        for r in rs:
            assert r.status_code == 200, r.text
            j = r.json()
            assert "final_energy" in j and "positions" in j

        # MD x10
        tasks = [client.post(base + "/serve/md", json=md) for _ in range(10)]
        rs = await asyncio.gather(*tasks)
        for r in rs:
            assert r.status_code == 200, r.text
            j = r.json()
            assert "final_energy" in j and "positions" in j and "velocities" in j
