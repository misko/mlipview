# tests_py/test_serve_batches_simple.py

import asyncio
import copy
import math
import time
from time import perf_counter

import httpx
import pytest
from ase.build import molecule


@pytest.mark.timeout(240)
@pytest.mark.asyncio
async def test_serve_batches_simple_requests(ray_serve_app):
    """
    Evidence-of-batching test:
      - Send N identical /serve/simple requests sequentially, measure total time.
      - Send the same N requests concurrently, measure total time.
      - With a single UMA replica, a meaningful speedup in the concurrent case
        implies Serve batched (coalesced) them into fewer GPU calls.
    """
    N = 20  # keep <= default UMA_BATCH_MAX(16)
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})

    payload = {
        "atomic_numbers": [int(z) for z in atoms.get_atomic_numbers()],
        "coordinates": atoms.get_positions().tolist(),
        "properties": ["energy", "forces"],  # non-trivial request
        "charge": int(atoms.info.get("charge", 0)),
        "spin_multiplicity": int(atoms.info.get("spin", 1)),
        "cell": [
            [15.0, 0.0, 0.0],
            [0.0, 15.0, 0.0],
            [0.0, 0.0, 15.0],
        ],
        "pbc": True,
        "calculator": "uma",
    }

    base = ray_serve_app.rstrip("/")
    url = f"{base}/serve/simple"

    async with httpx.AsyncClient(timeout=60.0) as client:
        # Warm-up: first request often includes model/graph init & CUDA warmup
        warm = await client.post(url, json=payload)
        assert warm.status_code == 200, warm.text

        payloads = []
        for _ in range(N):
            # Slightly vary input to avoid any caching effects
            payload["coordinates"][0][0] += 1
            payloads.append(copy.deepcopy(payload))
        # --- Sequential N requests ---
        t0 = perf_counter()
        for i in range(N):
            r = await client.post(url, json=payloads[i])
            assert r.status_code == 200, r.text
            out = r.json().get("results", {})
            # sanity: energy numeric, forces shaped correctly
            assert isinstance(out.get("energy"), (int, float)) and math.isfinite(
                out["energy"]
            )
            forces = out.get("forces", [])
            assert isinstance(forces, list) and len(forces) == len(atoms)
        t_seq = perf_counter() - t0

        # --- Concurrent N requests (expected to batch) ---
        t1 = perf_counter()
        tasks = [client.post(url, json=payloads[i]) for i in range(N)]
        rs = await asyncio.gather(*tasks)
        t_conc = perf_counter() - t1

    # All responses must be OK
    for r in rs:
        assert r.status_code == 200, r.text
        out = r.json().get("results", {})
        assert isinstance(out.get("energy"), (int, float)) and math.isfinite(
            out["energy"]
        )

    # Heuristic batching confirmation:
    # With batching, concurrent wall-time should be much less than sequential.
    # We require at least ~2x speedup to avoid flakiness across environments.
    speedup = t_seq / t_conc if t_conc > 0 else float("inf")
    # Helpful log for CI output
    print(
        f"\n[BATCHING] sequential={t_seq:.3f}s concurrent={t_conc:.3f}s speedup={speedup:.2f}x\n"
    )

    assert speedup >= 2.0, (
        "Concurrent requests did not achieve the expected speedup; "
        "this suggests Serve may not be batching. "
        f"(seq={t_seq:.3f}s conc={t_conc:.3f}s speedup={speedup:.2f}x)"
    )
