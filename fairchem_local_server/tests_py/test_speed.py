"""Simple throughput/latency sanity test for /serve/simple."""

import asyncio
import math
import os
import statistics
import time

import httpx
import pytest
from ase.build import molecule

# Tunables via env to avoid flakiness in CI
REQUESTS = int(os.environ.get("SPEED_N", "16"))  # total requests
CONCURRENCY = int(os.environ.get("SPEED_CONC", "8"))  # in-flight requests
TIMEOUT_S = float(os.environ.get("SPEED_TIMEOUT", "60"))
# Sanity ceiling per-request (very generous to avoid flakes on slow CI)
MAX_LATENCY_S = float(os.environ.get("SPEED_MAX_LAT", "10"))


def _payload_for_water():
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})
    return {
        "atomic_numbers": [int(z) for z in atoms.get_atomic_numbers()],
        "coordinates": atoms.get_positions().tolist(),
        "properties": ["energy", "forces"],  # keep it lighter than stress
        "charge": 0,
        "spin_multiplicity": 1,
        "cell": None,
        "pbc": False,
        "calculator": "uma",
    }


async def _worker(base_url: str, n: int):
    """Send n sequential requests and collect latencies."""
    payload = _payload_for_water()
    latencies = []
    ok = 0
    async with httpx.AsyncClient(timeout=TIMEOUT_S) as client:
        for _ in range(n):
            t0 = time.perf_counter()
            r = await client.post(f"{base_url}/serve/simple", json=payload)
            dt = time.perf_counter() - t0
            latencies.append(dt)
            if r.status_code == 200:
                data = r.json()
                res = data.get("results", {})
                # minimal correctness check
                e = res.get("energy", None)
                ok += int(isinstance(e, (int, float)) and math.isfinite(e))
            else:
                # still record latency; response code validated by the caller
                pass
    return ok, latencies


@pytest.mark.timeout(300)
@pytest.mark.asyncio
async def test_simple_calculate_speed(ray_serve_app):
    """
    Fires REQUESTS /serve/simple calls with CONCURRENCY workers.
    Asserts:
      - All responses are 200 with numeric energy
      - 95th percentile latency under MAX_LATENCY_S (very lenient by default)
    Prints a small timing summary for visibility in CI logs.
    """
    # Split the total across workers
    per_worker = [REQUESTS // CONCURRENCY] * CONCURRENCY
    for i in range(REQUESTS % CONCURRENCY):
        per_worker[i] += 1

    tasks = [_worker(ray_serve_app, n) for n in per_worker if n > 0]
    results = await asyncio.gather(*tasks)

    all_lat = []
    ok_total = 0
    for ok, lat in results:
        ok_total += ok
        all_lat.extend(lat)

    assert (
        ok_total == REQUESTS
    ), f"{REQUESTS - ok_total} requests failed correctness checks"

    # All should be 200 + numeric energy, but we also bound latency to catch regressions
    all_lat.sort()
    p95 = all_lat[int(0.95 * (len(all_lat) - 1))]
    avg = statistics.mean(all_lat) if all_lat else 0.0
    p50 = all_lat[len(all_lat) // 2] if all_lat else 0.0
    print(
        f"[/serve/simple] n={REQUESTS} conc={CONCURRENCY} "
        f"avg={avg*1000:.1f}ms p50={p50*1000:.1f}ms p95={p95*1000:.1f}ms"
    )
    assert p95 < MAX_LATENCY_S, f"p95 {p95:.2f}s exceeds cap {MAX_LATENCY_S:.2f}s"
