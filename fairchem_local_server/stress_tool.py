"""
Async stress/load test utility for UMA Serve API.

Sends a configurable mix of /serve/simple, /serve/relax, and /serve/md requests
at a target base URL, measuring throughput and latency statistics.

Usage examples:
  python -m fairchem_local_server.stress_tool --base http://127.0.0.1:8000 \
    --concurrency 32 --duration 60 --mix 6:2:2 --relax-steps 10 \
    --md-steps 10

All defaults are conservative; tune mix, steps, and concurrency to your
machine.
"""

from __future__ import annotations

import argparse
import asyncio
import os
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import httpx

# --- Simple built-in molecules (atomic_numbers, coordinates) -----------------


def _water() -> Tuple[List[int], List[List[float]]]:
    # Roughly centered gas-phase H2O
    Z = [8, 1, 1]
    xyz = [
        [0.000000, 0.000000, 0.000000],
        [0.957200, 0.000000, 0.000000],
        [-0.239987, 0.926627, 0.000000],
    ]
    return Z, xyz


# --- Request builders --------------------------------------------------------


def _default_molecules_dir() -> str:
    # Project root = parent of this file's directory
    root = os.path.dirname(os.path.dirname(__file__))
    return os.path.join(root, "public", "molecules")


def _load_molecules(
    mol_dir: str,
) -> List[Tuple[List[int], List[List[float]], str]]:
    """Load .xyz molecules from a directory via ASE: returns (Z, xyz, name)."""
    mols: List[Tuple[List[int], List[List[float]], str]] = []
    try:
        from ase.io import read  # type: ignore
    except Exception:
        # ASE not available: fallback to built-in water only
        Z, xyz = _water()
        mols.append((Z, xyz, "water"))
        return mols

    if not os.path.isdir(mol_dir):
        Z, xyz = _water()
        mols.append((Z, xyz, "water"))
        return mols

    for fn in os.listdir(mol_dir):
        if not fn.lower().endswith(".xyz"):
            continue
        path = os.path.join(mol_dir, fn)
        try:
            atoms = read(path)
            Z = [int(z) for z in atoms.get_atomic_numbers()]
            xyz = atoms.get_positions().tolist()
            name = os.path.splitext(fn)[0]
            mols.append((Z, xyz, name))
        except Exception:
            # Skip unreadable files
            continue
    if not mols:
        Z, xyz = _water()
        mols.append((Z, xyz, "water"))
    return mols


def build_payloads_from_mol(
    Z: List[int], xyz: List[List[float]], *, relax_steps: int, md_steps: int
) -> Dict[str, dict]:
    base = {
        "atomic_numbers": Z,
        "coordinates": xyz,
        "calculator": "uma",
        "charge": 0,
        "spin_multiplicity": 1,
        "cell": [
            [15.0, 0.0, 0.0],
            [0.0, 15.0, 0.0],
            [0.0, 0.0, 15.0],
        ],
        # PBC must be a list[bool,bool,bool] for RelaxIn/MDIn
        # "pbc": [False, False, False],
    }
    return {
        "simple": {
            **{k: base[k] for k in ("atomic_numbers", "coordinates", "cell")},
            "properties": ["energy", "forces"],
            "charge": base["charge"],
            "spin_multiplicity": base["spin_multiplicity"],
            "calculator": "uma",
        },
        "relax": {
            **base,
            "steps": int(relax_steps),
            "return_trace": False,
        },
        "md": {
            **base,
            "steps": int(md_steps),
            "temperature": 298.0,
            "timestep_fs": 1.0,
            "friction": 0.02,
            "return_trajectory": False,
        },
    }


# --- Runner -----------------------------------------------------------------


@dataclass
class Metrics:
    count_ok: int = 0
    count_err: int = 0
    latencies_ms: List[float] | None = None

    def __post_init__(self):
        if self.latencies_ms is None:
            self.latencies_ms = []

    def record(self, ms: float, ok: bool):
        self.latencies_ms.append(ms)
        if ok:
            self.count_ok += 1
        else:
            self.count_err += 1


async def _one(client: httpx.AsyncClient, url: str, payload: dict, m: Metrics):
    t0 = time.perf_counter()
    ok = False
    try:
        r = await client.post(url, json=payload)
        ok = r.status_code == 200
    except Exception:
        ok = False
    finally:
        dt = (time.perf_counter() - t0) * 1000.0
        m.record(dt, ok)


async def run_stress(
    base_url: str,
    duration_s: int,
    concurrency: int,
    mix: Tuple[int, int, int],
    relax_steps: int,
    md_steps: int,
    route_style: str = "serve",
    molecules_dir: str | None = None,
):
    # Build payload pool from available molecules
    mol_dir = molecules_dir or _default_molecules_dir()
    mols = _load_molecules(mol_dir)
    payload_pool = [
        build_payloads_from_mol(Z, xyz, relax_steps=relax_steps, md_steps=md_steps)
        for (Z, xyz, _name) in mols
    ]
    if route_style == "plain":
        endpoints = {
            "simple": "/simple_calculate",
            "relax": "/relax",
            "md": "/md",
        }
    else:
        endpoints = {
            "simple": "/serve/simple",
            "relax": "/serve/relax",
            "md": "/serve/md",
        }

    weights = [max(0, mix[0]), max(0, mix[1]), max(0, mix[2])]
    kinds = ["simple", "relax", "md"]
    choices = [k for k, w in zip(kinds, weights) for _ in range(int(w) or 0)] or [
        "simple"
    ]

    m = Metrics()
    started = time.time()
    deadline = started + float(duration_s)
    sem = asyncio.Semaphore(max(1, int(concurrency)))
    # Periodic reporting state
    report_every = 5.0
    last_report_t = started
    last_report_ok = 0

    async def worker():
        async with sem:
            kind = random.choice(choices)
            url = base_url.rstrip("/") + endpoints[kind]
            pset = random.choice(payload_pool)
            await _one(client, url, pset[kind], m)

    async with httpx.AsyncClient(timeout=120.0) as client:
        tasks: List[asyncio.Task] = []
        while time.time() < deadline:
            # Maintain roughly 'concurrency' in-flight tasks
            # Clean up any finished tasks
            tasks = [t for t in tasks if not t.done()]
            while len(tasks) < concurrency:
                tasks.append(asyncio.create_task(worker()))
            # Periodic stats output every ~5s
            now = time.time()
            if now - last_report_t >= report_every:
                interval = now - last_report_t
                d_ok = m.count_ok - last_report_ok
                rps = d_ok / interval if interval > 0 else 0.0
                lat = list(m.latencies_ms)
                p50 = _pct(lat, 0.50)
                p90 = _pct(lat, 0.90)
                elapsed = now - started
                print(
                    (
                        "[stress] t+%.1fs ok=%d err=%d rps~%.2f "
                        "lat_ms p50=%.1f p90=%.1f"
                    )
                    % (elapsed, m.count_ok, m.count_err, rps, p50, p90)
                )
                last_report_t = now
                last_report_ok = m.count_ok
                # no need to track last_report_err, count_err included above
            await asyncio.sleep(0.005)
        # Drain remaining
        if tasks:
            await asyncio.gather(*tasks, return_exceptions=True)

    return m


def _pct(values: List[float], p: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    idx = max(0, min(len(s) - 1, int(p * (len(s) - 1))))
    return s[idx]


def main():
    ap = argparse.ArgumentParser(description="UMA Serve stress utility")
    ap.add_argument("--base", default="http://127.0.0.1:8000", help="Base URL")
    ap.add_argument("--duration", type=int, default=30, help="Duration seconds")
    ap.add_argument("--concurrency", type=int, default=16, help="Concurrent tasks")
    ap.add_argument(
        "--mix",
        default="6:2:2",
        help="Ratio simple:relax:md, e.g. 6:2:2",
    )
    ap.add_argument("--relax-steps", type=int, default=5)
    ap.add_argument("--md-steps", type=int, default=5)
    ap.add_argument(
        "--route",
        choices=["serve", "plain"],
        default="serve",
        help="Use Serve ingress routes (/serve/*) or plain FastAPI routes",
    )
    ap.add_argument(
        "--molecules-dir",
        default=None,
        help=(
            "Directory with .xyz files to sample (default: project " "public/molecules)"
        ),
    )
    args = ap.parse_args()

    try:
        parts = [int(x) for x in str(args.mix).split(":")]
        while len(parts) < 3:
            parts.append(0)
        mix = (parts[0], parts[1], parts[2])
    except Exception:
        mix = (6, 2, 2)

    m: Metrics = asyncio.run(
        run_stress(
            base_url=args.base,
            duration_s=args.duration,
            concurrency=args.concurrency,
            mix=mix,
            relax_steps=args.relax_steps,
            md_steps=args.md_steps,
            route_style=args.route,
            molecules_dir=args.molecules_dir,
        )
    )

    total = m.count_ok + m.count_err
    duration = max(1.0, float(args.duration))
    rps = m.count_ok / duration
    lat = m.latencies_ms
    p50 = _pct(lat, 0.50)
    p90 = _pct(lat, 0.90)
    p99 = _pct(lat, 0.99)
    avg = statistics.mean(lat) if lat else 0.0

    print("\n=== UMA Serve Stress Results ===")
    print(f"base= {args.base} route_style= {args.route}")
    print(f"duration_s= {duration:.1f} concurrency= {args.concurrency}")
    print(f"req_total= {total} ok= {m.count_ok} err= {m.count_err}")
    print(f"throughput_rps= {rps:.2f}")
    print("latency_ms avg= %.1f p50= %.1f p90= %.1f p99= %.1f" % (avg, p50, p90, p99))


if __name__ == "__main__":
    main()
