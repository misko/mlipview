import math
import os
from pathlib import Path

import httpx
import pytest
import torch
from ase.build import molecule
from ase.io import read  # for loading ROY xyz
from ase.optimize import BFGS

try:
    from ase.calculators.lj import LennardJones as _LJ
except Exception:  # pragma: no cover
    _LJ = None


def LennardJonesCalculator():  # shim to match previous test usage
    if _LJ is None:
        raise RuntimeError("ASE LennardJones unavailable")
    return _LJ(rc=3.0)


STEPS = 10  # keep short for test runtime; adjust if needed
FMAX_TOL = 1e-4
ENERGY_TOL = 1e-6  # direct vs API should be nearly identical

# Extended test constants for ROY long relaxation
ROY_STEPS = (
    int(os.environ.get("MLIPVIEW_ROY_STEPS", "400"))
    if "MLIPVIEW_ROY_STEPS" in os.environ
    else 400
)
ROY_INIT_ENERGY_ABS_TOL = 5e-3  # UMA model slight nondeterminism tolerance (initial)
ROY_FINAL_ENERGY_ABS_TOL = 1e-2  # final energy tolerance after 400 steps


def load_roy_atoms():
    """Load ROY molecule from the repository xyz file by ascending to project root."""
    here = Path(__file__).resolve()
    # Try ascending up to 5 levels to find 'public/molecules/roy.xyz'
    for up in range(1, 6):
        candidate_root = here.parents[up]
        xyz_path = candidate_root / "public" / "molecules" / "roy.xyz"
        if xyz_path.is_file():
            return read(xyz_path)
    raise FileNotFoundError("ROY xyz not found via ascending search from test file")


def run_direct_bfgs(atoms, steps: int, fmax_stop: float | None = None):
    """Run BFGS for up to `steps`, optionally stopping early if max force < fmax_stop.

    Returns (energies, positions, forces).
    """
    energies = [float(atoms.get_potential_energy())]
    opt = BFGS(atoms, logfile=None)
    for _ in range(steps):
        opt.step()
        energies.append(float(atoms.get_potential_energy()))
        if fmax_stop is not None:
            # max force magnitude per-atom
            fmax = (atoms.get_forces() ** 2).sum(axis=1) ** 0.5
            if float(fmax.max()) < fmax_stop:
                break
    return (energies, atoms.get_positions().tolist(), atoms.get_forces().tolist())


@pytest.mark.asyncio
async def test_relax_lj_water(ray_serve_app):
    atoms = molecule("H2O")
    atoms.calc = LennardJonesCalculator()
    _energies, _positions, _forces = run_direct_bfgs(atoms, STEPS)

    # Build API payload from a fresh copy of H2O at its initial geometry
    atoms2 = molecule("H2O")
    payload = {
        "atomic_numbers": atoms2.get_atomic_numbers().tolist(),
        "coordinates": atoms2.get_positions().tolist(),
        "steps": STEPS,
        "calculator": "lj",
    }

    async with httpx.AsyncClient(timeout=60.0) as client:
        resp = await client.post(f"{ray_serve_app}/serve/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Direct initial energy for fresh geometry (LJ)
    atoms_init = molecule("H2O")
    atoms_init.calc = LennardJonesCalculator()
    e0_direct = float(atoms_init.get_potential_energy())
    assert math.isclose(
        data["initial_energy"], e0_direct, rel_tol=0, abs_tol=ENERGY_TOL
    )

    # Re-run direct on fresh atoms to match step count for final comparison
    atoms_cmp = molecule("H2O")
    atoms_cmp.calc = LennardJonesCalculator()
    energies_cmp, pos_cmp, _ = run_direct_bfgs(atoms_cmp, STEPS)
    assert math.isclose(data["final_energy"], energies_cmp[-1], rel_tol=0, abs_tol=1e-5)

    # Position comparison
    for a, b in zip(data["positions"], pos_cmp):
        for x, y in zip(a, b):
            assert abs(x - y) < 1e-5


@pytest.mark.asyncio
async def test_relax_uma_water(ray_serve_app, local_uma_calculator):
    # UMA model may introduce stochasticity; fix torch seed
    torch.manual_seed(0)

    # Direct UMA relaxation baseline
    atoms = molecule("H2O")
    atoms.calc = local_uma_calculator
    energies_direct, pos_direct, _ = run_direct_bfgs(atoms, STEPS)

    # API payload from fresh H2O
    atoms2 = molecule("H2O")
    payload = {
        "atomic_numbers": atoms2.get_atomic_numbers().tolist(),
        "coordinates": atoms2.get_positions().tolist(),
        "steps": STEPS,
        "calculator": "uma",
    }

    async with httpx.AsyncClient(timeout=60.0) as client:
        resp = await client.post(f"{ray_serve_app}/serve/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Initial energy parity (allow small absolute tolerance)
    atoms_init = molecule("H2O")
    atoms_init.calc = local_uma_calculator
    e0_direct = float(atoms_init.get_potential_energy())
    assert math.isclose(data["initial_energy"], e0_direct, rel_tol=0, abs_tol=5e-4)

    # Final energy parity (slightly looser tolerance)
    assert math.isclose(
        data["final_energy"], energies_direct[-1], rel_tol=1e-4, abs_tol=1e-3
    )

    # Positions roughly close
    for a, b in zip(data["positions"], pos_direct):
        for x, y in zip(a, b):
            assert abs(x - y) < 2e-3


@pytest.mark.asyncio
async def test_relax_uma_roy_400_steps_energy(ray_serve_app, local_uma_calculator):
    """Long (400-step) UMA relaxation parity test for ROY."""
    torch.manual_seed(0)

    atoms = load_roy_atoms()
    atoms.calc = local_uma_calculator
    direct_energies, _pos_direct, _forces_direct = run_direct_bfgs(
        atoms, ROY_STEPS, fmax_stop=1e-3
    )

    # Fresh atoms for initial energy (direct) reference
    atoms_init = load_roy_atoms()
    atoms_init.calc = local_uma_calculator
    e0_direct = float(atoms_init.get_potential_energy())

    # API payload from fresh geometry
    atoms_api = load_roy_atoms()
    payload = {
        "atomic_numbers": atoms_api.get_atomic_numbers().tolist(),
        "coordinates": atoms_api.get_positions().tolist(),
        "steps": ROY_STEPS,
        "calculator": "uma",
    }

    async with httpx.AsyncClient(timeout=None) as client:
        resp = await client.post(f"{ray_serve_app}/serve/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Initial energy parity
    assert math.isclose(
        data["initial_energy"], e0_direct, rel_tol=0, abs_tol=ROY_INIT_ENERGY_ABS_TOL
    ), f"ROY initial energy mismatch: api={data['initial_energy']}, direct={e0_direct}"

    # Final energy parity
    api_final = data["final_energy"]
    direct_final = direct_energies[-1]
    assert math.isclose(
        api_final, direct_final, rel_tol=1e-5, abs_tol=ROY_FINAL_ENERGY_ABS_TOL
    ), f"ROY final energy mismatch: api={api_final}, direct={direct_final}"

    # --- Optional: sequential 1-step API sampling to obtain per-step energy trace ---
    step_atoms = load_roy_atoms()
    seq_payload = {
        "atomic_numbers": step_atoms.get_atomic_numbers().tolist(),
        "calculator": "uma",
        "steps": 1,
    }
    step_positions = step_atoms.get_positions().tolist()
    api_step_energies = []
    async with httpx.AsyncClient(timeout=None) as client:
        for i in range(min(ROY_STEPS, 120)):
            seq_payload["coordinates"] = step_positions
            r = await client.post(f"{ray_serve_app}/serve/relax", json=seq_payload)
            assert (
                r.status_code == 200
            ), f"step {i} status {r.status_code} body={r.text}"
            jd = r.json()
            api_step_energies.append(jd["final_energy"])
            step_positions = jd["positions"]

    # Debug prints for inspection (pytest -s)
    print("\nROY_API_STEP_ENERGIES", api_step_energies)
    print(
        "ROY_DIRECT_FINAL",
        direct_final,
        "ROY_API_AGG_FINAL",
        api_final,
        "SEQ_LEN",
        len(api_step_energies),
    )

    if len(api_step_energies) > 2:
        diffs = [
            api_step_energies[i + 1] - api_step_energies[i]
            for i in range(len(api_step_energies) - 1)
        ]
        up = sum(1 for d in diffs if d > 0)
        down = sum(1 for d in diffs if d < 0)
        print(f"ROY_API_TRACE_DIFFS up={up} down={down} totalSteps={len(diffs)}")


@pytest.mark.asyncio
async def test_relax_uma_water_bfgs_linesearch_trace(ray_serve_app):
    """Verify plain BFGS optimizer selection + energy trace return for small system."""
    torch.manual_seed(0)
    atoms = molecule("H2O")
    payload = {
        "atomic_numbers": atoms.get_atomic_numbers().tolist(),
        "coordinates": atoms.get_positions().tolist(),
        "steps": 5,
        "calculator": "uma",
        "optimizer": "bfgs",  # server implements plain BFGS
        "optimizer_params": {"maxstep": 0.15},
        "return_trace": True,
        "fmax": 0.0,
    }

    async with httpx.AsyncClient(timeout=60.0) as client:
        resp = await client.post(f"{ray_serve_app}/serve/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    trace = data.get("trace_energies")
    assert trace is not None
    assert len(trace) == data["steps_completed"]
    assert 1 <= data["steps_completed"] <= 5
    for e in trace:
        assert isinstance(e, (float, int))
    assert data["initial_energy"] < 0 and data["final_energy"] < 0


@pytest.mark.asyncio
async def test_relax_uma_roy_400_steps_bfgs_linesearch_trace(ray_serve_app):
    """Run ROY long UMA relaxation with BFGS and capture trace (sanity checks only)."""
    torch.manual_seed(0)
    steps = ROY_STEPS

    atoms_api = load_roy_atoms()
    payload = {
        "atomic_numbers": atoms_api.get_atomic_numbers().tolist(),
        "coordinates": atoms_api.get_positions().tolist(),
        "steps": steps,
        "calculator": "uma",
        "optimizer": "bfgs",
        "optimizer_params": {"maxstep": 0.15},
        "return_trace": True,
        "fmax": 0.0,
    }

    async with httpx.AsyncClient(timeout=None) as client:
        resp = await client.post(f"{ray_serve_app}/serve/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    trace = data.get("trace_energies") or []
    print(
        f"ROY_LS_TRACE_LEN {len(trace)} steps_completed {data['steps_completed']} requested {steps}"
    )
    if trace:
        diffs = [trace[i + 1] - trace[i] for i in range(len(trace) - 1)]
        up = sum(1 for d in diffs if d > 0)
        down = sum(1 for d in diffs if d < 0)
        print(
            f"ROY_LS_TRACE_DIFFS up={up} down={down} total={len(diffs)} "
            f"firstE={trace[0]} lastE={trace[-1]}"
        )
        assert len(trace) == data["steps_completed"]
        assert data["steps_completed"] <= steps
        assert isinstance(trace[0], (float, int)) and isinstance(
            trace[-1], (float, int)
        )
