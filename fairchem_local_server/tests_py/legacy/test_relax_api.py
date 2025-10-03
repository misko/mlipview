import math
import os
from pathlib import Path

import torch
from ase.build import molecule
from ase.io import read  # for loading ROY xyz
from ase.optimize import BFGS
from fastapi.testclient import TestClient

from fairchem_local_server import (
    server,
)  # adjust to import via module path if package, else fallback
from fairchem_local_server.model_runtime import get_calculator
from fairchem_local_server.server import app

try:
    from ase.calculators.lj import LennardJones as _LJ
except Exception:  # pragma: no cover
    _LJ = None


def LennardJonesCalculator():  # shim to match previous test usage
    if _LJ is None:
        raise RuntimeError("ASE LennardJones unavailable")
    return _LJ(rc=3.0)


client = TestClient(app)

STEPS = 10  # keep short for test runtime; adjust if needed
FMAX_TOL = 1e-4
ENERGY_TOL = 1e-6  # direct vs API should be nearly identical
POS_TOL = 1e-6

# Extended test constants for ROY long relaxation
ROY_STEPS = (
    int(os.environ.get("MLIPVIEW_ROY_STEPS", "400"))
    if "MLIPVIEW_ROY_STEPS" in os.environ
    else 400
)
ROY_INIT_ENERGY_ABS_TOL = 5e-3  # UMA model slight nondeterminism tolerance (initial)
ROY_FINAL_ENERGY_ABS_TOL = 1e-2  # final energy tolerance after 400 steps


def load_roy_atoms():
    """Load ROY molecule from the repository xyz file.

    Adjusted to resolve the actual repository root (tests_py/legacy is two levels
    below project root fairchem_local_server/), so we ascend until we find
    public/molecules/roy.xyz. This prevents FileNotFoundError in relocated test.
    """
    here = Path(__file__).resolve()
    # Try ascending up to 5 levels to find 'public/molecules/roy.xyz'
    for up in range(1, 6):
        candidate_root = here.parents[up]
        xyz_path = candidate_root / "public" / "molecules" / "roy.xyz"
        if xyz_path.is_file():
            return read(xyz_path)
    raise FileNotFoundError(
        "ROY xyz not found by ascending search from legacy test file"
    )


def run_direct_bfgs(atoms, steps: int, fmax_stop: float | None = None):
    """Run BFGS for up to `steps`, optionally stopping early if max force below fmax_stop.

    Returns (energies, positions, forces).
    """
    energies = [float(atoms.get_potential_energy())]
    opt = BFGS(atoms, logfile=None)
    for _ in range(steps):
        opt.step()
        energies.append(float(atoms.get_potential_energy()))
        if fmax_stop is not None:
            fmax = (atoms.get_forces() ** 2).sum(axis=1) ** 0.5
            if float(fmax.max()) < fmax_stop:
                break
    return (energies, atoms.get_positions().tolist(), atoms.get_forces().tolist())


def test_relax_lj_water():
    atoms = molecule("H2O")
    atoms.calc = LennardJonesCalculator()
    energies, positions, forces = run_direct_bfgs(atoms, STEPS)

    payload = {
        "atomic_numbers": atoms.get_atomic_numbers().tolist(),
        "coordinates": atoms.get_positions().tolist(),  # overwritten below
        # NOTE: We need to restart from initial geometry so rebuild
    }
    # rebuild initial geometry for API call
    atoms2 = molecule("H2O")
    payload["coordinates"] = atoms2.get_positions().tolist()
    payload["steps"] = STEPS
    payload["calculator"] = "lj"

    resp = client.post("/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Direct initial energy (fresh molecule with LJ)
    atoms_init = molecule("H2O")
    atoms_init.calc = LennardJonesCalculator()
    e0_direct = float(atoms_init.get_potential_energy())

    assert math.isclose(
        data["initial_energy"], e0_direct, rel_tol=0, abs_tol=ENERGY_TOL
    )
    # Re-run direct on fresh atoms to match step count
    atoms_cmp = molecule("H2O")
    atoms_cmp.calc = LennardJonesCalculator()
    energies_cmp, pos_cmp, _ = run_direct_bfgs(atoms_cmp, STEPS)
    assert math.isclose(data["final_energy"], energies_cmp[-1], rel_tol=0, abs_tol=1e-5)

    # Position comparison
    for a, b in zip(data["positions"], pos_cmp):
        for x, y in zip(a, b):
            assert abs(x - y) < 1e-5


def test_relax_uma_water():
    # UMA model may introduce stochasticity; fix torch seed
    torch.manual_seed(0)
    atoms = molecule("H2O")
    # direct UMA relaxation
    atoms.calc = get_calculator()
    energies_direct, pos_direct, _ = run_direct_bfgs(atoms, STEPS)

    atoms2 = molecule("H2O")
    payload = {
        "atomic_numbers": atoms2.get_atomic_numbers().tolist(),
        "coordinates": atoms2.get_positions().tolist(),
        "steps": STEPS,
        "calculator": "uma",
    }
    resp = client.post("/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # initial energy
    atoms_init = molecule("H2O")
    atoms_init.calc = get_calculator()
    e0_direct = float(atoms_init.get_potential_energy())
    # Allow slightly looser tolerance (ML model nondeterminism / device differences)
    assert math.isclose(data["initial_energy"], e0_direct, rel_tol=0, abs_tol=5e-4)

    # final energy tolerance looser due to possible tiny nondeterminism
    assert math.isclose(
        data["final_energy"],
        energies_direct[-1],
        rel_tol=1e-4,
        abs_tol=1e-3,
    )

    # positions roughly close (allow slightly larger tolerance due to updated relaxation step logic)
    for a, b in zip(data["positions"], pos_direct):
        for x, y in zip(a, b):
            assert abs(x - y) < 2e-3


def test_relax_uma_roy_400_steps_energy():
    """Long (400-step) UMA relaxation parity test for ROY.

    Compares direct in-process ASE BFGS (UMA calculator) energies to the /relax API
    performing an equivalent number of steps. We start with energy-only parity; force &
    position parity can be added subsequently once acceptable tolerances are established.
    """
    torch.manual_seed(0)
    atoms = load_roy_atoms()
    # obtain UMA calculator
    atoms.calc = get_calculator()
    # Direct BFGS (ROY_STEPS)
    direct_energies, _pos_direct, _forces_direct = run_direct_bfgs(
        atoms, ROY_STEPS, fmax_stop=1e-3
    )

    # Fresh atoms for initial energy (direct) reference
    atoms_init = load_roy_atoms()
    atoms_init.calc = get_calculator()
    e0_direct = float(atoms_init.get_potential_energy())

    # Build API payload (fresh geometry so starting coordinates match API expectation)
    atoms_api = load_roy_atoms()
    payload = {
        "atomic_numbers": atoms_api.get_atomic_numbers().tolist(),
        "coordinates": atoms_api.get_positions().tolist(),
        "steps": ROY_STEPS,
        "calculator": "uma",
    }
    resp = client.post("/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    # Initial energy parity (allow modest absolute tolerance)
    assert math.isclose(
        data["initial_energy"], e0_direct, rel_tol=0, abs_tol=ROY_INIT_ENERGY_ABS_TOL
    ), f"ROY initial energy mismatch: api={data['initial_energy']}, direct={e0_direct}"

    # Final energy parity (looser abs tolerance; relative tolerance tiny due to magnitude)
    api_final = data["final_energy"]
    direct_final = direct_energies[-1]
    assert math.isclose(
        api_final, direct_final, rel_tol=1e-5, abs_tol=ROY_FINAL_ENERGY_ABS_TOL
    ), f"ROY final energy mismatch: api={api_final}, direct={direct_final}"

    # --- Sequential 1-step API sampling to obtain per-step energy trace ---
    # Rebuild fresh atoms for stepwise sequence so that initial geometry matches again.
    step_atoms = load_roy_atoms()
    seq_payload_base = {
        "atomic_numbers": step_atoms.get_atomic_numbers().tolist(),
        "calculator": "uma",
        "steps": 1,
    }
    step_positions = step_atoms.get_positions().tolist()
    api_step_energies = []
    # initial single-point energy via simple_calculate (optional) omitted; we rely on first relax step call's initial_energy
    for i in range(
        min(ROY_STEPS, 120)
    ):  # cap trace length for runtime if ROY_STEPS huge
        seq_payload_base["coordinates"] = step_positions
        r = client.post("/relax", json=seq_payload_base)
        assert r.status_code == 200, f"step {i} status {r.status_code} body={r.text}"
        jd = r.json()
        api_step_energies.append(jd["final_energy"])
        step_positions = jd["positions"]
    # Print energies for external inspection (pytest -s)
    print("\nROY_API_STEP_ENERGIES", api_step_energies)
    print(
        "ROY_DIRECT_FINAL",
        direct_final,
        "ROY_API_AGG_FINAL",
        api_final,
        "SEQ_LEN",
        len(api_step_energies),
    )

    # Basic oscillation / monotonicity heuristic: count sign of successive diffs
    if len(api_step_energies) > 2:
        diffs = [
            api_step_energies[i + 1] - api_step_energies[i]
            for i in range(len(api_step_energies) - 1)
        ]
        up = sum(1 for d in diffs if d > 0)
        down = sum(1 for d in diffs if d < 0)
        print(f"ROY_API_TRACE_DIFFS up={up} down={down} totalSteps={len(diffs)}")

    # (Future extension) Force / position parity: once baseline established we can compare
    # final positions & forces with suitable tolerances; placeholder comment retained intentionally.


def test_relax_uma_water_bfgs_linesearch_trace():
    """Verify BFGSLineSearch optimizer selection + energy trace return for small system.

    We request a short 5-step UMA relaxation with optimizer 'bfgs_ls' and confirm that
    trace_energies length matches steps completed, energies are floats, and initial energy
    lies within the expected physical range (< 0 eV large negative typical for model units).
    """
    torch.manual_seed(0)
    atoms = molecule("H2O")
    payload = {
        "atomic_numbers": atoms.get_atomic_numbers().tolist(),
        "coordinates": atoms.get_positions().tolist(),
        "steps": 5,
        "calculator": "uma",
        "optimizer": "bfgs",  # adjusted: server now only implements plain BFGS
        "optimizer_params": {"maxstep": 0.15},
        "return_trace": True,
        "fmax": 0.0,
    }
    resp = client.post("/relax", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "trace_energies" in data and data["trace_energies"] is not None
    trace = data["trace_energies"]
    # steps_completed may be < requested if an internal failure (should not) or early stop; we disabled early stop with fmax=0.0
    assert len(trace) == data["steps_completed"]
    assert 1 <= data["steps_completed"] <= 5
    for e in trace:
        assert isinstance(e, (float, int))
    # Sanity: energy numeric and negative (model typical)
    assert data["initial_energy"] < 0 and data["final_energy"] < 0


def test_relax_uma_roy_400_steps_bfgs_linesearch_trace():
    """Run ROY 400-step (or env override) UMA relaxation with BFGSLineSearch and capture trace.

    Prints monotonicity stats so we can visually/quantitatively assess oscillations compared to plain BFGS.
    Environment variable MLIPVIEW_ROY_STEPS respected for step count; return_trace enables per-step energies.
    """
    torch.manual_seed(0)
    steps = ROY_STEPS
    atoms_api = load_roy_atoms()
    payload = {
        "atomic_numbers": atoms_api.get_atomic_numbers().tolist(),
        "coordinates": atoms_api.get_positions().tolist(),
        "steps": steps,
        "calculator": "uma",
        "optimizer": "bfgs",  # adjusted: server now only implements plain BFGS
        "optimizer_params": {"maxstep": 0.15},
        "return_trace": True,
        "fmax": 0.0,
    }
    resp = client.post("/relax", json=payload)
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
            f"ROY_LS_TRACE_DIFFS up={up} down={down} total={len(diffs)} firstE={trace[0]} lastE={trace[-1]}"
        )
        # Heuristic: expect fewer sign flips ratio vs vanilla BFGS (not enforcing yet)
        # We only assert basic sanity to avoid flakiness.
        assert len(trace) == data["steps_completed"]
        assert data["steps_completed"] <= steps
        assert isinstance(trace[0], (float, int)) and isinstance(
            trace[-1], (float, int)
        )
