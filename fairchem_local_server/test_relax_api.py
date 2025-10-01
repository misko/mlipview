import math
import torch
from fastapi.testclient import TestClient
from ase.build import molecule
from ase.optimize import BFGS

from fairchem_local_server import server  # adjust to import via module path if package, else fallback
from fairchem_local_server.server import app, get_cached_unit_and_calculator
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


def run_direct_bfgs(atoms, steps: int):
    energies = [float(atoms.get_potential_energy())]
    opt = BFGS(atoms, logfile=None)
    for _ in range(steps):
        opt.step()
        energies.append(float(atoms.get_potential_energy()))
    return (
        energies,
        atoms.get_positions().tolist(),
        atoms.get_forces().tolist(),
    )


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
    assert math.isclose(
        data["final_energy"], energies_cmp[-1], rel_tol=0, abs_tol=1e-5
    )

    # Position comparison
    for a, b in zip(data["positions"], pos_cmp):
        for x, y in zip(a, b):
            assert abs(x - y) < 1e-5


def test_relax_uma_water():
    # UMA model may introduce stochasticity; fix torch seed
    torch.manual_seed(0)
    atoms = molecule("H2O")
    # direct UMA relaxation
    _pu, calc = get_cached_unit_and_calculator(
        atoms.get_atomic_numbers().tolist()
    )
    atoms.calc = calc
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
    _pu2, calc2 = get_cached_unit_and_calculator(
        atoms_init.get_atomic_numbers().tolist()
    )
    atoms_init.calc = calc2
    e0_direct = float(atoms_init.get_potential_energy())
    # Allow slightly looser tolerance (ML model nondeterminism / device differences)
    assert math.isclose(
        data["initial_energy"], e0_direct, rel_tol=0, abs_tol=5e-4
    )

    # final energy tolerance looser due to possible tiny nondeterminism
    assert math.isclose(
        data["final_energy"],
        energies_direct[-1],
        rel_tol=1e-4,
        abs_tol=1e-3,
    )

    # positions roughly close
    for a, b in zip(data["positions"], pos_direct):
        for x, y in zip(a, b):
            assert abs(x - y) < 5e-4
