import math
import time

import httpx
import pytest
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS

from fairchem_local_server import server as _server
from fairchem_local_server.server import MDIn, RelaxCalculatorName, RelaxIn, SimpleIn

RELAX_STEPS = 50
MD_STEPS = 50
TEMP_K = 298.0


@pytest.mark.parametrize("mol_name", ["water", "benzene", "roy"])
@pytest.mark.asyncio
async def test_single_point_relax_md_equivalence(
    molecules, local_uma_calculator, ray_serve_app, mol_name
):
    atoms = molecules[mol_name].copy()
    # Ensure pure python ints for JSON serialization
    Z = [int(z) for z in atoms.get_atomic_numbers()]
    atoms.calc = local_uma_calculator
    # Baseline initial energy
    e0 = float(atoms.get_potential_energy())
    # Relax baseline
    opt = BFGS(atoms, logfile=None)
    for _ in range(RELAX_STEPS):
        opt.step()
    e_relax = float(atoms.get_potential_energy())

    # MD baseline (reinitialize velocities)
    MaxwellBoltzmannDistribution(atoms, temperature_K=TEMP_K)
    dyn = Langevin(
        atoms, 1.0 * units.fs, temperature_K=TEMP_K, friction=0.02, logfile=None
    )
    for _ in range(MD_STEPS):
        dyn.run(1)
    e_md = float(atoms.get_potential_energy())

    async with httpx.AsyncClient(timeout=30.0) as client:
        # Relax via serve (batched) using /serve/relax
        relax_payload = {
            "atomic_numbers": Z,
            "coordinates": molecules[mol_name].get_positions().tolist(),
            "steps": RELAX_STEPS,
            "calculator": "uma",
            "return_trace": False,
        }
        r = await client.post(f"{ray_serve_app}/serve/relax", json=relax_payload)
        assert r.status_code == 200, r.text
        rj = r.json()
        serve_e0 = float(rj["initial_energy"])
        serve_efin = float(rj["final_energy"])

        assert math.isclose(e0, serve_e0, rel_tol=1e-6, abs_tol=1e-4)
        assert math.isclose(e_relax, serve_efin, rel_tol=1e-5, abs_tol=1e-4)

        # MD via serve
        md_payload = {
            "atomic_numbers": [int(z) for z in atoms.get_atomic_numbers()],
            "coordinates": molecules[mol_name].get_positions().tolist(),
            "steps": MD_STEPS,
            "temperature": TEMP_K,
            "timestep_fs": 1.0,
            "calculator": "uma",
        }
        r2 = await client.post(f"{ray_serve_app}/serve/md", json=md_payload)
        assert r2.status_code == 200, r2.text
        r2j = r2.json()
        serve_e0_md = float(r2j["initial_energy"])
        serve_efin_md = float(r2j["final_energy"])
        assert math.isclose(e0, serve_e0_md, rel_tol=1e-6, abs_tol=1e-4)
        # MD energy may drift; allow looser tolerance
        assert math.isclose(e_md, serve_efin_md, rel_tol=5e-3, abs_tol=5e-3)
        # Blow-up check: positions shouldn't have NaN or huge displacement
        for pos in r2j["positions"]:
            assert all(math.isfinite(p) for p in pos)
