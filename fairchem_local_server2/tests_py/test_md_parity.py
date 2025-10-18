from __future__ import annotations

try:
    import ray  # type: ignore
except Exception:  # pragma: no cover
    ray = None

import numpy as np
import pytest

from fairchem_local_server.atoms_utils import build_atoms
from fairchem_local_server.services import _md_run
from fairchem_local_server.models import RelaxCalculatorName
from fairchem_local_server2.worker_pool import ASEWorker


def test_md_single_step_parity_0k_worker_vs_direct():
    if ray is None:
        pytest.skip("ray not installed")
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True)

    # Simple dimer with zero initial velocities
    Z = [6, 6]
    R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
    V0 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    atoms = build_atoms(Z, R, cell=None)
    calc = RelaxCalculatorName.lj

    direct = _md_run(
        atoms,
        steps=5,
        temperature=0.0,
        timestep_fs=1.0,
        friction=0.02,
        calculator=calc,
        return_trajectory=False,
        precomputed=None,
        velocities_in=V0,
    ).dict()

    worker = ASEWorker.remote()
    fut = worker.run_md.remote(
        atomic_numbers=Z,
        positions=R,
        velocities=V0,
        cell=None,
        steps=5,
        temperature=0.0,
        timestep_fs=1.0,
        friction=0.02,
        calculator="lj",
    )
    result = ray.get(fut)

    np.testing.assert_allclose(
        np.array(direct["positions"]),
        np.array(result["positions"]),
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.array(direct["velocities"]),
        np.array(result["velocities"]),
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.array(direct["forces"]),
        np.array(result["forces"]),
        rtol=0,
        atol=1e-12,
    )
    assert pytest.approx(direct["final_energy"], rel=0, abs=1e-12) == result[
        "final_energy"
    ]
