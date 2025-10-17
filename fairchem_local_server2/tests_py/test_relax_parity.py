from __future__ import annotations

try:
    import ray  # type: ignore
except Exception:  # pragma: no cover
    ray = None

import numpy as np
import pytest

from fairchem_local_server.atoms_utils import build_atoms
from fairchem_local_server.services import _relax_run
from fairchem_local_server.models import RelaxCalculatorName
from fairchem_local_server2.worker_pool import ASEWorker


def test_relax_single_step_parity_cpu_vs_worker():
    if ray is None:
        pytest.skip("ray not installed")
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True)

    # Simple dimer
    Z = [6, 6]
    R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]

    atoms = build_atoms(Z, R, cell=None)
    calc = RelaxCalculatorName.lj

    direct = _relax_run(
        atoms,
        steps=1,
        calculator=calc,
        max_step=0.2,
        return_trace=False,
        precomputed=None,
    ).dict()

    worker = ASEWorker.remote()
    fut = worker.run_relax.remote(
        atomic_numbers=Z,
        positions=R,
        cell=None,
        steps=1,
        fmax=0.05,
        max_step=0.2,
        calculator="lj",
    )
    result = ray.get(fut)

    np.testing.assert_allclose(
        np.array(direct["positions"]),
        np.array(result["positions"]),
        rtol=0,
        atol=1e-6,
    )
