from __future__ import annotations

try:
    import ray  # type: ignore
except Exception:  # pragma: no cover
    ray = None

import numpy as np
import pytest

from fairchem_local_server2.worker_pool import ASEWorker
from fairchem_local_server2.atoms_utils import build_atoms
from fairchem_local_server2.models import RelaxCalculatorName
from fairchem_local_server2.services import _relax_run


def _uma_available() -> bool:
    try:
        import torch  # type: ignore

        return bool(torch.cuda.is_available())
    except Exception:
        return False


@pytest.mark.parametrize("calculator", ["uma"])
def test_relax_single_step_parity_cpu_vs_worker(calculator: str, uma_handle):
    if ray is None:
        pytest.skip("ray not installed")
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True)

    if calculator == "uma" and not _uma_available():
        pytest.skip("UMA/GPU not available")

    # Simple dimer
    Z = [6, 6]
    R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]

    if calculator == "uma":
        # UMA server is already running via session fixture; use handle
        from fairchem_local_server2.model_runtime import install_predict_handle

        h = uma_handle
        install_predict_handle(h)

        atoms = build_atoms(Z, R, cell=None)
        direct = _relax_run(
            atoms,
            steps=1,
            calculator=RelaxCalculatorName.uma,
            max_step=0.2,
            return_trace=False,
            precomputed=None,
        ).dict()

        worker = ASEWorker.remote(h)
        fut = worker.run_relax.remote(
            atomic_numbers=Z,
            positions=R,
            cell=None,
            steps=1,
            fmax=0.05,
            max_step=0.2,
            calculator="uma",
        )
        result = ray.get(fut)
    else:
        atoms = build_atoms(Z, R, cell=None)
        direct = _relax_run(
            atoms,
            steps=1,
            calculator=RelaxCalculatorName.lj,
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
