from __future__ import annotations

try:
    import ray  # type: ignore
except Exception:  # pragma: no cover
    ray = None

import numpy as np
import pytest

from fairchem_local_server2.worker_pool import ASEWorker
from fairchem_local_server2.ws_app import deploy


def _uma_available() -> bool:
    try:
        import torch  # type: ignore

        return bool(torch.cuda.is_available())
    except Exception:
        return False


from fairchem_local_server2.atoms_utils import build_atoms
from fairchem_local_server2.models import RelaxCalculatorName
from fairchem_local_server2.services import _md_run


@pytest.mark.parametrize("calculator", ["uma"])
def test_md_single_step_parity_0k_worker_vs_direct(calculator: str):
    if ray is None:
        pytest.skip("ray not installed")
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True)

    if calculator == "uma" and not _uma_available():
        pytest.skip("UMA/GPU not available")

    # Simple dimer with zero initial velocities
    Z = [6, 6]
    R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
    V0 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    if calculator == "uma":
        deploy(ngpus=1, ncpus=1, nhttp=1)
        try:
            from ray import serve as _serve

            from fairchem_local_server2.model_runtime import (
                UMA_DEPLOYMENT_NAME,
                install_predict_handle,
            )

            # Ray Serve API compatibility: retrieve handle by name
            app = _serve.get_app_handle("ws_app")
            # The UMA deployment is part of this app; get a deployment handle
            h = _serve.get_deployment_handle(UMA_DEPLOYMENT_NAME, app_name="ws_app")
            install_predict_handle(h)

            atoms = build_atoms(Z, R, cell=None)
            direct = _md_run(
                atoms,
                steps=5,
                temperature=0.0,
                timestep_fs=1.0,
                friction=0.02,
                calculator=RelaxCalculatorName.uma,
                return_trajectory=False,
                precomputed=None,
                velocities_in=V0,
            ).dict()
            worker = ASEWorker.remote(h)
            fut = worker.run_md.remote(
                atomic_numbers=Z,
                positions=R,
                velocities=V0,
                cell=None,
                steps=5,
                temperature=0.0,
                timestep_fs=1.0,
                friction=0.02,
                calculator="uma",
            )
            result = ray.get(fut)
        finally:
            try:
                from ray import serve as _serve

                _serve.shutdown()
            except Exception:
                pass
    else:
        atoms = build_atoms(Z, R, cell=None)
        direct = _md_run(
            atoms,
            steps=5,
            temperature=0.0,
            timestep_fs=1.0,
            friction=0.02,
            calculator=RelaxCalculatorName.lj,
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

    atol_pos = 1e-12 if calculator == "lj" else 1e-5
    atol_vel = 1e-12 if calculator == "lj" else 1e-5
    atol_for = 1e-12 if calculator == "lj" else 1e-4
    atol_E = 1e-12 if calculator == "lj" else 1e-6

    np.testing.assert_allclose(
        np.array(direct["positions"]),
        np.array(result["positions"]),
        rtol=0,
        atol=atol_pos,
    )
    np.testing.assert_allclose(
        np.array(direct["velocities"]),
        np.array(result["velocities"]),
        rtol=0,
        atol=atol_vel,
    )
    np.testing.assert_allclose(
        np.array(direct["forces"]),
        np.array(result["forces"]),
        rtol=0,
        atol=atol_for,
    )
    assert (
        pytest.approx(direct["final_energy"], rel=0, abs=atol_E)
        == result["final_energy"]
    )
