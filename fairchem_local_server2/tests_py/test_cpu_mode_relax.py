from __future__ import annotations

import os

import pytest
import ray
from ray import serve

from fairchem_local_server2.model_runtime import (
    UMA_DEPLOYMENT_NAME,
    install_predict_handle,
)
from fairchem_local_server2.ws_app import deploy
from fairchem_local_server2.worker_pool import WorkerPool


@pytest.mark.timeout(60)
def test_cpu_mode_relaxation(monkeypatch):
    """Ensure MLIPVIEW_FORCE_CPU allows relaxing a dimer without CUDA."""

    monkeypatch.setenv("MLIPVIEW_FORCE_CPU", "1")
    # Use a non-default HTTP port to avoid clashes with other tests
    deploy(ngpus=None, ncpus=1, nhttp=1, http_port=8051)
    try:
        handle = serve.get_deployment_handle(UMA_DEPLOYMENT_NAME, app_name="ws_app")
        install_predict_handle(handle)

        pool = WorkerPool(size=1, uma_handle=handle)
        worker = pool.any()
        fut = worker.run_relax.remote(
            atomic_numbers=[6, 6],
            positions=[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            cell=None,
            steps=10,
            fmax=0.05,
            max_step=0.2,
            calculator="uma",
            precomputed=None,
        )
        result = ray.get(fut)
        assert int(result.get("steps_completed", 0)) >= 10
    finally:
        serve.shutdown()
        monkeypatch.delenv("MLIPVIEW_FORCE_CPU", raising=False)
