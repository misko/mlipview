from __future__ import annotations

import asyncio
import os
from typing import Generator

import httpx
import pytest
from ray import serve

from fairchem_local_server2.ws_app import deploy
from fairchem_local_server2.model_runtime import (
    UMA_DEPLOYMENT_NAME,
    install_predict_handle,
)


def _cuda_available() -> bool:
    try:
        import torch  # type: ignore

        return bool(torch.cuda.is_available())
    except Exception:
        return False


@pytest.fixture(scope="session")
def ws_base_url() -> Generator[str, None, None]:
    """Start UMA-backed WS server once per session on GPU and yield base URL.

    Fails early if CUDA is not available or UMA device isn't CUDA.
    """
    # Force UMA on GPU only
    assert _cuda_available(), "CUDA is required for UMA tests"

    # Optional: allow overriding HTTP port via env for CI
    base = os.environ.get("WS_BASE", "http://127.0.0.1:8000")
    ws_url = base.replace("http://", "ws://") + "/ws"

    # Deploy UMA with one GPU and minimal CPU/HTTP replicas
    deploy(ngpus=1, ncpus=1, nhttp=1)

    # Poll health endpoint until ready
    health_url = base + "/serve/health"
    stats_url = base + "/uma/stats"
    for _ in range(60):
        try:
            r = httpx.get(health_url, timeout=2.0)
            if r.status_code == 200:
                j = r.json()
                # Model runtime reports global device, should be 'cuda'
                assert (
                    j.get("cuda_available") is True
                ), "health indicates no CUDA availability"
                assert (
                    j.get("device") or ""
                ).lower() == "cuda", f"runtime device is {j.get('device')} not cuda"
                break
        except Exception:
            pass
        import time

        time.sleep(0.5)

    # Query UMA stats from the replica to assert device
    r2 = httpx.get(stats_url, timeout=5.0)
    assert r2.status_code == 200, f"UMA stats failed: {r2.status_code}"
    s = r2.json()
    assert s.get("status") == "ok", f"UMA stats error: {s}"
    stats = s.get("stats") or {}
    dev = (stats.get("device") or "").lower()
    assert dev == "cuda", f"UMA device must be CUDA, got {dev}"

    # Retrieve UMA deployment handle for this app and install into test proc
    from ray import serve as rserve  # lazy import to ensure Serve is up

    # get handle for the UMA deployment within the 'ws_app' application
    uma_h = rserve.get_deployment_handle(UMA_DEPLOYMENT_NAME, app_name="ws_app")
    install_predict_handle(uma_h)

    try:
        yield ws_url
    finally:
        try:
            serve.shutdown()
        except Exception:
            pass


@pytest.fixture(scope="session")
def uma_handle(ws_base_url: str):  # depends on server being up
    from ray import serve as rserve

    return rserve.get_deployment_handle(UMA_DEPLOYMENT_NAME, app_name="ws_app")
