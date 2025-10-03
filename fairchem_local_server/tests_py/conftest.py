import os
import sys
import time

import pytest
import requests
import torch
from ase.io import read
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.units.mlip_unit import InferenceSettings
from ray import serve

# --- Project path / PYTHONPATH ------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
    # Ensure worker processes can import your package
    os.environ["PYTHONPATH"] = (
        PROJECT_ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")
    )

# Prefer importing models directly from your package
from fairchem_local_server.model_runtime import MODEL_NAME, TASK_NAME
from fairchem_local_server.models import (  # noqa: F401  (used by tests via import)
    MDIn,
    RelaxIn,
    SimpleIn,
)

# --- Molecule fixtures --------------------------------------------------------
DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "..", "public", "molecules"
)
MOLECULE_FILES = {
    "water": "water.xyz",
    "benzene": "benzene.xyz",
    "roy": "roy.xyz",
}


@pytest.fixture(scope="session")
def local_uma_calculator():
    assert torch.cuda.is_available(), "UMA calculator fixture requires CUDA"
    return FAIRChemCalculator(
        pretrained_mlip.get_predict_unit(
            MODEL_NAME,
            device="cuda",
        ),
        task_name=TASK_NAME,
    )


@pytest.fixture(scope="session")
def molecules():
    out = {}
    for name, fn in MOLECULE_FILES.items():
        path = os.path.abspath(os.path.join(DATA_DIR, fn))
        out[name] = read(path)
    return out


# --- Molecule fixtures --------------------------------------------------------
DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "..", "public", "molecules"
)
MOLECULE_FILES = {
    "water": "water.xyz",
    "benzene": "benzene.xyz",
    "roy": "roy.xyz",
}


@pytest.fixture(scope="session")
def molecules():
    out = {}
    for name, fn in MOLECULE_FILES.items():
        path = os.path.abspath(os.path.join(DATA_DIR, fn))
        out[name] = read(path)
    return out


# --- Debug print for device per test -----------------------------------------
@pytest.fixture(autouse=True)
def debug_device_capability():
    import torch

    # DEVICE is defined in your model_runtime; this import is cheap and safe.
    try:
        from fairchem_local_server.model_runtime import DEVICE  # type: ignore
    except Exception:
        DEVICE = "unknown"
    print(
        "[TEST-DEVICE] cuda=%s UMA_DEVICE=%s srv=%s torch=%s"
        % (
            torch.cuda.is_available(),
            os.environ.get("UMA_DEVICE"),
            DEVICE,
            "cuda" if torch.cuda.is_available() else "cpu",
        )
    )
    yield


# --- Ray Serve app (HTTP) bootstrap ------------------------------------------
@pytest.fixture(scope="session")
def ray_serve_app():
    """
    Starts the Ray Serve application (HTTP ingress + UMA batched predictor) once
    for the entire test session by calling serve_app.deploy().

    Returns the base URL for tests that hit HTTP endpoints.
    """
    os.environ["UMA_BATCH_MAX"] = "64"
    os.environ["UMA_BATCH_WAIT_S"] = "0.05"  # 50 ms makes batching obvious in tests
    # Ray is started implicitly by serve.run() in deploy(); no explicit ray.init needed.
    from fairchem_local_server import serve_app  # your Serve-owns-HTTP entrypoint

    base_url = "http://127.0.0.1:8000"  # serve_app.deploy() binds to this

    # Deploy (idempotent across repeated calls)
    serve_app.deploy()

    # Probe readiness: /serve/health returns 200 when FastAPI ingress is live.
    deadline = time.time() + 30.0
    last_err = None
    while time.time() < deadline:
        try:
            r = requests.get(base_url + "/serve/health", timeout=1.0)
            if r.status_code == 200:
                break
            last_err = f"status={r.status_code}, text={r.text[:200]}"
        except Exception as e:
            last_err = str(e)
        time.sleep(0.2)
    else:
        pytest.fail(f"Serve app did not become ready: {last_err}")

    yield base_url

    # Teardown: shut down Ray to free the HTTP port for subsequent sessions.
    # (serve.shutdown() is implicit in ray.shutdown())
    import ray

    try:
        ray.shutdown()
    except Exception:
        pass
