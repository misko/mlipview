import json
import math
import os
import sys
import time

import httpx
import pytest
import ray
from ase import Atoms
from ase.io import read

# Ensure project root on sys.path (two levels up from this file)
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
    # Also propagate PYTHONPATH so Ray worker processes (separate Python interpreters)
    # can import the local package when unpickling deployment definitions.
    os.environ["PYTHONPATH"] = (
        PROJECT_ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")
    )

from fairchem_local_server.server import MDIn, RelaxIn, SimpleIn

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "..", "public", "molecules"
)

MOLECULE_FILES = {
    "water": "water.xyz",
    "benzene": "benzene.xyz",
    "roy": "roy.xyz",
}


@pytest.fixture(scope="session")  # molecules
def molecules():
    out = {}
    for name, fn in MOLECULE_FILES.items():
        path = os.path.abspath(os.path.join(DATA_DIR, fn))
        atoms = read(path)
        out[name] = atoms
    return out


@pytest.fixture(scope="session")  # calculator
def uma_calculator():
    # Acquire calculator from cache for a typical composition (water)
    read(os.path.join(DATA_DIR, "water.xyz"))  # warm IO path
    # Single global calculator now
    from fairchem_local_server import server as _server

    return _server._CALCULATOR


@pytest.fixture(autouse=True)
def debug_device_capability():
    """Print device diagnostics for each test to confirm CUDA vs CPU usage."""
    import torch as _torch

    from fairchem_local_server import server as _server

    dev = "cuda" if _torch.cuda.is_available() else "cpu"
    print(
        "[TEST-DEVICE] cuda=%s UMA_DEVICE=%s srv=%s torch=%s"
        % (
            _torch.cuda.is_available(),
            os.environ.get("UMA_DEVICE"),
            getattr(_server, "DEVICE", "n/a"),
            dev,
        )
    )
    yield


@pytest.fixture(scope="session")  # ray serve
def ray_serve_app():
    # Start ray serve deployment (imports serve_app.deploy)
    import fairchem_local_server.serve_app as serve_app

    if not ray.is_initialized():
        # Keep PYTHONPATH; skip working_dir to avoid packaging venv
        # Exclude common large dirs explicitly if Ray still attempts packaging.
        ray.init(
            include_dashboard=False,
            ignore_reinit_error=True,
            runtime_env={
                # Avoid uploading virtualenv & large artifacts.
                "excludes": [
                    "mlipview_venv",
                    "**/mlipview_venv/**",
                    "**/.git/**",
                    "**/__pycache__/**",
                ],
                "env_vars": {"PYTHONPATH": os.environ.get("PYTHONPATH", "")},
            },
        )
    serve_app.deploy()
    # Ray Serve default HTTP server binds to 8000 unless configured otherwise.
    base_url = "http://127.0.0.1:8000"
    # quick readiness probe
    import time as _t

    import requests

    for _ in range(30):
        try:
            r = requests.get(base_url + "/serve/relax")  # 405 means up
            if r.status_code in (404, 405):
                break
        except Exception:
            _t.sleep(0.2)
    yield base_url
    ray.shutdown()
