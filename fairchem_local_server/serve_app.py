"""Ray Serve ingress that reuses the pure function API implementations.

This file wires HTTP routes -> api.* functions and reports health.
Model inference batching is handled by the internal deployment created
in the batched_predict module.
"""

from __future__ import annotations

import os

import torch
from fastapi import FastAPI
from ray import serve

from .api import md_step, relax, simple_calculate
from .model_runtime import DEVICE, ensure_model_loaded, health_snapshot
from .models import MDIn, RelaxIn, SimpleIn

_ingress_app = FastAPI(title="UMA Ray Serve API (slim)")

_runtime_env = {}
if os.environ.get("CUDA_VISIBLE_DEVICES"):
    # Propagate explicit restriction if user set it
    _runtime_env["env_vars"] = {
        "CUDA_VISIBLE_DEVICES": os.environ["CUDA_VISIBLE_DEVICES"]
    }


@serve.deployment(ray_actor_options={"num_gpus": 0.1, "runtime_env": _runtime_env})
@serve.ingress(_ingress_app)
class Ingress:
    def __init__(self):
        # Log device context after actor init (visibility of CUDA env)
        try:
            cuda_ok = torch.cuda.is_available()
            cvd = os.environ.get("CUDA_VISIBLE_DEVICES")
            print(f"[serve:init] DEV={DEVICE} cuda={cuda_ok} CVD={cvd}")
        except Exception as e:  # pragma: no cover
            print("[serve:Ingress:init] device log failed", e)

    @_ingress_app.get("/serve/health")
    async def _health(self):  # type: ignore
        snap = health_snapshot()
        snap["status"] = "ok"
        return snap

    @_ingress_app.post("/serve/simple")
    async def _simple(self, body: dict):  # type: ignore
        sin = SimpleIn(**body)
        return simple_calculate(sin)

    @_ingress_app.post("/serve/relax")
    async def _relax(self, body: dict):  # type: ignore
        rin = RelaxIn(**body)
        return relax(rin).dict()

    @_ingress_app.post("/serve/relax/batch")
    async def _relax_batch(self, body: list):  # type: ignore
        out = []
        for item in body:
            try:
                rin = RelaxIn(**item)
                out.append(relax(rin).dict())
            except Exception as e:  # pragma: no cover
                out.append({"error": str(e)})
        return out

    @_ingress_app.post("/serve/md")
    async def _md(self, body: dict):  # type: ignore
        minp = MDIn(**body)
        return md_step(minp).dict()


def deploy():
    if DEVICE != "cuda":
        raise RuntimeError(f"Ray Serve deploy aborted: DEVICE={DEVICE} not cuda")

    import ray

    if not ray.is_initialized():
        ray.init(include_dashboard=False, ignore_reinit_error=True)

    # Pre-warm model so startup errors surface early.
    ensure_model_loaded()

    # Bring up Serve HTTP head on localhost:8000
    serve.start(detached=False, http_options={"host": "127.0.0.1", "port": 8000})

    dag = Ingress.bind()
    serve.run(dag, route_prefix="/")
    return dag


# Allow `python -m fairchem_local_server.serve_app` manual startup
if __name__ == "__main__":  # pragma: no cover
    deploy()
    import time

    print("Serve running at http://127.0.0.1:8000/serve/health")
    while True:
        time.sleep(60)
