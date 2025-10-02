"""Modern Ray Serve application with batching and FastAPI ingress.

Pattern:
    * Three compute deployments: SimpleCalc, RelaxCalc, MDCalc.
    * Each defines a batched method via @serve.batch and a single-item __call__ wrapper.
    * An Ingress deployment (with @serve.ingress(FastAPI)) exposes HTTP routes under /serve/*.
    * serve.build + serve.run deploy the graph in one call (recommended pattern Ray Serve >=2.4).

Backward compatibility: routes keep the /serve/ prefix.
"""

from __future__ import annotations

import asyncio
import os
from typing import Any, Dict, List

import torch  # local import so health can report cuda availability
from fastapi import FastAPI, HTTPException
from ray import serve

from .server import _MODEL_LOADED  # to report model load status
from .server import DEVICE  # expose selected device from server module
from .server import MDIn, RelaxIn, SimpleIn
from .server import md_step as md_fn
from .server import relax as relax_fn
from .server import simple_calculate as simple_calculate_fn

MAX_BATCH = 8
WAIT_S = 0.02


_ingress_app = FastAPI(title="UMA Ray Serve batched API")

_runtime_env = {}
if os.environ.get("CUDA_VISIBLE_DEVICES"):
    # Propagate explicit restriction if user set it
    _runtime_env["env_vars"] = {"CUDA_VISIBLE_DEVICES": os.environ["CUDA_VISIBLE_DEVICES"]}

@serve.deployment(ray_actor_options={"num_gpus": 0.1, "runtime_env": _runtime_env})
@serve.ingress(_ingress_app)
class Ingress:
    def __init__(self):
        # Log device context from inside the Ray worker (after actor init) so we see CUDA_VISIBLE_DEVICES
        try:
            print(
                f"[serve:Ingress:init] DEVICE={DEVICE} cuda_available={torch.cuda.is_available()} "
                f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')}"
            )
        except Exception as e:  # pragma: no cover
            print("[serve:Ingress:init] device log failed", e)

    @_ingress_app.get("/serve/health")
    async def _health(self):  # type: ignore
        # Report server-selected device plus runtime CUDA visibility and model load status.
        return {
            "status": "ok",
            "batched": True,
            "device": DEVICE,
            "cuda_available": torch.cuda.is_available(),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "model_loaded": bool(_MODEL_LOADED),
        }


    @_ingress_app.post("/serve/simple")
    async def _simple_one(self, body: Dict[str, Any]):  # type: ignore
        try:
            sin = SimpleIn(**body)
            return simple_calculate_fn(sin)
        except HTTPException as he:
            return {"error": he.detail}
        except Exception as e:  # pragma: no cover
            return {"error": str(e)}


    @_ingress_app.post("/serve/simple/batch")
    async def _simple_batch(self, body: List[Dict[str, Any]]):  # type: ignore
        out: List[Any] = []
        for b in body:
            out.append(await self._simple_one(b))
        return out

    @_ingress_app.post("/serve/relax")
    async def _relax_one(self, body: Dict[str, Any]):  # type: ignore
        try:
            rin = RelaxIn(**body)
            return relax_fn(rin).dict()
        except HTTPException as he:
            return {"error": he.detail}
        except Exception as e:
            return {"error": str(e)}


    @_ingress_app.post("/serve/relax/batch")
    async def _relax_batch(self, body: List[Dict[str, Any]]):  # type: ignore
        out: List[Any] = []
        for b in body:
            out.append(await self._relax_one(b))
        return out

    @_ingress_app.post("/serve/md")
    async def _md_one(self, body: Dict[str, Any]):  # type: ignore
        try:
            minp = MDIn(**body)
            return md_fn(minp).dict()
        except HTTPException as he:
            return {"error": he.detail}
        except Exception as e:
            return {"error": str(e)}


    @_ingress_app.post("/serve/md/batch")
    async def _md_batch(self, body: List[Dict[str, Any]]):  # type: ignore
        out: List[Any] = []
        for b in body:
            out.append(await self._md_one(b))
        return out


def deploy():
    if DEVICE != "cuda":  # sanity guard; server.py now enforces cuda-only
        raise RuntimeError(f"Ray Serve deploy aborted: DEVICE={DEVICE} is not cuda")
    import ray

    # Initialize Ray if not already (local mode). Do NOT limit GPUs here; let Ray detect.
    if not ray.is_initialized():
        ray.init(include_dashboard=False, ignore_reinit_error=True)
    try:
        cluster = ray.cluster_resources()
        avail = ray.available_resources()
        print(f"[serve:deploy] cluster_resources={cluster}")
        print(f"[serve:deploy] available_resources={avail}")
        g_total = cluster.get("GPU", 0.0)
        if g_total <= 0:
            raise RuntimeError(
                f"No GPU resources registered with Ray (cluster_resources GPU={g_total}); cannot deploy GPU-only service"
            )
    except Exception as e:
        print("[serve:deploy] resource inspection failed", e)
        raise
    serve.start(detached=False, http_options={"host": "127.0.0.1", "port": 8000})
    app_dag = Ingress.bind()
    # route_prefix "/" exposes /serve/* endpoints defined on the FastAPI app
    serve.run(app_dag, route_prefix="/")
    return app_dag


# Allow `python -m fairchem_local_server.serve_app` manual startup
if __name__ == "__main__":  # pragma: no cover
    deploy()
    import time
    print("Serve running at http://127.0.0.1:8000/serve/health")
    while True:
        time.sleep(60)
