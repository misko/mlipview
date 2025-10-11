import os

from ase.io.jsonio import encode
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from ray import serve

from .model_runtime import (
    MODEL_NAME,
    TASK_NAME,
    UMA_DEPLOYMENT_NAME,
    _PredictDeploy,
    health_snapshot,
    install_predict_handle,
)
from .models import (
    MDFromCacheIn,
    MDIn,
    MDResult,
    RelaxCalculatorName,
    RelaxFromCacheIn,
    RelaxIn,
    RelaxResult,
    SimpleFromCacheIn,
    SimpleIn,
)
from .services import (
    md_step,
    md_step_from_cache,
    relax,
    relax_from_cache,
    simple_calculate,
    simple_calculate_from_cache,
)

app = FastAPI(title="UMA-small ASE HTTP server (slim)", debug=True)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
async def _startup_install_uma():
    """Ensure UMA Ray Serve deployment is running and install its handle.

    This allows the plain FastAPI server to use the same batched GPU-backed
    predictor as the Serve ingress variant, scaling replicas to the number of
    GPUs (or UMA_NUM_GPUS env var if provided).
    """
    # Choose replica count: env override -> auto -> 1
    ngpus_env = os.environ.get("UMA_NUM_GPUS") or os.environ.get("NGPUS")
    try:
        ngpus = int(ngpus_env) if ngpus_env is not None else None
    except Exception:
        ngpus = None
    if ngpus is None:
        try:
            import torch  # type: ignore

            ngpus = torch.cuda.device_count() if torch.cuda.is_available() else 1
        except Exception:
            ngpus = 1
    if ngpus <= 0:
        ngpus = 1

    # Deploy or update UMA predictor with given replica count
    _PredictDeploy.options(
        name=UMA_DEPLOYMENT_NAME,
        num_replicas=int(ngpus),
        ray_actor_options={"num_gpus": 1},
    ).deploy(MODEL_NAME, TASK_NAME)

    # Acquire a handle from Serve and install it for service layer usage
    handle = serve.get_deployment(UMA_DEPLOYMENT_NAME).get_handle()
    install_predict_handle(handle)


@app.get("/health")
def health():
    snap = health_snapshot()
    snap["status"] = "ok"
    return snap


@app.post("/simple_calculate")
def simple_calculate_ep(inp: SimpleIn):
    return simple_calculate(inp)


@app.post("/simple_calculate_from_cache")
def simple_calculate_from_cache_ep(inp: SimpleFromCacheIn):
    return simple_calculate_from_cache(inp)


@app.post("/relax", response_model=RelaxResult)
def relax_ep(inp: RelaxIn):
    return relax(inp)


@app.post("/relax_from_cache", response_model=RelaxResult)
def relax_from_cache_ep(inp: RelaxFromCacheIn):
    return relax_from_cache(inp)


@app.post("/md", response_model=MDResult)
def md_ep(inp: MDIn):
    return md_step(inp)


@app.post("/md_from_cache", response_model=MDResult)
def md_from_cache_ep(inp: MDFromCacheIn):
    return md_step_from_cache(inp)


__all__ = [
    "app",
    "encode",
    "SimpleIn",
    "RelaxIn",
    "MDIn",
    "RelaxResult",
    "MDResult",
    "RelaxCalculatorName",
]
