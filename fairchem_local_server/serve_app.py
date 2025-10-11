import logging
import os
import time
from typing import Optional

from fastapi import FastAPI
from ray import serve
from ray.serve.schema import LoggingConfig

from .model_runtime import (
    MODEL_NAME,
    TASK_NAME,
    UMA_DEPLOYMENT_NAME,
    _PredictDeploy,
    health_snapshot,
    install_predict_handle,
)
from .models import MDIn, RelaxIn, SimpleIn
from .services import md_step, relax, simple_calculate

app = FastAPI(title="UMA Serve API", debug=False)


@serve.deployment(
    ray_actor_options={"num_gpus": 0},
    logging_config={"enable_access_log": False},
)
@serve.ingress(app)
class Ingress:
    def __init__(self, predict_handle):
        # Inject the UMA handle; builds FAIRChemCalculator.
        logging.getLogger("ray.serve").setLevel(logging.ERROR)
        install_predict_handle(predict_handle)

    @app.get("/serve/health")
    def health(self):
        snap = health_snapshot()
        snap["status"] = "ok"
        return snap

    # SYNC endpoints so Starlette runs them in a threadpool,
    # which makes the blocking FAIRChem/ASE code safe.
    @app.post("/serve/simple")
    def simple_ep(self, inp: SimpleIn):
        return simple_calculate(inp)

    @app.post("/serve/relax")
    def relax_ep(self, inp: RelaxIn):
        return relax(inp).dict()

    @app.post("/serve/md")
    def md_ep(self, inp: MDIn):
        return md_step(inp).dict()


def _detect_default_ngpus() -> int:
    """Choose a sensible default replica count based on environment.

    Priority:
      1) UMA_NUM_GPUS env var if set and valid (>0)
      2) torch.cuda.device_count() if torch is available and >0
      3) fallback to 1
    """
    # 1) Explicit env override
    try:
        if "UMA_NUM_GPUS" in os.environ:
            v = int(os.environ.get("UMA_NUM_GPUS", "1") or "1")
            return max(1, v)
    except Exception:
        pass
    # 2) Probe torch for GPU count
    try:
        import torch  # type: ignore

        n = torch.cuda.device_count() if torch.cuda.is_available() else 0
        if n and n > 0:
            return int(n)
    except Exception:
        pass
    # 3) Fallback
    return 1


def deploy(ngpus: Optional[int] = None):
    # Ray + Serve will be started implicitly by serve.run if needed.
    # Each UMA predictor replica consumes 1 GPU; scale replicas == number of
    # GPUs while keeping the ingress on CPU-only.
    replica_count = int(ngpus) if ngpus is not None else _detect_default_ngpus()
    if replica_count <= 0:
        replica_count = 1

    uma = _PredictDeploy.options(
        name=UMA_DEPLOYMENT_NAME,
        num_replicas=replica_count,
        ray_actor_options={"num_gpus": 1},
    ).bind(MODEL_NAME, TASK_NAME)
    dag = Ingress.bind(uma)
    serve.run(
        dag,
        name="http_app",
        route_prefix="/",
        logging_config=LoggingConfig(
            enable_access_log=False,  # hide frequent access logs
            log_level="ERROR",  # optional: raise log level for Serve internals
        ),
    )
    return dag


if __name__ == "__main__":
    # Lightweight CLI to launch with a specific number of GPU replicas
    import argparse

    parser = argparse.ArgumentParser(description="Start UMA Serve app")
    parser.add_argument(
        "--ngpus",
        type=int,
        default=None,
        help=(
            "Number of GPU-backed UMA replicas (defaults to UMA_NUM_GPUS or "
            "auto-detect)"
        ),
    )
    args = parser.parse_args()

    deploy(ngpus=args.ngpus)
    while True:
        time.sleep(60)
