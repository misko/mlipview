"""Ray Serve batched inference wrapper for UMA predict unit.

Single-GPU replica with @serve.batch on .predict. The local wrapper exposes
a synchronous .predict(...) that returns ONE result per call while Ray Serve
coalesces concurrent requests into batches under the hood.
"""

from __future__ import annotations

import os
from typing import Any, List, Tuple

import ray
from fairchem.core import pretrained_mlip
from fairchem.core.units.mlip_unit import InferenceSettings
from ray import serve
from ray.serve.handle import DeploymentResponse  # public in 2.49.x

# Tunables: how big and how long we wait to form a batch.
MAX_BATCH = int(os.environ.get("UMA_BATCH_MAX", 16))
WAIT_S = float(os.environ.get("UMA_BATCH_WAIT_S", 0.003))


@serve.deployment(ray_actor_options={"num_gpus": 1})
class _PredictDeploy:  # pragma: no cover (runs remotely)
    def __init__(self, model_name: str, task_name: str):
        print(f"[batched:init] loading model={model_name} task={task_name}")
        self._unit = pretrained_mlip.get_predict_unit(
            model_name,
            device="cuda",
            inference_settings=InferenceSettings(
                tf32=True,
                activation_checkpointing=False,
                merge_mole=False,
                compile=False,
                external_graph_gen=False,
                internal_graph_gen_version=2,
            ),
        )
        self._meta = {
            "dataset_to_tasks": getattr(self._unit, "dataset_to_tasks", {}),
            "device": "cuda",
        }

    # In Ray 2.49.x, the batched function MUST be async.
    @serve.batch(max_batch_size=MAX_BATCH, batch_wait_timeout_s=WAIT_S)
    async def predict(self, payloads: List[Tuple[tuple, dict]]):
        """payloads: [(args, kwargs), ...]; return one output per input in order."""
        # Keep order stable; let exceptions propagate for proper per-request errors.
        out: List[Any] = []
        for args, kwargs in payloads:
            out.append(self._unit.predict(*args, **kwargs))
        return out

    async def get_meta(self) -> dict:  # type: ignore[override]
        return self._meta


class BatchedPredictUnit:
    """Synchronous client wrapper. Each .predict(...) call returns ONE result.

    Ray Serve will batch multiple concurrent calls to _PredictDeploy.predict.
    """

    def __init__(
        self,
        handle,
        model_name: str | None = None,
        task_name: str | None = None,
    ):
        self._handle = handle
        self._model_name = model_name or os.environ.get("UMA_MODEL_NAME", "uma-s-1p1")
        self._task_name = task_name or os.environ.get("UMA_TASK_NAME", "omol")

        # Populate metadata once (best effort). In 2.49.x .remote(...) returns DeploymentResponse.
        try:
            m = self._handle.get_meta.remote()  # type: ignore[attr-defined]
            meta = m.result() if isinstance(m, DeploymentResponse) else ray.get(m)
        except Exception:  # pragma: no cover
            meta = {"dataset_to_tasks": {}, "device": "cuda"}

        self.dataset_to_tasks = meta.get("dataset_to_tasks", {}) or {}
        self.device = meta.get("device", "cuda")

        if not self.dataset_to_tasks:
            # Minimal stub for callers that expect FAIRChem-style task descriptors.
            class _TaskStub:
                def __init__(self, prop: str):
                    self.property = prop

            self.dataset_to_tasks = {
                "omol": [
                    _TaskStub("energy"),
                    _TaskStub("forces"),
                    _TaskStub("stress"),
                ]
            }

    def predict(self, *args, **kwargs):
        """Return SINGLE result for this call.

        We pass exactly one (args, kwargs) payload. Serve batches such payloads
        and maps each caller's result back by index.
        """
        res = self._handle.predict.remote((args, kwargs))  # type: ignore[attr-defined]
        # Ray 2.49.x returns a DeploymentResponse here; older code may return ObjectRef.
        return res.result() if isinstance(res, DeploymentResponse) else ray.get(res)

    def __getattr__(self, item):  # pragma: no cover
        if item in {"dataset_to_tasks", "device"}:
            return self.__dict__[item]
        raise AttributeError(item)


_DEPLOYED = False
_WRAPPER: BatchedPredictUnit | None = None


def get_batched_predict_unit(model_name: str, task_name: str) -> BatchedPredictUnit:
    """Create (once) and return the Ray Serve batched predict unit wrapper."""
    global _DEPLOYED, _WRAPPER
    if _WRAPPER is not None:
        return _WRAPPER

    if not ray.is_initialized():
        ray.init(include_dashboard=False, ignore_reinit_error=True)

    # Start Serve controller (no HTTP server here; your ingress manages HTTP).
    try:
        serve.start()
    except Exception:
        pass  # already started

    if not _DEPLOYED:
        # Public API path for 2.49.2: bind -> run -> get_deployment_handle
        app = _PredictDeploy.options(name="uma_predict").bind(model_name, task_name)
        serve.run(app)  # idempotent by deployment name within the default app

        # Get a *deployment* handle by name (public API).
        # If you later name the app via serve.run(..., name="myapp"), set app_name="myapp".
        handle = serve.get_deployment_handle("uma_predict", app_name="default")

        _WRAPPER = BatchedPredictUnit(handle, model_name, task_name)
        _DEPLOYED = True

    return _WRAPPER


class BatchPredictUnit(BatchedPredictUnit):  # backwards compatibility alias
    pass


__all__ = ["get_batched_predict_unit", "BatchedPredictUnit", "BatchPredictUnit"]
