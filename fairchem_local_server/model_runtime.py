"""Central runtime utilities for UMA model + ASE calculator.

Serve-owns-HTTP model:
- This module defines the UMA batched predictor deployment (no HTTP)
- The Serve ingress app constructs the DAG and injects a handle via
  install_predict_handle(...)
- We cache a FAIRChemCalculator backed by the remote predictor
- No serve.start / serve.run in this file
"""

from __future__ import annotations

import os
from typing import Any, List, Tuple

import ray
import torch
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch
from fairchem.core.units.mlip_unit import InferenceSettings
from ray import serve

# --- Config -----------------------------------------------------------------

DEVICE = "cuda"
MODEL_NAME = os.getenv("UMA_MODEL", "uma-s-1p1")
TASK_NAME = os.getenv("UMA_TASK", "omol")

# Batch tunables
MAX_BATCH = int(os.environ.get("UMA_BATCH_MAX", 16))
WAIT_S = float(os.environ.get("UMA_BATCH_WAIT_S", 0.003))

# Logical deployment name (must match the one you bind in the Serve DAG)
UMA_DEPLOYMENT_NAME = "uma_predict"

# --- State ------------------------------------------------------------------

_PU: BatchedPredictUnit | None = None
_CALC: FAIRChemCalculator | None = None


# --- Ray Serve UMA deployment (no HTTP route) --------------------------------


@serve.deployment(ray_actor_options={"num_gpus": 1})
class _PredictDeploy:  # runs on GPU replica
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

    # Batched inference (async, required by Serve)
    @serve.batch(max_batch_size=MAX_BATCH, batch_wait_timeout_s=WAIT_S)
    async def predict(self, payloads: List[Tuple[tuple, dict]]):
        # Preserve order; return one result per payload.
        out: List[Any] = []
        if len(payloads) == 1:
            return [
                {
                    k: v.detach().cpu()
                    for k, v in self._unit.predict(payloads[0][0][0]).items()
                }
            ]
        else:
            batch = atomicdata_list_to_batch([x[0][0] for x in payloads])  # warmup
            batch.dataset = [x[0] for x in batch.dataset]
            all_outputs = {
                k: v.detach().cpu() for k, v in self._unit.predict(batch).items()
            }
            forces_by_mol = all_outputs["forces"].split(batch["natoms"].tolist())
            # stress_by_mol = all_outputs["stress"].split(batch["num_atoms"])
            out = [
                {"energy": energy, "forces": forces, "stress": stress[0]}
                for energy, forces, stress in zip(
                    all_outputs["energy"].split(1),
                    forces_by_mol,
                    all_outputs["stress"].split(1),
                )
            ]
        return out


class BatchedPredictUnit:
    """Synchronous client wrapper (FAIRChem expects sync .predict)."""

    def __init__(self, handle, dataset_to_tasks: dict | None = None):
        self._handle = handle
        # Provide constant-enough metadata; avoids any remote call here.
        # If you prefer to fetch the real one, do it at Serve ingress __init__
        # and pass it in via 'dataset_to_tasks'.
        self.dataset_to_tasks = dataset_to_tasks or {
            "omol": [
                {"property": "energy"},
                {"property": "forces"},
                {"property": "stress"},
            ]
        }
        self.device = "cuda"

    def predict(self, *args, **kwargs):
        # print("RUNNING INFERENCE", flush=True)
        # NOTE: This is called from sync code (e.g., FastAPI sync handler running
        # in a thread pool). DeploymentResponse.result() is safe in that context.
        resp = self._handle.predict.remote((args, kwargs))  # type: ignore[attr-defined]
        r = resp.result()
        assert r.keys() == {"energy", "forces", "stress"}
        print(
            f"[BatchedPredictUnit:predict] done x", args, kwargs, r.keys(), flush=True
        )
        return r

    def __getattr__(self, item):  # pragma: no cover
        if item in {"dataset_to_tasks", "device"}:
            return self.__dict__[item]
        raise AttributeError(item)


# --- Public helpers used by the ingress / API layer --------------------------


def install_predict_handle(handle) -> None:
    """Install the UMA deployment handle and build the calculator.

    Call this ONCE from your Serve ingress deployment __init__, passing the
    bound handle for `_PredictDeploy`. Example (in your serve_app.py):

        uma = _PredictDeploy.options(name=UMA_DEPLOYMENT_NAME).bind(MODEL_NAME, TASK_NAME)
        ing = Ingress.bind(uma)
        serve.run(ing, name="http_app", route_prefix="/")

    And inside Ingress.__init__(predict_handle): install_predict_handle(predict_handle)
    """
    global _PU, _CALC

    dtt = pretrained_mlip.get_predict_unit(
        MODEL_NAME,
        device="cpu",
    ).dataset_to_tasks

    _PU = BatchedPredictUnit(handle, dataset_to_tasks=dtt)
    _CALC = FAIRChemCalculator(_PU, task_name=TASK_NAME)


def get_batched_predict_unit() -> BatchedPredictUnit:
    if _PU is None:
        raise RuntimeError(
            "UMA handle not installed. Ensure your Serve ingress called install_predict_handle(...)"
        )
    return _PU


def get_calculator() -> FAIRChemCalculator:
    return FAIRChemCalculator(_PU, task_name=TASK_NAME)
    if _CALC is None:
        raise RuntimeError(
            "Calculator not initialized. Ensure your Serve ingress called install_predict_handle(...)"
        )
    return _CALC


def health_snapshot():
    return {
        "model": MODEL_NAME,
        "task": TASK_NAME,
        "device": DEVICE,
        "cuda_available": torch.cuda.is_available(),
        "model_loaded": bool(_PU is not None),
    }


__all__ = [
    "MODEL_NAME",
    "TASK_NAME",
    "DEVICE",
    "UMA_DEPLOYMENT_NAME",
    "_PredictDeploy",
    "install_predict_handle",
    "get_batched_predict_unit",
    "get_calculator",
    "health_snapshot",
]
