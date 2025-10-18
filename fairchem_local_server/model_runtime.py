"""Central runtime utilities for UMA model + ASE calculator.

Serve-owns-HTTP model:
- This module defines the UMA batched predictor deployment (no HTTP)
- The Serve ingress app constructs the DAG and injects a handle via
    install_predict_handle(...)
- No serve.start / serve.run in this file

Notes:
- Server-side Atoms cache has been removed across the stack.
- We do not keep a global FAIRChemCalculator instance; a lightweight
    calculator is constructed on demand using the installed UMA handle.
"""

from __future__ import annotations

import os
from typing import Any, List, Tuple
import time

# 'ray' module reference is via ray.serve imported above; avoid unused
# direct import
import torch
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch
from fairchem.core.units.mlip_unit import InferenceSettings
from ray import serve

# --- Config -----------------------------------------------------------------

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
MODEL_NAME = os.getenv("UMA_MODEL", "uma-s-1p1")
TASK_NAME = os.getenv("UMA_TASK", "omol")

# Batch tunables
MAX_BATCH = int(os.environ.get("UMA_BATCH_MAX", 16))
WAIT_S = float(os.environ.get("UMA_BATCH_WAIT_S", 0.003))

# Logical deployment name (must match the one you bind in the Serve DAG)
UMA_DEPLOYMENT_NAME = "uma_predict"

# --- State ------------------------------------------------------------------

_PU: BatchedPredictUnit | None = None


# --- Ray Serve UMA deployment (no HTTP route) --------------------------------


@serve.deployment(ray_actor_options={"num_gpus": 1})
class _PredictDeploy:  # runs on GPU replica
    def __init__(self, model_name: str, task_name: str):
        print(f"[batched:init] loading model={model_name} task={task_name}")
        # Simple in-replica metrics for observability
        self._predict_calls = 0  # number of predict() invocations (batches)
        self._predict_items = 0  # total items processed across all batches
        self._unit = pretrained_mlip.get_predict_unit(
            model_name,
            device=DEVICE,
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
        print(
            "[UMA] predict called on device %s size=%d"
            % (DEVICE, len(payloads)),
            flush=True,
        )
        t0 = time.perf_counter()
        try:
            # update simple metrics
            self._predict_calls += 1
            self._predict_items += int(len(payloads))
            out: List[Any] = []
            if len(payloads) == 1:
                single = payloads[0][0][0]
                pred = self._unit.predict(single)
                res = [{k: v.detach().cpu() for k, v in pred.items()}]
                dt = time.perf_counter() - t0
                print(
                    (
                        "[UMA] predict finished size=1 "
                        f"wall={dt:.4f}s "
                        f"ms/item={(dt * 1000.0):.2f}"
                    ),
                    flush=True,
                )
                return res
            else:
                # warmup
                batch = atomicdata_list_to_batch([x[0][0] for x in payloads])
                batch.dataset = [x[0] for x in batch.dataset]
                pred = self._unit.predict(batch)
                all_outputs = {k: v.detach().cpu() for k, v in pred.items()}
                forces_by_mol = all_outputs["forces"].split(
                    batch["natoms"].tolist()
                )
                # stress_by_mol = all_outputs["stress"].split(
                #     batch["num_atoms"]
                # )
                out = [
                    {"energy": energy, "forces": forces, "stress": stress[0]}
                    for energy, forces, stress in zip(
                        all_outputs["energy"].split(1),
                        forces_by_mol,
                        all_outputs["stress"].split(1),
                    )
                ]
            dt = time.perf_counter() - t0
            size = len(payloads)
            print(
                (
                    f"[UMA] predict finished size={size} "
                    f"wall={dt:.4f}s "
                    f"ms/item={(dt * 1000.0 / max(1, size)):.2f}"
                ),
                flush=True,
            )
            return out
        except Exception as e:
            print("Unexpected error in UMA predictor: %s" % e)
            raise

    # Non-batched methods for simple stats/health observability
    async def get_stats(self) -> dict:
        return {
            "calls": int(self._predict_calls),
            "items": int(self._predict_items),
            "device": DEVICE,
            "model": MODEL_NAME,
            "task": TASK_NAME,
        }

    async def reset_stats(self) -> dict:
        self._predict_calls = 0
        self._predict_items = 0
        return {"ok": True}


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
        self.device = DEVICE

    def predict(self, *args, **kwargs):
        print("[UMA-client] calling remote predict", flush=True)
        t0 = time.perf_counter()
        # NOTE: Called from sync code (e.g., FastAPI handler in a thread pool).
        # DeploymentResponse.result() is safe in that context.
        resp = self._handle.predict.remote((args, kwargs))  # type: ignore
        r = resp.result()
        dt = time.perf_counter() - t0
        print(f"[UMA-client] predict wall={dt:.4f}s", flush=True)
        # Serve batched methods always return a list of results with the same
        # length as the input batch. For our client wrapper, unwrap a single
        # element list to a dict for downstream FAIRChemCalculator.
        if isinstance(r, list):
            if len(r) == 0:
                raise RuntimeError("UMA predict returned empty list")
            r0 = r[0]
        else:
            r0 = r
        if not isinstance(r0, dict):
            raise RuntimeError(
                f"UMA predict returned unexpected type: {type(r0)}"
            )
        # basic schema check
        for k in ("energy", "forces", "stress"):
            if k not in r0:
                raise RuntimeError(f"UMA predict missing key: {k}")
        return r0

    def __getattr__(self, item):  # pragma: no cover
        if item in {"dataset_to_tasks", "device"}:
            return self.__dict__[item]
        raise AttributeError(item)


# --- Public helpers used by the ingress / API layer --------------------------


def install_predict_handle(handle) -> None:
    """Install the UMA deployment handle and build the calculator.

    Call this ONCE from your Serve ingress deployment __init__, passing the
    bound handle for `_PredictDeploy`. Example (in your serve_app.py):

        uma = _PredictDeploy.options(name=UMA_DEPLOYMENT_NAME).bind(
            MODEL_NAME, TASK_NAME
        )
        ing = Ingress.bind(uma)
        serve.run(ing, name="http_app", route_prefix="/")

    And inside Ingress.__init__(predict_handle):
        install_predict_handle(predict_handle)
    """
    global _PU

    dtt = pretrained_mlip.get_predict_unit(
        MODEL_NAME,
        device="cpu",
    ).dataset_to_tasks

    _PU = BatchedPredictUnit(handle, dataset_to_tasks=dtt)


def get_batched_predict_unit() -> BatchedPredictUnit:
    if _PU is None:
        raise RuntimeError(
            "UMA handle not installed. Ensure your Serve ingress called "
            "install_predict_handle(...)"
        )
    return _PU


def get_calculator() -> FAIRChemCalculator:
    if _PU is None:
        raise RuntimeError(
            "UMA handle not installed. Ensure your Serve ingress called "
            "install_predict_handle(...)"
        )
    # Construct a lightweight calculator bound to the UMA predictor handle.
    # FAIRChemCalculator is inexpensive to instantiate and holds no heavy
    # state.
    return FAIRChemCalculator(_PU, task_name=TASK_NAME)


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
