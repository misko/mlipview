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
import time
from typing import Any, List, Tuple

import numpy as np

# 'ray' module reference is via ray.serve imported above; avoid unused
# direct import
import torch
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch
from fairchem.core.units.mlip_unit import InferenceSettings
from ray import serve

# --- Config -----------------------------------------------------------------

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

        self.DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
        print(
            f"[model_runtime] torch_version={getattr(torch, '__version__', None)} "
            f"torch_cuda_version={getattr(torch.version, 'cuda', None)} "
            f"cuda_available={torch.cuda.is_available()} "
            f"device={self.DEVICE}",
            flush=True,
        )
        print(f"[batched:init] loading model={model_name} task={task_name}")
        print(
            (
                f"[batched:init] DEVICE={self.DEVICE} "
                f"cuda_available={torch.cuda.is_available()} "
            ),
            flush=True,
        )
        # Simple in-replica metrics for observability
        self._predict_calls = 0  # number of predict() invocations (batches)
        self._predict_items = 0  # total items processed across all batches
        # Geometry debug toggle: when enabled, log positions and cell
        # for each predict call
        self._geom_debug = os.getenv("UMA_GEOM_DEBUG", "0") in (
            "1",
            "true",
            "TRUE",
            "True",
        )
        self._unit = pretrained_mlip.get_predict_unit(
            model_name,
            device=self.DEVICE,
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
        if self._geom_debug:
            print(
                "[UMA] predict called on device %s size=%d"
                % (self.DEVICE, len(payloads)),
                flush=True,
            )
        # Optional geometry debug logging
        if self._geom_debug:
            singles = [p[0][0] for p in payloads] if payloads else []
            for idx, item in enumerate(singles):
                pos = getattr(item, "pos", None)
                cell = getattr(item, "cell", None)
                ds = getattr(item, "dataset", None)
                zs = getattr(item, "atomic_numbers", None)
                n = pos.shape[0]
                print(
                    f"[UMA][geom] item={idx} natoms={n} "
                    f"pos={pos} cell={cell} "
                    f"Z={zs} dataset={ds}",
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
                # Geometry debug: log outputs (energy, forces)
                if self._geom_debug:
                    e0 = res[0].get("energy")
                    f0 = res[0].get("forces")
                    print(
                        "[UMA][geom][out] item=0 E=" + str(e0),
                        flush=True,
                    )
                dt = time.perf_counter() - t0
                if self._geom_debug:
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
                forces_by_mol = all_outputs["forces"].split(batch["natoms"].tolist())
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
                # Geometry debug: log outputs (energy, forces) for each item
                if self._geom_debug:
                    for idx, it in enumerate(out):
                        e0 = it.get("energy")
                        f0 = it.get("forces")
                        print(
                            f"[UMA][geom][out] item={idx} E=" + str(e0),
                            flush=True,
                        )
            dt = time.perf_counter() - t0
            size = len(payloads)
            if self._geom_debug:
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
        snap = {
            "calls": int(self._predict_calls),
            "items": int(self._predict_items),
            "device": self.DEVICE,
            "model": MODEL_NAME,
            "task": TASK_NAME,
        }
        print("[batched:get_stats] " + str(snap), flush=True)
        return snap

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

        self.DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
        print(
            f"[model_runtime] torch_version="
            f"{getattr(torch, '__version__', None)} "
            f"torch_cuda_version={getattr(torch.version, 'cuda', None)} "
            f"cuda_available={torch.cuda.is_available()} "
            f"device={self.DEVICE}",
            flush=True,
        )

    def predict(self, *args, **kwargs):
        # print("[UMA-client] calling remote predict", flush=True)
        t0 = time.perf_counter()
        # NOTE: Called from sync code (e.g., FastAPI handler in a thread pool).
        # DeploymentResponse.result() is safe in that context.
        resp = self._handle.predict.remote((args, kwargs))  # type: ignore
        r = resp.result()
        dt = time.perf_counter() - t0
        # print(f"[UMA-client] predict wall={dt:.4f}s", flush=True)
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
            raise RuntimeError(f"UMA predict returned unexpected type: {type(r0)}")
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

    print(
        "[install_predict_handle] building CPU metadata predict_unit for "
        "dataset_to_tasks",
        flush=True,
    )
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
    # Default to ingress/local process view
    snap = {
        "model": MODEL_NAME,
        "task": TASK_NAME,
        "device": "cuda" if torch.cuda.is_available() else "cpu",
        "cuda_available": torch.cuda.is_available(),
        "model_loaded": bool(_PU is not None),
        "ingress_device": "cuda" if torch.cuda.is_available() else "cpu",
    }
    # If UMA handle installed, attempt to read device from the GPU replica
    try:
        if _PU is not None:
            # Call replica's get_stats (async on deployment) synchronously
            resp = _PU._handle.get_stats.remote()  # type: ignore[attr-defined]
            stats = resp.result()
            if isinstance(stats, dict) and "device" in stats:
                snap["device"] = stats.get("device")
                snap["replica_device"] = stats.get("device")
                snap["replica_stats_calls"] = stats.get("calls")
                snap["replica_stats_items"] = stats.get("items")
                snap["cuda_available"] = bool(stats.get("device") == "cuda")
    except Exception as e:  # pragma: no cover - best-effort health path
        snap["replica_error"] = str(e)
    print("[health_snapshot] " + str(snap), flush=True)
    return snap


__all__ = [
    "MODEL_NAME",
    "TASK_NAME",
    "UMA_DEPLOYMENT_NAME",
    "_PredictDeploy",
    "install_predict_handle",
    "get_batched_predict_unit",
    "get_calculator",
    "health_snapshot",
]


def _to_list3(x):
    try:
        if hasattr(x, "detach") and callable(x.detach):
            x = x.detach()
        if hasattr(x, "cpu") and callable(getattr(x, "cpu")):
            x = x.cpu()
        if hasattr(x, "numpy") and callable(getattr(x, "numpy")):
            x = x.numpy()
    except Exception:
        pass
    try:
        arr = np.asarray(x, dtype=float)
        return arr.tolist()
    except Exception:
        # best effort for lists
        return x


def _to_list_int(x):
    """
    Convert tensor/array/list of ints to plain Python
    list[int].
    """
    try:
        if hasattr(x, "detach") and callable(x.detach):
            x = x.detach()
        if hasattr(x, "cpu") and callable(getattr(x, "cpu")):
            x = x.cpu()
        if hasattr(x, "numpy") and callable(getattr(x, "numpy")):
            x = x.numpy()
    except Exception:
        pass
    try:
        arr = np.asarray(x, dtype=int)
        return arr.tolist()
    except Exception:
        # best effort for lists
        return x


def _short(x, limit: int = 400):
    try:
        s = str(x)
    except Exception:
        try:
            s = repr(x)
        except Exception:
            s = f"<unprintable:{type(x).__name__}>"
    if len(s) <= limit:
        return s
    return s[:limit] + "...(truncated)"
