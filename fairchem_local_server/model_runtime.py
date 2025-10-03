"""Central runtime utilities for UMA model + ASE calculator.

Provides single point for:
  * Environment variable parsing (UMA_MODEL, UMA_TASK, UMA_DEVICE)
  * CUDA-only enforcement
  * Lazy acquisition of batched predict unit
  * Construction / caching of FAIRChemCalculator
  * Health snapshot
"""

from __future__ import annotations

import os
import threading
from typing import Optional

import torch
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator

from .batched_predict import get_batched_predict_unit

MODEL_NAME = os.environ.get("UMA_MODEL", "uma-s-1p1")
TASK_NAME = os.environ.get("UMA_TASK", "omol")
_requested_device = os.environ.get("UMA_DEVICE", "cuda")
if not _requested_device.startswith("cuda"):
    raise RuntimeError("UMA_DEVICE must be 'cuda' (GPU-only mode enforced)")
DEVICE = "cuda"

_lock = threading.Lock()
_predict_unit = None
_calculator: Optional[FAIRChemCalculator] = None
_loaded = False


def ensure_model_loaded():  # idempotent
    global _predict_unit, _calculator, _loaded
    if _loaded and _calculator is not None:
        return
    with _lock:
        if _loaded and _calculator is not None:
            return
        if not torch.cuda.is_available():  # defer until GPUs visible
            raise RuntimeError("CUDA not available (GPU-only mode)")
        print(f"[runtime:init] batched unit model={MODEL_NAME} task={TASK_NAME}")
        _predict_unit = get_batched_predict_unit(MODEL_NAME, TASK_NAME)
        _calculator = FAIRChemCalculator(_predict_unit, task_name=TASK_NAME)
        _loaded = True
        print("[runtime:init] calculator ready on", DEVICE)


def get_calculator() -> FAIRChemCalculator:
    ensure_model_loaded()
    assert _calculator is not None  # for type checkers
    return _calculator


def health_snapshot():
    return {
        "model": MODEL_NAME,
        "task": TASK_NAME,
        "device": DEVICE,
        "cuda_available": torch.cuda.is_available(),
        "model_loaded": bool(_loaded),
    }


__all__ = [
    "MODEL_NAME",
    "TASK_NAME",
    "DEVICE",
    "ensure_model_loaded",
    "get_calculator",
    "health_snapshot",
]
