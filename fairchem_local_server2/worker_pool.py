from __future__ import annotations

import itertools
from typing import Any, Dict, List, Optional
import ray

from fairchem_local_server.atoms_utils import build_atoms
from fairchem_local_server.model_runtime import (
    get_calculator,
    install_predict_handle,
)
from fairchem_local_server.models import RelaxCalculatorName
from fairchem_local_server.services import _md_run, _relax_run


@ray.remote(num_cpus=1)
class ASEWorker:
    """CPU worker for MD/Relax using UMA-backed calculator.

    Each call is synchronous on the worker. Caller should use ray.get for
    results.
    """

    def __init__(self, handle=None):
        # Each worker process must install the UMA handle in its own module
        # state so get_calculator() works here.
        try:
            if handle is not None:
                install_predict_handle(handle)
        except Exception:
            # If Serve is not running (e.g., unit tests with LJ), ignore.
            pass

    def run_md(
        self,
        *,
        atomic_numbers: List[int],
        positions: List[List[float]],
        velocities: Optional[List[List[float]]],
        cell: Optional[List[List[float]]],
        steps: int,
        temperature: float,
        timestep_fs: float,
        friction: float,
        calculator: str = "uma",
    ) -> Dict[str, Any]:
        calc_enum = RelaxCalculatorName(calculator)
        atoms = build_atoms(atomic_numbers, positions, cell=cell)
        atoms.calc = (
            get_calculator() if calc_enum == RelaxCalculatorName.uma else None
        )
        md_res = _md_run(
            atoms,
            steps=int(steps),
            temperature=float(temperature),
            timestep_fs=float(timestep_fs),
            friction=float(friction),
            calculator=calc_enum,
            return_trajectory=False,
            precomputed=None,
            velocities_in=velocities,
        )
        return md_res.dict()

    def run_relax(
        self,
        *,
        atomic_numbers: List[int],
        positions: List[List[float]],
        cell: Optional[List[List[float]]],
        steps: int,
        fmax: float,
        max_step: float,
        calculator: str = "uma",
    ) -> Dict[str, Any]:
        calc_enum = RelaxCalculatorName(calculator)
        atoms = build_atoms(atomic_numbers, positions, cell=cell)
        atoms.calc = (
            get_calculator() if calc_enum == RelaxCalculatorName.uma else None
        )
        rx = _relax_run(
            atoms,
            steps=int(steps),
            calculator=calc_enum,
            max_step=float(max_step),
            return_trace=False,
            precomputed=None,
        )
        return rx.dict()


class WorkerPool:
    def __init__(self, size: int, uma_handle=None):
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
        self._actors = [
            ASEWorker.remote(uma_handle) for _ in range(max(1, int(size)))
        ]
        self._rr = itertools.cycle(self._actors)

    def any(self):  # choose a worker (round-robin)
        return next(self._rr)
