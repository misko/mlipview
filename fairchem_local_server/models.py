"""Central data models & enums for UMA local server API.

Separated from server.py to allow reuse across direct FastAPI and Ray Serve
ingress without circular imports.
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel


class RelaxCalculatorName(str, Enum):
    uma = "uma"
    lj = "lj"


class SimpleIn(BaseModel):
    atomic_numbers: List[int]
    coordinates: List[List[float]]  # N x 3
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    properties: Optional[List[str]] = None
    cell: Optional[List[List[float]]] = None  # 3x3 cell matrix or None
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma


class RelaxIn(BaseModel):
    atomic_numbers: List[int]
    coordinates: List[List[float]]
    steps: int = 20
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma
    cell: Optional[List[List[float]]] = None
    pbc: Optional[List[bool]] = None
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    fmax: float = 0.05
    max_step: float = 0.2
    optimizer: Optional[str] = "bfgs"
    optimizer_params: Optional[Dict[str, Any]] = None
    return_trace: bool = False


class RelaxResult(BaseModel):
    initial_energy: float
    final_energy: float
    positions: List[List[float]]
    forces: List[List[float]]
    stress: Optional[List[float]] = None
    steps_completed: int
    calculator: RelaxCalculatorName
    trace_energies: Optional[List[float]] = None


class MDIn(BaseModel):
    atomic_numbers: List[int]
    coordinates: List[List[float]]
    steps: int = 1
    temperature: float = 298.0
    timestep_fs: float = 1.0
    friction: float = 0.02
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma
    cell: Optional[List[List[float]]] = None
    pbc: Optional[List[bool]] = None
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    return_trajectory: bool = False


class MDResult(BaseModel):
    initial_energy: float
    final_energy: float
    positions: List[List[float]]
    velocities: List[List[float]]
    forces: List[List[float]]
    steps_completed: int
    temperature: float
    energies: Optional[List[float]] = None
    calculator: RelaxCalculatorName


__all__ = [
    "RelaxCalculatorName",
    "SimpleIn",
    "RelaxIn",
    "RelaxResult",
    "MDIn",
    "MDResult",
]
