"""Central data models & enums for UMA local server API.

Separated from server.py to allow reuse across direct FastAPI and Ray Serve
ingress without circular imports.
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, root_validator, validator

# Hard cap for number of atoms accepted by any request.
# Enforced at the model layer so all ingress paths share the same check.
MAX_ATOMS_PER_REQUEST: int = 40


def _enforce_atom_limit(values: Dict[str, Any]) -> Dict[str, Any]:
    """Shared check for atom count across request models.

    Expects mapping with keys 'atomic_numbers' and 'coordinates'.
    """
    zs = values.get("atomic_numbers") or []
    coords = values.get("coordinates") or []
    n = max(len(zs), len(coords))
    if n > MAX_ATOMS_PER_REQUEST:
        raise ValueError(
            ("Too many atoms: " f"{n} > MAX_ATOMS_PER_REQUEST={MAX_ATOMS_PER_REQUEST}")
        )
    return values


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

    @root_validator(skip_on_failure=True)
    def _limit_natoms(cls, values):  # type: ignore
        return _enforce_atom_limit(values)


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
    precomputed: Optional["PrecomputedValues"] = None

    @root_validator(skip_on_failure=True)
    def _limit_natoms(cls, values):  # type: ignore
        return _enforce_atom_limit(values)


class RelaxResult(BaseModel):
    initial_energy: float
    final_energy: float
    positions: List[List[float]]
    forces: List[List[float]]
    stress: Optional[List[float]] = None
    steps_completed: int
    calculator: RelaxCalculatorName
    trace_energies: Optional[List[float]] = None
    precomputed_applied: Optional[List[str]] = None


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
    precomputed: Optional["PrecomputedValues"] = None
    # Optional initial velocities (N x 3)
    velocities: Optional[List[List[float]]] = None

    @root_validator(skip_on_failure=True)
    def _limit_natoms(cls, values):  # type: ignore
        return _enforce_atom_limit(values)

    @validator("velocities")
    def _validate_velocities(cls, v, values):  # type: ignore
        if v is None:
            return v
        if not all(len(row) == 3 for row in v):
            raise ValueError("velocities must be a list of [vx,vy,vz] vectors")
        atomic_numbers = values.get("atomic_numbers")
        if atomic_numbers is not None and len(v) != len(atomic_numbers):
            raise ValueError("velocities length must match number of atoms")
        return v


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
    precomputed_applied: Optional[List[str]] = None


class PrecomputedValues(BaseModel):
    """Optional precomputed results a client can supply to skip first calc.

    energy: potential energy (eV)
    forces: N x 3 forces (eV/Å)
    stress: Voigt 6-vector (eV/Å^3) consistent with ASE ordering
            [xx, yy, zz, yz, xz, xy]. We also accept length-9 (3x3 row-major)
            and convert to Voigt ordering.
    """

    energy: Optional[float] = None
    forces: Optional[List[List[float]]] = None
    stress: Optional[List[float]] = None

    @validator("forces")
    def _validate_forces(cls, v):  # type: ignore
        if v is None:
            return v
        if not all(len(row) == 3 for row in v):
            raise ValueError("forces must be a list of [x,y,z] vectors")
        return v

    @validator("stress")
    def _validate_stress(cls, v):  # type: ignore
        if v is None:
            return v
        if len(v) not in (6, 9):
            raise ValueError("stress must have length 6 (Voigt) or 9 (matrix)")
        return v


__all__ = [
    "RelaxCalculatorName",
    "SimpleIn",
    "RelaxIn",
    "RelaxResult",
    "MDIn",
    "MDResult",
    "PrecomputedValues",
    "MAX_ATOMS_PER_REQUEST",
]
