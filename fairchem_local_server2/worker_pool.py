from __future__ import annotations

import itertools
import time
from typing import Any, Dict, List, Optional

import numpy as np
import ray
from ase import units as _units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS as _BFGS
from fastapi import HTTPException

from fairchem_local_server.atoms_utils import (
    build_atoms,
    center_and_return_shift,
    compute_properties,
)
from fairchem_local_server.model_runtime import get_calculator, install_predict_handle
from fairchem_local_server.models import (
    MDResult,
    PrecomputedValues,
    RelaxCalculatorName,
    RelaxResult,
    SimpleIn,
)


def _validate_atomic_numbers_or_raise(atomic_numbers: List[int]) -> None:
    """Ensure all atomic numbers are positive integers."""
    bad_idx = [i for i, z in enumerate(atomic_numbers) if int(z) <= 0]
    if bad_idx:
        bad_vals = [int(atomic_numbers[i]) for i in bad_idx]
        raise ValueError(
            "[INVALID_ATOM_Z] atomic_numbers must be positive integers; "
            f"got invalid values at indices {bad_idx}: {bad_vals}"
        )


@ray.remote(num_cpus=0.1)
class ASEWorker:
    """
    CPU worker for MD/Relax using UMA-backed calculator.

    Each call is synchronous on the worker process. The caller should await
    via ray.get(ObjectRef).
    """

    def __init__(self, handle=None):
        # Each worker process must install the UMA handle in its own module
        # state so get_calculator() works here.
        try:
            if handle is not None:
                install_predict_handle(handle)
        except Exception:
            # If Serve is not running (e.g., unit tests), ignore.
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
        precomputed: Optional[PrecomputedValues] = None,
    ) -> Dict[str, Any]:
        print(f"[run_md] with temperature={temperature}K steps={steps}", flush=True)
        calc_enum = RelaxCalculatorName(calculator)
        _validate_atomic_numbers_or_raise(atomic_numbers)

        atoms = build_atoms(atomic_numbers, positions, cell=cell)

        # Attach calculator (UMA-only; raises if not UMA)
        _attach_calc(atoms, calc_enum)

        md_res = _md_run(
            atoms,
            steps=int(steps),
            temperature=float(temperature),
            timestep_fs=float(timestep_fs),
            friction=float(friction),
            calculator=calc_enum,
            return_trajectory=False,
            precomputed=precomputed,
            velocities_in=velocities,
        )
        out = md_res.dict()
        return out

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
        precomputed: Optional[PrecomputedValues] = None,
    ) -> Dict[str, Any]:
        t0 = time.perf_counter()

        calc_enum = RelaxCalculatorName(calculator)
        _validate_atomic_numbers_or_raise(atomic_numbers)

        atoms = build_atoms(atomic_numbers, positions, cell=cell)

        # Attach calculator (UMA-only; raises if not UMA)
        _attach_calc(atoms, calc_enum)

        rx = _relax_run(
            atoms,
            steps=int(steps),
            calculator=calc_enum,
            max_step=float(max_step),
            precomputed=precomputed,
        )
        out = rx.dict()
        return out

    def run_simple(
        self,
        *,
        atomic_numbers: List[int],
        positions: List[List[float]],
        cell: Optional[List[List[float]]],
        properties: Optional[List[str]] = None,
        calculator: str = "uma",
    ) -> Dict[str, Any]:
        _validate_atomic_numbers_or_raise(atomic_numbers)
        inp = SimpleIn(
            atomic_numbers=list(map(int, atomic_numbers)),
            coordinates=positions,
            cell=cell,
            properties=list(properties or ("energy", "forces")),
            calculator=RelaxCalculatorName(calculator),
        )

        res = _simple_calculate(inp)
        return res


# at top
import threading


class WorkerPool:
    def __init__(self, size: int, uma_handle=None):
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
        n = max(1, int(size))
        self._actors = [ASEWorker.remote(uma_handle) for _ in range(n)]
        self._idx = 0
        self._n = n
        self._lock = threading.Lock()  # protects _idx

    def any(self):
        # round-robin selection, lock-protected
        with self._lock:
            a = self._actors[self._idx]
            self._idx = (self._idx + 1) % self._n
            return a


def _md_run(
    atoms,
    *,
    steps: int,
    temperature: float,
    timestep_fs: float,
    friction: float,
    calculator: RelaxCalculatorName,
    return_trajectory: bool,
    precomputed: PrecomputedValues | None,
    velocities_in: Optional[List[List[float]]],
) -> MDResult:
    t_start = time.perf_counter()

    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")
    if steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")

    # Calculator is already attached by caller for clarity/consistency.
    shift = center_and_return_shift(atoms)

    # ---- Velocities: use provided; else use existing; else initialize from T ----
    v_existing = atoms.get_velocities()
    if v_existing is not None and (np.asarray(v_existing) > 0).any():
        # Keep existing velocities on the Atoms object.
        pass
    elif velocities_in is not None:
        try:
            v = np.asarray(velocities_in, dtype=float)
        except Exception as ve:  # pragma: no cover - defensive
            raise HTTPException(status_code=400, detail=f"invalid velocities: {ve}")
        if v.shape != (len(atoms), 3):
            raise HTTPException(status_code=400, detail="velocities shape mismatch")
        if not np.all(np.isfinite(v)):
            raise HTTPException(status_code=400, detail="velocities contain non-finite")
        atoms.set_velocities(v)
    else:
        # Initialize from temperature (standard Langevin setup)
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

    # Apply any precomputed values (if provided) before first energy access.
    pre_applied: list[str] = _maybe_apply_precomputed(atoms, precomputed, len(atoms))
    if "energy" in pre_applied:
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore[attr-defined]
    else:
        initial_energy = float(atoms.get_potential_energy())

    # Langevin MD
    print(
        f"[md] Starting MD: steps={steps} T={temperature}K dt={timestep_fs}fs",
        flush=True,
    )
    dyn = Langevin(
        atoms,
        float(timestep_fs) * _units.fs,
        temperature_K=temperature,
        friction=friction,
        logfile=None,
    )

    energies = [] if return_trajectory else None
    prev = atoms.get_positions().copy()

    for _ in range(int(steps)):
        dyn.run(1)
        if energies is not None:
            energies.append(float(atoms.get_potential_energy()))

        new_pos = atoms.get_positions()
        # Guard for instabilities (large single-step displacement)
        max_disp = float(np.sqrt(((new_pos - prev) ** 2).sum(axis=1)).max())
        if not np.isfinite(max_disp) or max_disp > 5.0:
            raise HTTPException(
                status_code=500,
                detail=f"MD instability detected (max step disp {max_disp:.2f} Å)",
            )
        prev[:] = new_pos

    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    velocities = atoms.get_velocities().tolist()
    KE = float(atoms.get_kinetic_energy())
    Tfinal = (2.0 * KE) / (3.0 * len(atoms) * _units.kB)

    # Undo centering shift to return caller-space coordinates
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        atoms.set_cell(None)

    res = MDResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        velocities=velocities,
        forces=forces,
        steps_completed=int(steps),
        temperature=Tfinal,
        kinetic=KE,
        energies=energies,
        calculator=calculator,
        precomputed_applied=(pre_applied if pre_applied else None),
    )

    dt = time.perf_counter() - t_start
    print(
        f"[timing] _md_run natoms={len(atoms)} steps={steps} calc={calculator} wall={dt:.4f}s",
        flush=True,
    )
    return res


def _relax_run(
    atoms,
    steps: int,
    calculator: RelaxCalculatorName,
    max_step: float,
    precomputed: PrecomputedValues | None,
) -> RelaxResult:
    t_start = time.perf_counter()

    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")
    if steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")

    # Calculator is already attached by caller for clarity/consistency.
    shift = center_and_return_shift(atoms)

    # Apply precomputed results (if any) before first energy access
    pre_applied: list[str] = _maybe_apply_precomputed(atoms, precomputed, len(atoms))

    if "energy" in pre_applied:
        # Honor client-provided energy without triggering a recalculation
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore[attr-defined]
    else:
        # Warm first energy evaluation (some calculators may have lazy init)
        initial_energy = float(atoms.get_potential_energy())

    opt = _BFGS(atoms, logfile=None, maxstep=float(max_step))

    steps_completed = 0
    for _ in range(int(steps)):
        opt.step()
        steps_completed += 1

    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    try:
        stress = atoms.get_stress().tolist()  # type: ignore[attr-defined]
    except Exception:
        stress = None

    # Undo centering shift to return caller-space coordinates
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        atoms.set_cell(None)

    res = RelaxResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        forces=forces,
        stress=stress,
        steps_completed=steps_completed,
        calculator=calculator,
        precomputed_applied=(pre_applied if pre_applied else None),
    )

    dt = time.perf_counter() - t_start
    print(
        f"[timing] _relax_run natoms={len(atoms)} steps={steps} wall={dt:.4f}s",
        flush=True,
    )
    return res


def _simple_calculate(inp: SimpleIn) -> Dict[str, Any]:
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    return _simple_run(
        atoms,
        list(inp.properties or ("energy", "forces")),
        inp.calculator,
    )


def _simple_run(
    atoms, properties: List[str], calculator: RelaxCalculatorName
) -> Dict[str, Any]:
    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")

    # UMA-only; raise otherwise
    _attach_calc(atoms, calculator)

    shift = center_and_return_shift(atoms)

    props = tuple(properties or ("energy", "forces"))
    results = compute_properties(atoms, props)

    # Undo centering shift to return caller-space coordinates
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        atoms.set_cell(None)

    return {"results": results}


def _attach_calc(atoms, which: RelaxCalculatorName) -> None:
    """
    Attach the proper calculator to the Atoms instance.

    UMA-only enforcement: anything else throws a clear, stable error string that
    the frontend already knows how to surface.
    """
    if which == RelaxCalculatorName.uma:
        atoms.calc = get_calculator()
        return

    # Explicit UMA-only guard (keeps workers aligned with WS ingress rules)
    raise HTTPException(
        status_code=400,
        detail=f"CALCULATOR_NOT_SUPPORTED: UMA_ONLY (requested: {which.value!r})",
    )


def _maybe_apply_precomputed(
    atoms, pre: PrecomputedValues | None, natoms: int
) -> list[str]:
    """Populate atoms.calc.results with client-provided precomputed values.

    Returns list of applied keys. Validation beyond basic structural checks is
    done here (finite numbers, shapes). Stress: support Voigt length-6 or
    length-9 (interpreted row-major 3x3).
    """
    if pre is None:
        return []
    applied: list[str] = []

    calc = getattr(atoms, "calc", None)
    if calc is None:
        return []

    if pre.energy is not None:
        e = float(pre.energy)
        if not np.isfinite(e):
            raise HTTPException(status_code=400, detail="precomputed.energy not finite")
        calc.results["energy"] = e
        calc.results["free_energy"] = e
        applied.append("energy")

    if pre.forces is not None:
        f = np.asarray(pre.forces, dtype=float)
        if f.shape != (natoms, 3):
            raise HTTPException(
                status_code=400, detail="precomputed.forces shape mismatch"
            )
        if not np.all(np.isfinite(f)):
            raise HTTPException(
                status_code=400, detail="precomputed.forces contain non-finite"
            )
        calc.results["forces"] = f
        applied.append("forces")
        # Some ASE calculators expect 'atoms' to match results to avoid reuse issues.
        calc.atoms = atoms.copy()

    if pre.stress is not None:
        s = np.asarray(pre.stress, dtype=float)
        if s.shape == (6,):
            pass  # already Voigt
        elif s.shape == (9,):
            m = s.reshape(3, 3)
            # ASE Voigt order: xx, yy, zz, yz, xz, xy
            s = np.asarray(
                [m[0, 0], m[1, 1], m[2, 2], m[1, 2], m[0, 2], m[0, 1]],
                dtype=float,
            )
        else:
            raise HTTPException(
                status_code=400,
                detail="precomputed.stress must have length 6 or 9",
            )
        if not np.all(np.isfinite(s)):
            raise HTTPException(
                status_code=400, detail="precomputed.stress contain non-finite"
            )
        calc.results["stress"] = s
        applied.append("stress")

    return applied
