from __future__ import annotations

import asyncio  # NEW: for least-in-flight completion watcher
import itertools
import os

# at top
import threading
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
        # Avoid CPU over-subscription from BLAS in each worker process
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

        # Each worker process must install the UMA handle in its own module state
        # so get_calculator() works here. Also cache UMA calculator once per worker.
        try:
            if handle is not None:
                install_predict_handle(handle)
            try:
                self._uma_calc = get_calculator()
            except Exception:
                self._uma_calc = None
        except Exception:
            # If Serve is not running (e.g., unit tests), ignore.
            self._uma_calc = None

    def _attach_calc_cached(self, atoms, which: RelaxCalculatorName) -> None:
        """Attach UMA calculator cheaply; enforce UMA-only."""
        if which == RelaxCalculatorName.uma:
            if getattr(self, "_uma_calc", None) is not None:
                atoms.calc = self._uma_calc
                return
            # Fallback if cache is absent for any reason.
            atoms.calc = get_calculator()
            return
        # Explicit UMA-only guard (aligned with WS ingress rules)
        raise HTTPException(
            status_code=400,
            detail=f"CALCULATOR_NOT_SUPPORTED: UMA_ONLY (requested: {which.value!r})",
        )

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
        t0 = time.perf_counter()

        calc_enum = RelaxCalculatorName(calculator)
        _validate_atomic_numbers_or_raise(atomic_numbers)

        atoms = build_atoms(atomic_numbers, positions, cell=cell)

        # Attach calculator (UMA-only), using cached UMA instance when possible.
        self._attach_calc_cached(atoms, calc_enum)

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
        out = md_res.dict()
        dt = time.perf_counter() - t0
        if os.environ.get("WS_WORKER_TIMING", "0") in ("1", "true", "TRUE", "True"):
            print(
                f"[timing] ASEWorker.run_md natoms={len(atomic_numbers)} "
                f"calc={calc_enum.value} wall={dt:.4f}s",
                flush=True,
            )
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
    ) -> Dict[str, Any]:
        t0 = time.perf_counter()

        calc_enum = RelaxCalculatorName(calculator)
        _validate_atomic_numbers_or_raise(atomic_numbers)

        atoms = build_atoms(atomic_numbers, positions, cell=cell)

        # Attach calculator (UMA-only), using cached UMA instance when possible.
        self._attach_calc_cached(atoms, calc_enum)

        rx = _relax_run(
            atoms,
            steps=int(steps),
            calculator=calc_enum,
            max_step=float(max_step),
            precomputed=None,
        )
        out = rx.dict()
        dt = time.perf_counter() - t0
        if os.environ.get("WS_WORKER_TIMING", "0") in ("1", "true", "TRUE", "True"):
            print(
                f"[timing] ASEWorker.run_relax natoms={len(atomic_numbers)} "
                f"calc={calc_enum.value} wall={dt:.4f}s",
                flush=True,
            )
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
        t0 = time.perf_counter()

        _validate_atomic_numbers_or_raise(atomic_numbers)
        inp = SimpleIn(
            atomic_numbers=list(map(int, atomic_numbers)),
            coordinates=positions,
            cell=cell,
            properties=list(properties or ("energy", "forces")),
            calculator=RelaxCalculatorName(calculator),
        )

        res = _simple_calculate(inp)
        dt = time.perf_counter() - t0
        if os.environ.get("WS_WORKER_TIMING", "0") in ("1", "true", "TRUE", "True"):
            print(
                f"[timing] ASEWorker.run_simple natoms={len(atomic_numbers)} "
                f"calc={calculator} wall={dt:.4f}s",
                flush=True,
            )
        return res


# -------------------- Least-in-flight WorkerPool --------------------


class _ActorProxy:
    """
    Lightweight proxy around a single Ray actor that keeps the pool's
    per-actor in-flight counter correct while preserving call sites:
      worker = pool.any()
      ref = worker.run_md.remote(...)
    """

    def __init__(self, pool: "WorkerPool", idx: int):
        self._pool = pool
        self._idx = idx
        self._actor = pool._actors[idx]

    # --- inflight accounting ---
    def _inc(self) -> None:
        with self._pool._lock:
            self._pool._inflight[self._idx] += 1

    def _dec(self) -> None:
        with self._pool._lock:
            self._pool._inflight[self._idx] -= 1

    def _watch_done(self, ref) -> None:
        """
        Decrement in-flight when `ref` completes.
        Uses asyncio if a running loop exists, otherwise a daemon thread.
        """
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            loop = None

        async def _aw():
            try:
                # run blocking ray.get off the event loop
                await asyncio.to_thread(ray.get, ref)
            except Exception:
                pass
            finally:
                self._dec()

        def _tw():
            try:
                ray.get(ref)
            except Exception:
                pass
            finally:
                self._dec()

        if loop is not None and loop.is_running():
            loop.create_task(_aw())
        else:
            t = threading.Thread(target=_tw, daemon=True)
            t.start()

    # --- expose .run_md/.run_relax/.run_simple with a .remote(...) method ---
    def __getattr__(self, name):
        if name in ("run_md", "run_relax", "run_simple"):

            def _submit_factory(method_name: str):
                def _submit(*args, **kwargs):
                    self._inc()
                    ref = getattr(self._actor, method_name).remote(*args, **kwargs)
                    self._watch_done(ref)
                    return ref

                class _RW:
                    def remote(self_inner, *a, **k):  # provides .remote(...)
                        return _submit(*a, **k)

                return _RW()

            return _submit_factory(name)
        # Fallback: raw actor attrs (NOTE: direct use won't update inflight)
        return getattr(self._actor, name)


class WorkerPool:
    def __init__(self, size: int, uma_handle=None):
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
        n = max(1, int(size))
        self._actors = [ASEWorker.remote(uma_handle) for _ in range(n)]
        self._n = n
        self._lock = threading.Lock()
        self._inflight: List[int] = [0] * n  # outstanding tasks per actor

    def any(self) -> _ActorProxy:
        """
        Pick the actor with the fewest outstanding tasks (least in-flight).
        Returns a proxy whose run_md/run_relax/run_simple expose .remote(...)
        identical to Ray's method handles.
        """
        with self._lock:
            # Tie-break by lowest index
            idx = min(range(self._n), key=lambda i: self._inflight[i])
        return _ActorProxy(self, idx)


# --------------------------------------------------------------------


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
    # Capture original cell before the centering routine mutates state.
    orig_cell = atoms.get_cell().copy() if atoms.cell is not None else None
    shift = center_and_return_shift(atoms)

    # ---- Velocities: use provided; else use existing; else initialize from T ----
    v_existing = atoms.get_velocities()
    if v_existing is not None and np.asarray(v_existing).shape == (len(atoms), 3):
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

    # Langevin MD (single-step intended; loop retained for API consistency)
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

    # Undo centering shift to return caller-space coordinates and restore cell.
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        if orig_cell is not None:
            atoms.set_cell(orig_cell)

    res = MDResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        velocities=velocities,
        forces=forces,
        steps_completed=int(steps),
        temperature=Tfinal,
        energies=energies,
        calculator=calculator,
        precomputed_applied=(pre_applied if pre_applied else None),
    )

    dt = time.perf_counter() - t_start
    if os.environ.get("WS_WORKER_TIMING", "0") in ("1", "true", "TRUE", "True"):
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
    orig_cell = atoms.get_cell().copy() if atoms.cell is not None else None
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

    # Undo centering shift and restore cell for caller.
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        if orig_cell is not None:
            atoms.set_cell(orig_cell)

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
    if os.environ.get("WS_WORKER_TIMING", "0") in ("1", "true", "TRUE", "True"):
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
        # Unlike MD/Relax, we don't know/care about the original cell here; leave as-is.

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
