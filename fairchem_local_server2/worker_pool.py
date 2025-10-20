from __future__ import annotations

import itertools
import time
from typing import Any, Dict, List, Optional

import numpy as np
import ray
from ase import units as _units
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS as _BFGS
from fastapi import HTTPException

from fairchem_local_server.atoms_utils import (
    build_atoms,
    center_and_return_shift,
    compute_properties,
)
from fairchem_local_server.log import log_event
from fairchem_local_server.model_runtime import get_calculator, install_predict_handle
from fairchem_local_server.models import (
    MDResult,
    PrecomputedValues,
    RelaxCalculatorName,
    RelaxResult,
    SimpleIn,
)


def _validate_atomic_numbers_or_raise(atomic_numbers: List[int]) -> None:
    bad_idx = [i for i, z in enumerate(atomic_numbers) if int(z) <= 0]
    if bad_idx:
        bad_vals = [int(atomic_numbers[i]) for i in bad_idx]
        raise ValueError(
            ("[INVALID_ATOM_Z] atomic_numbers must be positive " "integers; ")
            + f"got invalid values at indices {bad_idx}: {bad_vals}"
        )


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
        print("RUNNING MD WITH VELOCITIES:", velocities, flush=True)
        t0 = time.perf_counter()
        calc_enum = RelaxCalculatorName(calculator)
        _validate_atomic_numbers_or_raise(atomic_numbers)
        atoms = build_atoms(atomic_numbers, positions, cell=cell)
        atoms.calc = get_calculator() if calc_enum == RelaxCalculatorName.uma else None
        print("RUNNING MD WITH VELOCITIES2:", velocities, flush=True)
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
        print(
            (
                f"[timing] ASEWorker.run_md natoms={len(atomic_numbers)} "
                f"calc={calc_enum.value} wall={dt:.4f}s"
            ),
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
        atoms.calc = get_calculator() if calc_enum == RelaxCalculatorName.uma else None
        rx = _relax_run(
            atoms,
            steps=int(steps),
            calculator=calc_enum,
            max_step=float(max_step),
            precomputed=None,
        )
        out = rx.dict()
        dt = time.perf_counter() - t0
        print(
            (
                f"[timing] ASEWorker.run_relax natoms={len(atomic_numbers)} "
                f"calc={calc_enum.value} wall={dt:.4f}s"
            ),
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
        from fairchem_local_server.models import SimpleIn

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
        out = res
        dt = time.perf_counter() - t0
        print(
            (
                f"[timing] ASEWorker.run_simple natoms="
                f"{len(atomic_numbers)} calc={calculator} wall={dt:.4f}s"
            ),
            flush=True,
        )
        return out


class WorkerPool:
    def __init__(self, size: int, uma_handle=None):
        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
        self._actors = [ASEWorker.remote(uma_handle) for _ in range(max(1, int(size)))]
        self._rr = itertools.cycle(self._actors)

    def any(self):  # choose a worker (round-robin)
        return next(self._rr)


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
    velocities_in,
) -> MDResult:
    t_start = time.perf_counter()
    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")
    if steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")
    _attach_calc(atoms, calculator)
    shift = center_and_return_shift(atoms)

    # If client supplied velocities, use them directly; else initialize from
    # temperature.
    if (atoms.get_velocities() > 0.0).any():
        pass  # already has velocities on atoms
    elif velocities_in is not None:
        import numpy as _np  # type: ignore

        try:
            v = _np.array(velocities_in, dtype=float)
        except Exception as ve:  # pragma: no cover - defensive
            raise HTTPException(
                status_code=400,
                detail=f"invalid velocities: {ve}",
            )
        if v.shape != (len(atoms), 3):
            raise HTTPException(
                status_code=400,
                detail="velocities shape mismatch",
            )
        if not _np.all(_np.isfinite(v)):
            raise HTTPException(
                status_code=400,
                detail="velocities contain non-finite",
            )
        atoms.set_velocities(v)
    else:
        # print("MD: initializing velocities from temperature", temperature)
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        # MaxwellBoltzmannDistribution(atoms, temperature_K=0.0)
        # print(atoms.get_velocities())
    # velocities newly initialized from temperature

    pre_applied: list[str] = _maybe_apply_precomputed(atoms, precomputed, len(atoms))

    if "energy" in pre_applied:
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore
    else:
        initial_energy = float(atoms.get_potential_energy())

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
        max_disp = float(np.sqrt(((new_pos - prev) ** 2).sum(axis=1)).max())
        if not np.isfinite(max_disp) or max_disp > 5.0:
            raise HTTPException(
                status_code=500,
                detail=(f"MD instability detected (max step disp {max_disp:.2f} Å)"),
            )
        prev[:] = new_pos

    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    velocities = atoms.get_velocities().tolist()
    KE = float(atoms.get_kinetic_energy())
    Tfinal = (2.0 * KE) / (3.0 * len(atoms) * _units.kB)

    # log_event(
    #     "md_done",
    #     natoms=len(atoms),
    #     steps=steps,
    #     initial=initial_energy,
    #     final=final_energy,
    #     T=Tfinal,
    #     calc=calculator,
    #     precomputed=bool(pre_applied),
    #     precomputed_keys=pre_applied,
    #     reused_velocities=reused_velocities,
    # )

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
        energies=energies,
        calculator=calculator,
        precomputed_applied=(pre_applied if pre_applied else None),
    )
    dt = time.perf_counter() - t_start
    nat = len(atoms)
    print(
        (
            f"[timing] _md_run natoms={nat} steps={steps} "
            f"calc={calculator} wall={dt:.4f}s"
        ),
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
    _attach_calc(atoms, calculator)

    shift = center_and_return_shift(atoms)

    # Apply precomputed results (if any) before first energy access
    pre_applied: list[str] = _maybe_apply_precomputed(atoms, precomputed, len(atoms))

    if "energy" in pre_applied:
        # Honor client-provided energy without triggering a recalculation
        # that would overwrite calc.results. Forces may also be injected and
        # should persist until first calculator call by optimizer.
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore
    else:
        # Warm first energy evaluation (some calculators may have lazy init)
        initial_energy = float(atoms.get_potential_energy())

    opt = _BFGS(atoms, logfile=None, maxstep=float(max_step))

    steps_completed = 0

    # Execute the full requested number of BFGS steps.
    # (Previously capped at 10 to limit runtime; parity tests need full count.)
    for _ in range(int(steps)):
        opt.step()
        steps_completed += 1

    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    try:
        stress = atoms.get_stress().tolist()  # type: ignore
    except Exception:
        stress = None

    # log_event(
    #     "relax_done",
    #     natoms=len(atoms),
    #     steps=steps_completed,
    #     initial=initial_energy,
    #     final=final_energy,
    #     trace_len=(len(trace) if trace_enabled else 0),
    #     calc=calculator,
    #     precomputed=bool(pre_applied),
    #     precomputed_keys=pre_applied,
    # )

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
    nat = len(atoms)
    print(
        f"[timing] _relax_run natoms={nat} steps={steps} wall={dt:.4f}s",
        flush=True,
    )
    return res


def _simple_calculate(inp: SimpleIn):
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


def _simple_run(atoms, properties: List[str], calculator: RelaxCalculatorName):
    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")
    _attach_calc(atoms, calculator)

    shift = center_and_return_shift(atoms)

    props = tuple(properties or ("energy", "forces"))
    results = compute_properties(atoms, props)
    if shift is not None:
        atoms.set_positions(atoms.get_positions() - shift)
        atoms.set_cell(None)
    # log_event(
    #     "simple_calc",
    #     natoms=len(atoms),
    #     props=list(props),
    #     stress=bool(results.get("stress") is not None),
    # )
    return {"results": results}


def _attach_calc(atoms, which: RelaxCalculatorName):
    if which == RelaxCalculatorName.uma:
        atoms.calc = get_calculator()
    elif which == RelaxCalculatorName.lj:
        atoms.calc = LennardJones(rc=3.0)
    else:  # pragma: no cover
        raise HTTPException(status_code=400, detail="Unknown calculator")


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
    import numpy as np

    calc = getattr(atoms, "calc", None)
    if calc is None:
        return []

    if pre.energy is not None:
        e = float(pre.energy)
        if not np.isfinite(e):
            raise HTTPException(
                status_code=400,
                detail="precomputed.energy not finite",
            )
        calc.results["energy"] = e
        calc.results["free_energy"] = e
        applied.append("energy")

    if pre.forces is not None:
        f = np.array(pre.forces, dtype=float)
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
        calc.atoms = atoms.copy()

    if pre.stress is not None:
        s = np.array(pre.stress, dtype=float)
        if s.shape == (6,):
            pass  # already Voigt
        elif s.shape == (9,):
            m = s.reshape(3, 3)
            # ASE Voigt order: xx, yy, zz, yz, xz, xy
            s = np.array(
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
                status_code=400,
                detail="precomputed.stress contain non-finite",
            )
        calc.results["stress"] = s
        applied.append("stress")
    return applied
