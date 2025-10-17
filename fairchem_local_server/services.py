"""Pure service functions (no FastAPI / Ray specifics)."""

from __future__ import annotations

from typing import List

import numpy as np
from ase import units as _units
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import BFGS as _BFGS
from fastapi import HTTPException

from .atoms_utils import (
    build_atoms,
    center_and_return_shift,
    compute_properties,
)
from .log import log_event
from .models import (
    MDIn,
    MDResult,
    PrecomputedValues,
    RelaxCalculatorName,
    RelaxIn,
    RelaxResult,
    SimpleIn,
)


def _attach_calc(atoms, which: RelaxCalculatorName):
    if which == RelaxCalculatorName.uma:
        # Lazy import to avoid importing UMA stack unless needed
        from .model_runtime import get_calculator  # type: ignore

        atoms.calc = get_calculator()
    elif which == RelaxCalculatorName.lj:
        atoms.calc = LennardJones(rc=3.0)
    else:  # pragma: no cover
        raise HTTPException(status_code=400, detail="Unknown calculator")


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


def simple_calculate(inp: SimpleIn):
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


# Cache-removed: simple_from_cache endpoint/function deleted


def _relax_run(
    atoms,
    steps: int,
    calculator: RelaxCalculatorName,
    max_step: float,
    return_trace: bool,
    precomputed: PrecomputedValues | None,
) -> RelaxResult:
    if len(atoms) == 0:
        raise HTTPException(status_code=400, detail="No atoms provided")
    if steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")
    _attach_calc(atoms, calculator)

    shift = center_and_return_shift(atoms)

    # Apply precomputed results (if any) before first energy access
    pre_applied: list[str] = _maybe_apply_precomputed(
        atoms, precomputed, len(atoms)
    )

    if "energy" in pre_applied:
        # Honor client-provided energy without triggering a recalculation
        # that would overwrite calc.results. Forces may also be injected and
        # should persist until first calculator call by optimizer.
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore
    else:
        # Warm first energy evaluation (some calculators may have lazy init)
        try:
            initial_energy = float(atoms.get_potential_energy())
        except Exception as ee:
            raise HTTPException(
                status_code=500,
                detail=f"Initial energy failed: {ee}",
            )

    opt = _BFGS(atoms, logfile=None, maxstep=float(max_step))

    trace_enabled = bool(return_trace)
    trace: List[float] = []
    steps_completed = 0

    # Execute the full requested number of BFGS steps.
    # (Previously capped at 10 to limit runtime; parity tests need full count.)
    for step_idx in range(int(steps)):
        try:
            opt.step()
            steps_completed += 1

            if trace_enabled:
                try:
                    trace.append(float(atoms.get_potential_energy()))
                except Exception:
                    # Keep trace length consistent even on transient failures
                    trace.append(float("nan"))

        except Exception as ste:
            log_event("relax_step_error", step=step_idx + 1, error=str(ste))
            break

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

    return RelaxResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        forces=forces,
        stress=stress,
        steps_completed=steps_completed,
        calculator=calculator,
        trace_energies=(trace if trace_enabled else None),
        precomputed_applied=(pre_applied if pre_applied else None),
    )


def relax(inp: RelaxIn) -> RelaxResult:
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    return _relax_run(
        atoms,
        int(inp.steps),
        inp.calculator,
        float(inp.max_step),
        bool(inp.return_trace),
        inp.precomputed,
    )


# Cache-removed: relax_from_cache deleted


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

    pre_applied: list[str] = _maybe_apply_precomputed(
        atoms, precomputed, len(atoms)
    )

    if "energy" in pre_applied:
        initial_energy = float(atoms.calc.results["energy"])  # type: ignore
    else:
        try:
            initial_energy = float(atoms.get_potential_energy())
        except Exception as ee:
            raise HTTPException(
                status_code=500,
                detail=f"Initial energy failed: {ee}",
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
            try:
                energies.append(float(atoms.get_potential_energy()))
            except Exception:
                energies.append(float("nan"))

        new_pos = atoms.get_positions()
        max_disp = float(np.sqrt(((new_pos - prev) ** 2).sum(axis=1)).max())
        if not np.isfinite(max_disp) or max_disp > 5.0:
            raise HTTPException(
                status_code=500,
                detail=(
                    f"MD instability detected (max step disp {max_disp:.2f} Å)"
                ),
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

    return MDResult(
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


def md_step(inp: MDIn) -> MDResult:
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    return _md_run(
        atoms,
        steps=int(inp.steps),
        temperature=float(inp.temperature),
        timestep_fs=float(inp.timestep_fs),
        friction=float(inp.friction),
        calculator=inp.calculator,
        return_trajectory=bool(inp.return_trajectory),
        precomputed=inp.precomputed,
        velocities_in=inp.velocities,
    )


# Cache-removed: md_step_from_cache deleted


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


__all__ = [
    "simple_calculate",
    "relax",
    "md_step",
]
