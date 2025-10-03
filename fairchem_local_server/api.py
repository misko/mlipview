"""Pure function API implementations (no FastAPI or Ray specifics)."""

from __future__ import annotations

import os
from typing import List

from ase.calculators.lj import LennardJones
from ase.optimize import BFGS as _BFGS
from fastapi import HTTPException

from .atoms_utils import build_atoms, compute_properties
from .log import log_event
from .model_runtime import ensure_model_loaded, get_calculator
from .models import MDIn, MDResult, RelaxCalculatorName, RelaxIn, RelaxResult, SimpleIn


def _attach_calc(atoms, which: RelaxCalculatorName):
    ensure_model_loaded()
    if which == RelaxCalculatorName.uma:
        atoms.calc = get_calculator()
    elif which == RelaxCalculatorName.lj:
        atoms.calc = LennardJones(rc=3.0)
    else:  # pragma: no cover
        raise HTTPException(status_code=400, detail="Unknown calculator")


def simple_calculate(inp: SimpleIn):
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    _attach_calc(atoms, inp.calculator)
    props = tuple(inp.properties or ("energy", "forces"))
    results = compute_properties(atoms, props)
    log_event(
        "simple_calc",
        natoms=len(inp.atomic_numbers),
        props=list(props),
        stress=bool(results.get("stress") is not None),
    )
    return {"results": results}


def relax(inp: RelaxIn) -> RelaxResult:
    if inp.steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        pbc=inp.pbc,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    _attach_calc(atoms, inp.calculator)

    try:
        initial_energy = float(atoms.get_potential_energy())
    except Exception as ee:
        raise HTTPException(status_code=500, detail=f"Initial energy failed: {ee}")
    opt = _BFGS(atoms, logfile=None, maxstep=float(inp.max_step or 0.2))
    trace: List[float] = []
    steps_completed = 0
    for step_idx in range(int(inp.steps)):
        try:
            opt.step()
            steps_completed += 1
            if inp.return_trace:
                try:
                    trace.append(float(atoms.get_potential_energy()))
                except Exception:
                    pass
            if os.environ.get("RELAX_VERBOSE") == "1":  # optional debug
                ene = float(atoms.get_potential_energy())
                log_event(
                    "relax_step",
                    step=step_idx + 1,
                    total=inp.steps,
                    energy=ene,
                    calc=inp.calculator,
                )
        except Exception as ste:
            log_event("relax_step_error", step=step_idx + 1, error=str(ste))
            break
    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    try:
        stress = atoms.get_stress().tolist()  # type: ignore
    except Exception:
        stress = None
    log_event(
        "relax_done",
        natoms=len(inp.atomic_numbers),
        steps=steps_completed,
        initial=initial_energy,
        final=final_energy,
        trace_len=len(trace),
        calc=inp.calculator,
    )
    return RelaxResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        forces=forces,
        stress=stress,
        steps_completed=steps_completed,
        calculator=inp.calculator,
        trace_energies=(trace if inp.return_trace else None),
    )


def md_step(inp: MDIn) -> MDResult:
    if inp.steps <= 0:
        raise HTTPException(status_code=400, detail="steps must be >0")
    atoms = build_atoms(
        inp.atomic_numbers,
        inp.coordinates,
        cell=inp.cell,
        pbc=inp.pbc,
        charge=inp.charge or 0,
        spin=inp.spin_multiplicity or 1,
    )
    _attach_calc(atoms, inp.calculator)
    import numpy as np
    from ase import units as _units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    MaxwellBoltzmannDistribution(atoms, temperature_K=inp.temperature)
    try:
        initial_energy = float(atoms.get_potential_energy())
    except Exception as ee:
        raise HTTPException(status_code=500, detail=f"Initial energy failed: {ee}")
    dyn = Langevin(
        atoms,
        float(inp.timestep_fs) * _units.fs,
        temperature_K=inp.temperature,
        friction=inp.friction,
        logfile=None,
    )
    energies = [] if inp.return_trajectory else None
    prev = atoms.get_positions().copy()
    for _ in range(int(inp.steps)):
        dyn.run(1)
        if energies is not None:
            try:
                energies.append(float(atoms.get_potential_energy()))
            except Exception:
                energies.append(float("nan"))
        new_pos = atoms.get_positions()
        max_disp = np.sqrt(((new_pos - prev) ** 2).sum(axis=1)).max()
        if not np.isfinite(max_disp) or max_disp > 5.0:
            raise HTTPException(
                status_code=500,
                detail=("MD instability detected (max step disp " f"{max_disp:.2f} Å)"),
            )
        prev[:] = new_pos
    final_energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces().tolist()
    velocities = atoms.get_velocities().tolist()
    KE = float(atoms.get_kinetic_energy())
    Tfinal = (2.0 * KE) / (3.0 * len(inp.atomic_numbers) * _units.kB)
    log_event(
        "md_done",
        natoms=len(inp.atomic_numbers),
        steps=inp.steps,
        initial=initial_energy,
        final=final_energy,
        T=Tfinal,
        calc=inp.calculator,
    )
    return MDResult(
        initial_energy=initial_energy,
        final_energy=final_energy,
        positions=atoms.get_positions().tolist(),
        velocities=velocities,
        forces=forces,
        steps_completed=int(inp.steps),
        temperature=Tfinal,
        energies=energies,
        calculator=inp.calculator,
    )


__all__ = [
    "simple_calculate",
    "relax",
    "md_step",
]
