import os
import time
import traceback
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import torch
from ase.io.jsonio import decode, encode
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from fairchem.core.units.mlip_unit import InferenceSettings
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel


class RelaxCalculatorName(str, Enum):
    uma = "uma"
    lj = "lj"

try:
    from ase.calculators.lj import LennardJones  # standard ASE LJ
except Exception:  # pragma: no cover
    LennardJones = None  # type: ignore

try:
    from ase.optimize import BFGS  # ASE optimizer
except Exception:
    BFGS = None  # type: ignore

# --- configure model ---
MODEL_NAME = os.environ.get("UMA_MODEL", "uma-s-1p1")  # UMA small
TASK_NAME = os.environ.get("UMA_TASK", "omol")  # task: omol|omat|oc20|odac|omc
"""Server providing UMA model inference + simple LJ baseline.

Device selection: single variable DEVICE determined once from UMA_DEVICE env
or CUDA availability. No secondary _MODEL_DEVICE variable needed.
"""
_requested_device = os.environ.get("UMA_DEVICE", "cuda")
if not _requested_device.startswith("cuda"):
    raise RuntimeError("UMA_DEVICE must be 'cuda' (GPU-only mode enforced)")
DEVICE = "cuda"  # logical device target (single model design)

# Defer model load if CUDA not yet visible (e.g. importing inside Ray controller process)
_PREDICT_UNIT = None  # type: ignore
_CALCULATOR = None  # type: ignore
_MODEL_LOADED = False

def ensure_model_loaded():  # callable from endpoints / serve_app actor
    global _PREDICT_UNIT, _CALCULATOR, _MODEL_LOADED
    if _MODEL_LOADED and _CALCULATOR is not None:
        return
    if not torch.cuda.is_available():
        # Strict GPU-only mode: fail here (not at import) so Ray controller can still import module
        raise RuntimeError("CUDA not available but GPU-only mode enforced (deferred load).")
    print(f"[model:init] loading UMA model on device={DEVICE} (requested={_requested_device})")
    _PREDICT_UNIT = pretrained_mlip.get_predict_unit(
        MODEL_NAME,
        device=DEVICE,
        inference_settings=InferenceSettings(
            tf32=True,
            activation_checkpointing=False,
            merge_mole=False,
            compile=False,
            external_graph_gen=False,
            internal_graph_gen_version=2,
        ),
    )
    from fairchem.core.calculate.ase_calculator import \
        FAIRChemCalculator as _FC
    _CALCULATOR = _FC(_PREDICT_UNIT, task_name=TASK_NAME)
    _MODEL_LOADED = True
    print("[model:init] UMA model loaded and ready on", DEVICE)

# Attempt eager load if CUDA visible now (normal python tests path)
try:  # pragma: no cover (import-time branch)
    if torch.cuda.is_available():
        ensure_model_loaded()
    else:
        print("[model:init] CUDA not visible at import; deferring model load until first use")
except Exception as _e:  # logged; will re-attempt on demand
    print("[model:init] deferred load exception:", _e)

## Removed legacy cache helpers (_composition_key, get_cached_unit_and_calculator)

## Removed legacy lazy-init & device resolver; single eager load above.

########################
# Property helper
########################
def _compute_properties(atoms, props: tuple[str, ...]):
    results: Dict[str, Any] = {}
    if "energy" in props or "free_energy" in props:
        results["energy"] = float(atoms.get_potential_energy())
        results["free_energy"] = results["energy"]
    if "forces" in props:
        results["forces"] = atoms.get_forces().tolist()
    if "stress" in props:
        try:
            results["stress"] = atoms.get_stress().tolist()
        except Exception:
            results["stress"] = None
    return results


app = FastAPI(title="UMA-small ASE HTTP server")
# Re-export encode for test modules
__all__ = [
    "app",
    "encode",
    # removed helper exports
    # models for tests
    "RelaxIn",
    "RelaxResult",
    "MDIn",
    "MDResult",
    "SimpleIn",
    "RelaxCalculatorName",
]

# Allow local browser dev (adjust origins in production)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class CalcIn(BaseModel):
    atoms_json: Dict[str, Any]  # result of ase.io.jsonio.encode(Atoms)
    properties: Optional[List[str]] = None  # e.g. ["energy","forces","stress"]
    info: Optional[Dict[str, Any]] = None  # optional atoms.info overrides


class SimpleIn(BaseModel):
    atomic_numbers: list[int]
    coordinates: list[list[float]]  # N x 3
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    properties: Optional[List[str]] = None
    # Optional 3x3 cell matrix (a,b,c) row vectors
    cell: Optional[list[list[float]]] = None
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma  # uma|lj


## Duplicate RelaxCalculatorName definition removed (was previously here)

@app.get("/health")
def health():
    return {"status": "ok", "model": MODEL_NAME, "device": DEVICE, "cuda_available": torch.cuda.is_available(), "model_loaded": _MODEL_LOADED}


class RelaxIn(BaseModel):
    atomic_numbers: list[int]
    coordinates: list[list[float]]  # N x 3
    steps: int = 20  # number of BFGS steps to perform
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma
    cell: Optional[list[list[float]]] = None  # 3x3 or None
    pbc: Optional[list[bool]] = None  # length 3 or None -> infer from cell
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    fmax: float = 0.05  # convergence threshold; if <=0 disable early break
    max_step: float = 0.2  # not directly used by ASE BFGS but kept for parity
    optimizer: Optional[str] = "bfgs"  # bfgs|bfgs_ls|lbfgs|fire (extensible)
    optimizer_params: Optional[Dict[str, Any]] = None  # forwarded filtered params
    return_trace: bool = False  # whether to return per-step energy list


class RelaxResult(BaseModel):
    initial_energy: float
    final_energy: float
    positions: list[list[float]]
    forces: list[list[float]]
    stress: Optional[list[float]] = None
    steps_completed: int
    calculator: RelaxCalculatorName
    trace_energies: Optional[list[float]] = None


# ----------------------------- MD API ---------------------------------
class MDIn(BaseModel):
    atomic_numbers: list[int]
    coordinates: list[list[float]]
    steps: int = 1  # number of MD integration steps to advance
    temperature: float = 298.0  # Kelvin target (NVT)
    timestep_fs: float = 1.0  # femtoseconds per step
    friction: float = 0.02  # Langevin friction gamma (1/fs) – slightly higher for stability
    calculator: RelaxCalculatorName = RelaxCalculatorName.uma
    cell: Optional[list[list[float]]] = None
    pbc: Optional[list[bool]] = None
    charge: Optional[int] = 0
    spin_multiplicity: Optional[int] = 1
    return_trajectory: bool = False  # if true, return per-step energies & positions (lightweight)


class MDResult(BaseModel):
    initial_energy: float
    final_energy: float
    positions: list[list[float]]
    velocities: list[list[float]]
    forces: list[list[float]]
    steps_completed: int
    temperature: float  # instantaneous final temperature (kinetic)
    energies: Optional[list[float]] = None  # per-step potential energies if requested
    calculator: RelaxCalculatorName
    # No fallback indicator now: MD uses the requested calculator directly.


###############################
# NOTE: Using ASE's built-in LennardJones calculator now. The previous
# custom implementation has been removed for simplicity & correctness.
###############################


def _build_atoms_from_relax(inp: RelaxIn):  # helper
    import numpy as _np
    from ase import Atoms as _Atoms

    Z = inp.atomic_numbers
    xyz = inp.coordinates
    if len(Z) != len(xyz):
        raise HTTPException(
            status_code=400,
            detail="Length mismatch atomic_numbers vs coordinates",
        )
    pbc = (
        inp.pbc
        if inp.pbc is not None
        else ([True, True, True] if inp.cell else [False, False, False])
    )
    atoms = _Atoms(
        numbers=Z,
        positions=_np.array(xyz, dtype=float),
        cell=inp.cell,
        pbc=pbc,
    )
    atoms.info.update(
        {"charge": inp.charge or 0, "spin": inp.spin_multiplicity or 1}
    )
    return atoms


@app.post("/simple_calculate")
def simple_calculate(inp: SimpleIn):
    try:
        import time as _time

        import numpy as _np
        from ase import Atoms as _Atoms

        t0 = _time.time()
        Z = inp.atomic_numbers
        xyz = inp.coordinates
        if len(Z) != len(xyz):
            raise HTTPException(
                status_code=400,
                detail="Length mismatch atomic_numbers vs coordinates",
            )
        # Debug print cell if provided
        if inp.cell:
            try:
                rows = []
                if isinstance(inp.cell, list) and len(inp.cell) == 3:
                    for r in inp.cell:
                        if isinstance(r, list) and len(r) == 3:
                            rows.append(f"({r[0]:.3f},{r[1]:.3f},{r[2]:.3f})")
                        else:
                            rows.append(str(r))
                    print("[simple_calculate] cell=", ", ".join(rows))
            except Exception as _ce:  # noqa: F841
                print("[simple_calculate] cell print failed", _ce)
            pbc = [True, True, True]
        else:
            pbc = [False, False, False]
        atoms = _Atoms(
            numbers=Z,
            positions=_np.array(xyz, dtype=float),
            cell=inp.cell,
            pbc=pbc,
        )
        atoms.info.update(
            {
                "charge": inp.charge or 0,
                "spin": inp.spin_multiplicity or 1,
            }
        )
        props = tuple(inp.properties or ("energy", "forces"))
        # Choose calculator
        ensure_model_loaded()
        if inp.calculator == RelaxCalculatorName.lj:
            if LennardJones is None:
                raise HTTPException(
                    status_code=500, detail="ASE LJ unavailable"
                )
            atoms.calc = LennardJones(rc=3.0)
        else:  # default UMA (uma) lazy load
            atoms.calc = _CALCULATOR
        st = time.time()
        results = _compute_properties(atoms, props)
        print("[simple_calculate] inference_time_s=", round(time.time()-st, 4))
        dt = (_time.time() - t0) * 1000.0
        per_atom = dt / max(1, len(Z))
        has_stress = "stress" in results and results.get("stress") is not None
        stress_shape = (
            (
                len(results.get("stress", []))
                if isinstance(results.get("stress"), list)
                else None
            )
            if has_stress
            else None
        )
        print(
            "[simple_calculate] natoms="
            f"{len(Z)} props={props} time_ms={dt:.1f} "
            f"ms_per_atom={per_atom:.2f} "
            f"stress={'yes' if has_stress else 'no'}"
            + (f" stress_len={stress_shape}" if has_stress else "")
        )
        return {"results": results}
    except HTTPException:
        raise
    except Exception as e:
        tb = traceback.format_exc()
        print("[ERROR] /simple_calculate exception:\n", tb)
        raise HTTPException(status_code=500, detail=str(e))
    # (CalcIn fields moved above)


@app.post("/calculate")
def calculate(inp: CalcIn):
    try:
        import time as _time

        # Reconstruct Atoms (accept either JSON string or pre-parsed dict)
        t0 = _time.time()
        raw = inp.atoms_json
        if isinstance(raw, dict):
            import json as _json

            atoms = decode(_json.dumps(raw))
        else:
            atoms = decode(raw)
        if inp.info:
            atoms.info.update(inp.info)

        # Select properties to compute
        props = tuple(inp.properties or ("energy", "forces"))

        # Attach calculator, run single point
        # Determine composition for caching (Atoms with numbers attr)
        try:
            Z = list(getattr(atoms, "numbers", getattr(atoms, "get_atomic_numbers")()))
        except Exception:  # fallback if Atoms style accessor fails
            Z = []
        ensure_model_loaded()
        atoms.calc = _CALCULATOR
        results = _compute_properties(atoms, props)
        dt = (_time.time() - t0) * 1000.0
        natoms = (
            len(atoms.get_atomic_numbers())
            if hasattr(atoms, "get_atomic_numbers")
            else "n/a"
        )
        per_atom = dt / max(1, (natoms if isinstance(natoms, int) else 1))
        has_stress = "stress" in results and results.get("stress") is not None
        stress_shape = (
            (
                len(results.get("stress", []))
                if isinstance(results.get("stress"), list)
                else None
            )
            if has_stress
            else None
        )
        print(
            "[calculate] natoms="
            f"{natoms} props={props} time_ms={dt:.1f} "
            f"ms_per_atom={per_atom:.2f} "
            f"stress={'yes' if has_stress else 'no'}"
            + (f" stress_len={stress_shape}" if has_stress else "")
        )
        return {"results": results, "info": atoms.info}
    except Exception as e:
        tb = traceback.format_exc()
        print("[ERROR] /calculate exception:\n", tb)
        try:
            with open(
                os.path.join(os.path.dirname(__file__), "last_error.log"), "w"
            ) as f:
                f.write(tb)
        except Exception:
            pass
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/relax", response_model=RelaxResult)
def relax(inp: RelaxIn):
    """Run a fixed-number BFGS atomic relaxation starting from given positions.

    Always performs exactly `steps` BFGS steps (unless a step raises).
    Convergence via fmax is not enforced early (could be extended later).
    """
    try:
        if BFGS is None:  # pragma: no cover
            raise HTTPException(
                status_code=500, detail="ASE BFGS not available"
            )
        atoms = _build_atoms_from_relax(inp)
        # Select calculator
        ensure_model_loaded()
        if inp.calculator == RelaxCalculatorName.uma:
            atoms.calc = _CALCULATOR
        elif inp.calculator == RelaxCalculatorName.lj:
            if LennardJones is None:  # pragma: no cover
                raise HTTPException(
                    status_code=500, detail="ASE LJ unavailable"
                )
            # Default sigma=1, epsilon=1; cutoff 3*sigma similar to prior logic
            atoms.calc = LennardJones(rc=3.0)
        else:  # pragma: no cover - Enum guards this
            raise HTTPException(status_code=400, detail="Unknown calculator")

        # Capture initial energy (single point)
        try:
            initial_energy = float(atoms.get_potential_energy())
        except Exception as ee:
            raise HTTPException(status_code=500, detail=f"Initial energy failed: {ee}")

        # Optimizer: only BFGS kept
        from ase.optimize import BFGS as _BFGS
        optimizer = _BFGS(atoms, logfile=None, maxstep=float(inp.max_step or 0.2))
        opt_name = "bfgs"
        trace: list[float] = []

        # Run fixed number of steps (no early stopping to preserve parity)
        steps_completed = 0
        forces = None
        stress = None
        max_steps = max(0, int(inp.steps))
        per_step_debug = os.environ.get("RELAX_VERBOSE", "0") == "1" or bool(inp.return_trace)
        for step_idx in range(max_steps):
            try:
                optimizer.step()
                steps_completed += 1
                forces = atoms.get_forces()
                if inp.return_trace:
                    try:
                        trace.append(float(atoms.get_potential_energy()))
                    except Exception:
                        pass
                try:
                    stress = atoms.get_stress()
                except Exception:
                    stress = None
                if per_step_debug:
                    try:
                        E_step = float(atoms.get_potential_energy())
                    except Exception:
                        E_step = float('nan')
                    print(f"[relax:step] idx={step_idx+1}/{max_steps} calc={inp.calculator} E={E_step:.6f}")
            except Exception as ste:
                print(f"[relax] step failed after {steps_completed} steps: {ste}")
                break

        final_energy = float(atoms.get_potential_energy())
        if forces is None:
            forces = atoms.get_forces()
        print(
            f"[relax] calc={inp.calculator} opt={opt_name} natoms={len(inp.atomic_numbers)} "
            f"steps={steps_completed}/{inp.steps} E0={initial_energy:.6f} Efin={final_energy:.6f} "
            f"trace_len={len(trace)}"
        )
        return RelaxResult(
            initial_energy=initial_energy,
            final_energy=final_energy,
            positions=atoms.get_positions().tolist(),
            forces=forces.tolist(),
            stress=(stress.tolist() if stress is not None else None),
            steps_completed=steps_completed,
            calculator=inp.calculator,
            trace_energies=(trace if inp.return_trace else None),
        )
    except HTTPException:
        raise
    except Exception as e:
        tb = traceback.format_exc()
        print("[ERROR] /relax exception:\n", tb)
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/md", response_model=MDResult)
def md_step(inp: MDIn):
    """Run a basic NVT MD segment using ASE's Langevin integrator.

    Advances the provided coordinates by `steps` integration steps with
    timestep `timestep_fs` (fs) toward a target temperature using a Langevin
    thermostat (friction gamma = inp.friction 1/fs). No calculator fallback:
    the specified calculator is used directly.
    """
    try:
        import numpy as _np
        from ase import Atoms as _Atoms
        from ase import units as _units
        from ase.md.langevin import Langevin
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

        if inp.steps <= 0:
            raise HTTPException(status_code=400, detail="steps must be >0")

        # Build Atoms (reuse relax helper style)
        Z = inp.atomic_numbers
        if len(Z) != len(inp.coordinates):
            raise HTTPException(status_code=400, detail="Length mismatch atomic_numbers vs coordinates")
        pbc = (
            inp.pbc
            if inp.pbc is not None
            else ([True, True, True] if inp.cell else [False, False, False])
        )
        atoms = _Atoms(
            numbers=Z,
            positions=_np.array(inp.coordinates, dtype=float),
            cell=inp.cell,
            pbc=pbc,
        )
        atoms.info.update({"charge": inp.charge or 0, "spin": inp.spin_multiplicity or 1})

        # Attach the requested calculator directly
        ensure_model_loaded()
        if inp.calculator == RelaxCalculatorName.uma:
            atoms.calc = _CALCULATOR
        elif inp.calculator == RelaxCalculatorName.lj:
            if LennardJones is None:
                raise HTTPException(status_code=500, detail="ASE LJ unavailable")
            atoms.calc = LennardJones(rc=3.0)
        else:  # pragma: no cover
            raise HTTPException(status_code=400, detail="Unknown calculator")

        # Initialize velocities with Maxwell-Boltzmann distribution (Å/fs units internally)
        MaxwellBoltzmannDistribution(atoms, temperature_K=inp.temperature)

        # Baseline energy (any error surfaces directly)
        try:
            initial_energy = float(atoms.get_potential_energy())
        except Exception as ee:
            raise HTTPException(status_code=500, detail=f"Initial energy failed: {ee}")
        natoms = len(Z)
        mass = atoms.get_masses()  # amu (not directly used but kept for potential extensions)

        timestep = float(inp.timestep_fs) * _units.fs
        dyn = Langevin(atoms, timestep, temperature_K=inp.temperature, friction=inp.friction, logfile=None)

        energies = [] if inp.return_trajectory else None
        prev_positions = atoms.get_positions().copy()

        for _ in range(int(inp.steps)):
            dyn.run(1)
            if energies is not None:
                try:
                    energies.append(float(atoms.get_potential_energy()))
                except Exception:
                    energies.append(float('nan'))
            new_pos = atoms.get_positions()
            max_disp = _np.sqrt(((new_pos - prev_positions) ** 2).sum(axis=1)).max()
            if not _np.isfinite(max_disp) or max_disp > 5.0:
                raise HTTPException(status_code=500, detail=f"MD instability detected (max step disp {max_disp:.2f} Å)")
            prev_positions[:] = new_pos

        final_energy = float(atoms.get_potential_energy())
        final_forces = atoms.get_forces()
        final_vel = atoms.get_velocities()
        # Kinetic energy & temperature using ASE helper (atoms.get_kinetic_energy already in eV)
        KE = float(atoms.get_kinetic_energy())  # eV
        # Equipartition: KE = (3/2) N kB T  -> T = 2*KE / (3 N kB)
        final_T = (2.0 * KE) / (3.0 * natoms * _units.kB) if natoms > 0 else 0.0
        print(f"[md] calc={inp.calculator} natoms={natoms} steps={inp.steps} T_target={inp.temperature:.1f}K T_final={final_T:.1f}K E0={initial_energy:.6f} Efin={final_energy:.6f}")
        return MDResult(
            initial_energy=initial_energy,
            final_energy=final_energy,
            positions=atoms.get_positions().tolist(),
            velocities=final_vel.tolist(),
            forces=final_forces.tolist(),
            steps_completed=int(inp.steps),
            temperature=float(final_T),
            energies=energies,
            calculator=inp.calculator,
        )
    except HTTPException:
        raise
    except Exception as e:
        tb = traceback.format_exc()
        print("[ERROR] /md exception:\n", tb)
        raise HTTPException(status_code=500, detail=str(e))


# 1) env with PyTorch + fairchem-core + fastapi + uvicorn
#    See FAIRChem install docs (v2.x) & model names (e.g., uma-s-1p1).
#    You must have HF credentials locally if weights are gated:
#      huggingface-cli login
# 2) start
#    export UMA_MODEL=uma-s-1p1
#    export UMA_TASK=omol  # or omat/oc20/odac/omc depending on domain
#    uvicorn server:app --host 0.0.0.0 --port 8000
