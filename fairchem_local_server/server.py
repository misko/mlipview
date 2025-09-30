# server.py
import os
import threading
import traceback
from typing import Any, Dict, List, Optional, Tuple

import torch
from ase.io.jsonio import decode, encode
from fairchem.core import pretrained_mlip  # fairchem-core >= 2.x
from fairchem.core import FAIRChemCalculator
from fairchem.core.units.mlip_unit import InferenceSettings
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# --- configure model ---
MODEL_NAME = os.environ.get("UMA_MODEL", "uma-s-1p1")  # UMA small
TASK_NAME = os.environ.get("UMA_TASK", "omol")  # task: omol|omat|oc20|odac|omc
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Load a default predict unit at startup (fallback / first use)
_default_predict_unit = pretrained_mlip.get_predict_unit(
    MODEL_NAME,
    device=DEVICE,
    inference_settings=InferenceSettings(
        tf32=True,
        activation_checkpointing=False,
        merge_mole=True,
        compile=False,
        external_graph_gen=False,
        internal_graph_gen_version=2,
    ),
)
_default_calc = FAIRChemCalculator(_default_predict_unit, task_name=TASK_NAME)

# ----------------------------------------------------------------------------
# Composition-based predict unit + calculator cache
# ----------------------------------------------------------------------------
# Some UMA model internal graph building / compilation work can benefit from
# reusing a predict unit + calculator for identical compositions (multiset of
# atomic numbers). We store a cache keyed by a hashable composition signature.

_cache_lock = threading.RLock()
_composition_cache: Dict[
    Tuple[Tuple[int, int], ...], Tuple[Any, FAIRChemCalculator]
] = {}


def _composition_key(atomic_numbers: List[int]) -> Tuple[Tuple[int, int], ...]:
    """Return a canonical hashable key for a list of atomic numbers.

    Key is a tuple of (Z, count) pairs sorted by Z. Order-invariant.
    """
    from collections import Counter

    c = Counter(int(z) for z in atomic_numbers)
    return tuple(sorted(c.items()))


def get_cached_unit_and_calculator(atomic_numbers: List[int]):
    """Return (predict_unit, calculator) for the given atomic numbers.

    Reuses cached objects for identical compositions to avoid repeated model
    instantiation / graph setup cost. Falls back to a default global instance
    when atomic_numbers is empty or too small to justify a separate entry.
    """
    if not atomic_numbers:
        return _default_predict_unit, _default_calc
    key = _composition_key(atomic_numbers)
    with _cache_lock:
        hit = _composition_cache.get(key)
        if hit is not None:
            return hit
        # Miss: create a new predict unit + calculator.
        # NOTE: We reuse the same MODEL_NAME / TASK_NAME but could extend this
        # to adapt hyperparameters based on composition size if needed.
        pu = pretrained_mlip.get_predict_unit(
            MODEL_NAME,
            device=DEVICE,
            inference_settings=InferenceSettings(
                tf32=True,
                activation_checkpointing=False,
                merge_mole=True,
                compile=False,
                external_graph_gen=False,
                internal_graph_gen_version=2,
            ),
        )
        calculator = FAIRChemCalculator(pu, task_name=TASK_NAME)
        _composition_cache[key] = (pu, calculator)
        print(
            "[cache] created new predict_unit for composition="
            f"{key} cache_size={len(_composition_cache)}"
        )
        return pu, calculator


app = FastAPI(title="UMA-small ASE HTTP server")
# Re-export encode for test modules
__all__ = ["app", "encode", "get_cached_unit_and_calculator"]

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
        # Use / create cached calculator based on composition
        _pu, _calc = get_cached_unit_and_calculator(Z)
        atoms.calc = _calc
        results = {}
        # --- compute selected properties ---
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
            Z = list(
                getattr(
                    atoms,
                    "numbers",
                    getattr(atoms, "get_atomic_numbers")(),
                )
            )
        except Exception:  # fallback if Atoms style accessor fails
            Z = []
        _pu, _calc = get_cached_unit_and_calculator(Z)
        atoms.calc = _calc
        results = {}
        if "energy" in props or "free_energy" in props:
            results["energy"] = float(atoms.get_potential_energy())
            # some ASE front-ends expect free_energy = energy for MLIPs
            results["free_energy"] = results["energy"]
        if "forces" in props:
            results["forces"] = atoms.get_forces().tolist()
        if "stress" in props:
            try:
                results["stress"] = atoms.get_stress().tolist()
            except Exception:
                # stress may be unavailable for some tasks/geometries
                results["stress"] = None
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


# 1) env with PyTorch + fairchem-core + fastapi + uvicorn
#    See FAIRChem install docs (v2.x) & model names (e.g., uma-s-1p1).
#    You must have HF credentials locally if weights are gated:
#      huggingface-cli login
# 2) start
#    export UMA_MODEL=uma-s-1p1
#    export UMA_TASK=omol  # or omat/oc20/odac/omc depending on domain
#    uvicorn server:app --host 0.0.0.0 --port 8000
