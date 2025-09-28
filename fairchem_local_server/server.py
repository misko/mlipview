# server.py
import os
import traceback
from typing import Any, Dict, List, Optional

import torch
from ase.io.jsonio import decode, encode
from fairchem.core import FAIRChemCalculator, pretrained_mlip  # fairchem-core >= 2.x
from fairchem.core.units.mlip_unit import InferenceSettings
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# --- configure model ---
MODEL_NAME = os.environ.get("UMA_MODEL", "uma-s-1p1")  # UMA small
TASK_NAME = os.environ.get("UMA_TASK", "omol")  # omol / omat / oc20 / odac / omc
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Load once at startup (kept in memory)
predictor = pretrained_mlip.get_predict_unit(
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
calc = FAIRChemCalculator(predictor, task_name=TASK_NAME)

app = FastAPI(title="UMA-small ASE HTTP server")

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
        atoms = _Atoms(numbers=Z, positions=_np.array(xyz, dtype=float))
        atoms.info.update(
            {
                "charge": inp.charge or 0,
                "spin": inp.spin_multiplicity or 1,
            }
        )
        props = tuple(inp.properties or ("energy", "forces"))
        atoms.calc = calc
        results = {}
        # --- compute selected properties ---
        if "energy" in props or "free_energy" in props:
            results["energy"] = float(atoms.get_potential_energy())
            results["free_energy"] = results["energy"]
        if "forces" in props:
            results["forces"] = atoms.get_forces().tolist()
        dt = (_time.time() - t0) * 1000.0
        per_atom = dt / max(1, len(Z))
        print(
            "[simple_calculate] natoms="
            f"{len(Z)} props={props} time_ms={dt:.1f} "
            f"ms_per_atom={per_atom:.2f}"
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
        atoms.calc = calc
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
        print(
            "[calculate] natoms="
            f"{natoms} props={props} time_ms={dt:.1f} "
            f"ms_per_atom={per_atom:.2f}"
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
