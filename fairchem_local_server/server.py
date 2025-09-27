# server.py
import os
import traceback
from typing import Any, Dict, List, Optional

import torch
from ase.io.jsonio import decode, encode
from fairchem.core import FAIRChemCalculator, pretrained_mlip  # fairchem-core >= 2.x
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

# --- configure model ---
MODEL_NAME = os.environ.get("UMA_MODEL", "uma-s-1p1")  # UMA small
TASK_NAME = os.environ.get("UMA_TASK", "omol")  # omol / omat / oc20 / odac / omc
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Load once at startup (kept in memory)
predictor = pretrained_mlip.get_predict_unit(MODEL_NAME, device=DEVICE)
calc = FAIRChemCalculator(predictor, task_name=TASK_NAME)

app = FastAPI(title="UMA-small ASE HTTP server")


class CalcIn(BaseModel):
    atoms_json: Dict[str, Any]  # result of ase.io.jsonio.encode(Atoms)
    properties: Optional[List[str]] = None  # e.g. ["energy","forces","stress"]
    info: Optional[Dict[str, Any]] = None  # optional atoms.info overrides


@app.post("/calculate")
def calculate(inp: CalcIn):
    try:
        # Reconstruct Atoms (accept either JSON string or pre-parsed dict)
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


## 1) env with PyTorch + fairchem-core + fastapi + uvicorn
#    See FAIRChem install docs for v2.x (pip/uv) and note model names (e.g., uma-s-1p1).
#    You must have HF credentials locally if the weights are gated.
#    huggingface-cli login
# 2) start
# export UMA_MODEL=uma-s-1p1
# export UMA_TASK=omol    # or omat/oc20/odac/omc depending on your domain
# uvicorn server:app --host 0.0.0.0 --port 8000
