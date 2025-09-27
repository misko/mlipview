# uma_http_calculator.py
import json

import numpy as np
import requests
from ase.calculators.calculator import Calculator
from ase.io.jsonio import decode, encode


class UMAHTTPCalculator(Calculator):
    implemented_properties = ["energy", "free_energy", "forces", "stress"]

    def __init__(self, url: str, task_name: str = "omol", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.url = url
        self.task_name = task_name

    def calculate(
        self,
        atoms=None,
        properties=("energy", "forces"),
        system_changes=("positions",),
    ):
        Calculator.calculate(self, atoms, properties, system_changes)

        # encode() returns a JSON string; convert to dict for FastAPI model
        encoded = encode(atoms)
        if isinstance(encoded, str):
            atoms_json_obj = json.loads(encoded)
        else:  # older versions might already give a dict
            atoms_json_obj = encoded

        payload = {
            # dict matching CalcIn.atoms_json schema
            "atoms_json": atoms_json_obj,
            "properties": list(properties),
            # optional metadata like charge/spin
            "info": getattr(atoms, "info", {}),
        }
        r = requests.post(
            self.url.rstrip("/") + "/calculate",
            json=payload,
            timeout=300,
        )
        r.raise_for_status()
        out = r.json()["results"]

        # Pack into ASE-style results dict
        res = {}
        if "energy" in out and out["energy"] is not None:
            res["energy"] = out["energy"]
            res["free_energy"] = out.get("free_energy", out["energy"])
        if "forces" in out and out["forces"] is not None:
            res["forces"] = np.array(out["forces"])
        if "stress" in out:
            if out["stress"] is not None:
                res["stress"] = np.array(out["stress"])
        self.results = res
