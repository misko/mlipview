import json
import os
from typing import Any, Optional, Sequence

import requests
from ase.calculators.calculator import CalculationFailed, Calculator, all_changes
from ase.io.jsonio import encode

FAIRCHEM_SERVER_URL = os.environ.get("FAIRCHEM_SERVER_URL", "http://127.0.0.1:8000")


class FairChemHTTP(Calculator):
    implemented_properties = ["energy", "forces", "stress", "free_energy"]

    def __init__(self, url: Optional[str] = None, timeout: float = 120.0, **kwargs):
        super().__init__(**kwargs)
        self.url = url or FAIRCHEM_SERVER_URL
        self.timeout = timeout

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):  # type: ignore[override]
        Calculator.calculate(self, atoms, properties, system_changes)
        if atoms is None:
            raise CalculationFailed("No atoms provided")
        # Prepare atoms JSON via ase.io.jsonio.encode
        encoded = encode(atoms)
        if isinstance(encoded, str):
            atoms_json = json.loads(encoded)
        else:
            atoms_json = encoded
        payload = {
            "atoms_json": atoms_json,
            "properties": list(properties),
            "info": atoms.info or {},
        }
        try:
            resp = requests.post(
                f"{self.url}/calculate", json=payload, timeout=self.timeout
            )
        except requests.RequestException as e:  # pragma: no cover
            raise CalculationFailed(f"HTTP request failed: {e}")
        if resp.status_code != 200:
            raise CalculationFailed(f"Server error {resp.status_code}: {resp.text}")
        data = resp.json()
        results = data.get("results", {})
        energy = results.get("energy")
        if energy is not None:
            self.results["energy"] = float(energy)
            self.results["free_energy"] = float(energy)
        if "forces" in results and results["forces"] is not None:
            import numpy as np

            self.results["forces"] = np.array(results["forces"], dtype=float)
        if "stress" in results and results["stress"] is not None:
            import numpy as np

            self.results["stress"] = np.array(results["stress"], dtype=float)
        return self.results
