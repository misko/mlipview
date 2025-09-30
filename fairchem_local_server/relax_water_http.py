import json
from pathlib import Path

from ase.build import molecule
from ase.optimize import BFGS
from http_calculator import FairChemHTTP

OUT = Path(__file__).parent / "water_fairchem_bfgs_trace.json"


def main():
    atoms = molecule("H2O")
    atoms.info.update({"charge": 0, "spin": 1})
    calc = FairChemHTTP()
    atoms.calc = calc
    trace = []
    # record initial
    e0 = atoms.get_potential_energy()
    trace.append(e0)
    opt = BFGS(atoms, logfile=None)
    # Run fixed number of steps (20) like JS parity
    for _ in range(20):
        opt.step()
        trace.append(atoms.get_potential_energy())
    data = {"energies": trace}
    OUT.write_text(json.dumps(data, indent=2))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
