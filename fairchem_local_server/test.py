# run_client.py
from ase.build import molecule
from uma_http_calculator import UMAHTTPCalculator

atoms = molecule("H2O")
atoms.info.update({"charge": 0, "spin": 1})  # omol accepts charge/spin

atoms.calc = UMAHTTPCalculator("http://localhost:8000", task_name="omol")
print("E (eV):", atoms.get_potential_energy())
print("F (eV/Å):\n", atoms.get_forces())
