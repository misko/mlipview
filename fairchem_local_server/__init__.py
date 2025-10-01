"""Package init for fairchem_local_server.

Exposes key server objects for tests and external imports.
"""

from .server import app, encode, get_cached_unit_and_calculator, RelaxCalculatorName  # noqa: F401
