import math
from typing import List

import pytest

from fairchem_local_server.models import (
    MAX_ATOMS_PER_REQUEST,
    MDIn,
    RelaxCalculatorName,
    RelaxIn,
    SimpleIn,
)
from fairchem_local_server.services import md_step, relax, simple_calculate


def _dummy_coords(n: int) -> List[List[float]]:
    # Simple cubic lattice to avoid overlaps; 2.0 Å spacing
    coords = []
    for i in range(n):
        x = (i % 4) * 2.0
        y = ((i // 4) % 4) * 2.0
        z = (i // 16) * 2.0
        coords.append([float(x), float(y), float(z)])
    return coords


@pytest.mark.parametrize("n", [MAX_ATOMS_PER_REQUEST])
def test_models_accept_at_limit(n):
    zs = [1] * n
    coords = _dummy_coords(n)
    # Construction should succeed
    s = SimpleIn(atomic_numbers=zs, coordinates=coords)
    r = RelaxIn(atomic_numbers=zs, coordinates=coords, steps=1)
    m = MDIn(atomic_numbers=zs, coordinates=coords, steps=1)
    assert len(s.atomic_numbers) == n
    assert len(r.atomic_numbers) == n
    assert len(m.atomic_numbers) == n


@pytest.mark.parametrize("n", [MAX_ATOMS_PER_REQUEST + 1])
def test_models_reject_over_limit(n):
    zs = [1] * n
    coords = _dummy_coords(n)
    with pytest.raises(ValueError):
        SimpleIn(atomic_numbers=zs, coordinates=coords)
    with pytest.raises(ValueError):
        RelaxIn(atomic_numbers=zs, coordinates=coords, steps=1)
    with pytest.raises(ValueError):
        MDIn(atomic_numbers=zs, coordinates=coords, steps=1)


def test_services_accept_at_limit_cpu_lj():
    n = MAX_ATOMS_PER_REQUEST
    zs = [1] * n
    coords = _dummy_coords(n)

    # Use LJ calculator to avoid UMA/GPU dependency
    s = SimpleIn(
        atomic_numbers=zs,
        coordinates=coords,
        properties=["energy", "forces"],
        calculator=RelaxCalculatorName.lj,
    )
    out = simple_calculate(s)
    assert "results" in out
    assert set(out["results"].keys()) >= {"energy", "forces"}

    r = RelaxIn(
        atomic_numbers=zs,
        coordinates=coords,
        steps=2,
        calculator=RelaxCalculatorName.lj,
    )
    rout = relax(r)
    assert rout.steps_completed >= 1
    assert math.isfinite(rout.final_energy)

    m = MDIn(
        atomic_numbers=zs,
        coordinates=coords,
        steps=2,
        calculator=RelaxCalculatorName.lj,
        temperature=10.0,
        timestep_fs=1.0,
        friction=0.05,
    )
    mout = md_step(m)
    assert mout.steps_completed == 2
    assert len(mout.positions) == n


def test_services_reject_over_limit_via_model():
    n = MAX_ATOMS_PER_REQUEST + 1
    zs = [1] * n
    coords = _dummy_coords(n)
    # Model construction should raise before service is invoked
    with pytest.raises(ValueError):
        SimpleIn(
            atomic_numbers=zs,
            coordinates=coords,
            calculator=RelaxCalculatorName.lj,
        )
    with pytest.raises(ValueError):
        RelaxIn(
            atomic_numbers=zs,
            coordinates=coords,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
    with pytest.raises(ValueError):
        MDIn(
            atomic_numbers=zs,
            coordinates=coords,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
