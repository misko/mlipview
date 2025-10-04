import math
from fairchem_local_server.models import (
    RelaxIn,
    MDIn,
    RelaxCalculatorName,
    PrecomputedValues,
)
from fairchem_local_server.services import relax, md_step

# Lightweight tests using LJ calculator so no UMA / Ray dependency

WATER_Z = [8, 1, 1]
WATER_POS = [
    [0.0, 0.0, 0.0],
    [0.9575, 0.0, 0.0],
    [-0.2399872, 0.92662721, 0.0],
]


def test_relax_precomputed_energy_single_step():
    base = relax(
        RelaxIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
    )
    injected = base.initial_energy + 0.123456
    res = relax(
        RelaxIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
            precomputed=PrecomputedValues(energy=injected),
        )
    )
    assert math.isclose(res.initial_energy, injected, rel_tol=0, abs_tol=1e-10)
    assert res.precomputed_applied == ["energy"]


def test_md_precomputed_energy_single_step():
    base = md_step(
        MDIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
    )
    injected = base.initial_energy + 0.222222
    res = md_step(
        MDIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
            precomputed=PrecomputedValues(energy=injected),
        )
    )
    assert math.isclose(res.initial_energy, injected, rel_tol=0, abs_tol=1e-10)
    assert res.precomputed_applied == ["energy"]


def test_md_precomputed_forces_only():
    zero_forces = [[0.0, 0.0, 0.0] for _ in WATER_Z]
    res = md_step(
        MDIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
            precomputed=PrecomputedValues(forces=zero_forces),
        )
    )
    assert res.precomputed_applied == ["forces"]
    assert len(res.forces) == 3


def test_relax_precomputed_stress_voigt():
    # Provide an arbitrary Voigt stress vector of length 6
    stress = [0.1, 0.2, 0.3, 0.01, 0.02, 0.03]
    res = relax(
        RelaxIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
            precomputed=PrecomputedValues(stress=stress),
        )
    )
    assert res.precomputed_applied == ["stress"]


def test_relax_precomputed_stress_matrix():
    # 3x3 row-major flattened matrix -> convert to Voigt
    matrix = [
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
    ]
    res = relax(
        RelaxIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
            precomputed=PrecomputedValues(stress=matrix),
        )
    )
    # Expect Voigt: xx=1, yy=5, zz=9, yz=6, xz=3, xy=2
    assert res.precomputed_applied == ["stress"]
