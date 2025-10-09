import pytest

from fairchem_local_server.atoms_cache import cache_clear
from fairchem_local_server.models import (
    MDFromCacheIn,
    MDIn,
    RelaxCalculatorName,
    RelaxFromCacheIn,
    RelaxIn,
    SimpleFromCacheIn,
    SimpleIn,
)
from fairchem_local_server.services import (
    md_step,
    md_step_from_cache,
    relax,
    relax_from_cache,
    simple_calculate,
    simple_calculate_from_cache,
)


@pytest.fixture(autouse=True)
def _clear_cache_each():
    cache_clear()
    yield
    cache_clear()


WATER_Z = [8, 1, 1]
WATER_POS = [
    [0.0, 0.0, 0.0],
    [0.9575, 0.0, 0.0],
    [-0.2399872, 0.92662721, 0.0],
]


def test_simple_calculate_returns_cache_key():
    res = simple_calculate(
        SimpleIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            calculator=RelaxCalculatorName.lj,
            properties=["energy", "forces"],
        )
    )
    assert "cache_key" in res and isinstance(res["cache_key"], str)
    assert "results" in res and "energy" in res["results"]


def test_simple_calculate_from_cache_roundtrip():
    first = simple_calculate(
        SimpleIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            calculator=RelaxCalculatorName.lj,
            properties=["energy", "forces"],
        )
    )
    key = first["cache_key"]

    second = simple_calculate_from_cache(
        SimpleFromCacheIn(cache_key=key, calculator=RelaxCalculatorName.lj)
    )
    assert "cache_key" in second and isinstance(second["cache_key"], str)
    assert "energy" in second["results"]


def test_relax_from_cache_works():
    base = relax(
        RelaxIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
    )
    assert base.cache_key

    again = relax_from_cache(
        RelaxFromCacheIn(
            cache_key=base.cache_key,
            steps=1,
            calculator=base.calculator,
        )
    )
    assert again.steps_completed >= 1
    assert again.cache_key


def test_md_from_cache_works():
    base = md_step(
        MDIn(
            atomic_numbers=WATER_Z,
            coordinates=WATER_POS,
            steps=1,
            calculator=RelaxCalculatorName.lj,
        )
    )
    assert base.cache_key

    again = md_step_from_cache(
        MDFromCacheIn(
            cache_key=base.cache_key,
            steps=1,
            calculator=base.calculator,
            return_trajectory=True,
        )
    )
    assert again.steps_completed == 1
    assert again.cache_key
    assert isinstance(again.energies, list)


def test_invalid_cache_key_raises():
    with pytest.raises(Exception):
        relax_from_cache(
            RelaxFromCacheIn(
                cache_key="missing", steps=1, calculator=RelaxCalculatorName.lj
            )
        )
