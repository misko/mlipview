"""Helpers for building ASE Atoms objects and extracting properties."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Sequence

from fastapi import HTTPException


def build_atoms(
    numbers: Sequence[int],
    coords: Sequence[Sequence[float]],
    *,
    cell=None,
    pbc=None,
    charge=0,
    spin=1,
):
    import numpy as np
    from ase import Atoms

    if len(numbers) != len(coords):
        raise HTTPException(
            status_code=400,
            detail="Length mismatch atomic_numbers vs coordinates",
        )
    if pbc is None:
        pbc = [True, True, True] if cell is not None else [False, False, False]
    atoms = Atoms(
        numbers=numbers,
        positions=np.array(coords, dtype=float),
        cell=cell,
        pbc=pbc,
    )
    atoms.info.update({"charge": charge or 0, "spin": spin or 1})
    return atoms


def compute_properties(atoms, props: Iterable[str]):
    res: Dict[str, Any] = {}
    props_tuple = tuple(props)
    if "energy" in props_tuple or "free_energy" in props_tuple:
        e = float(atoms.get_potential_energy())
        res["energy"] = e
        res["free_energy"] = e
    if "forces" in props_tuple:
        res["forces"] = atoms.get_forces().tolist()
    if "stress" in props_tuple:
        try:
            res["stress"] = atoms.get_stress().tolist()
        except Exception:
            res["stress"] = None
    return res


__all__ = ["build_atoms", "compute_properties"]
