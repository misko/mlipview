"""Helpers for building ASE Atoms objects and extracting properties."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Sequence

import numpy as np
from fastapi import HTTPException


def center_and_return_shift(atoms, buffer=10.0):
    if atoms.cell is None or np.abs(atoms.cell.array).sum() < 1e-6:
        old_pos = atoms.get_positions().copy()[0][:]
        atoms.center(vacuum=buffer)
        return atoms.get_positions()[0] - old_pos
    return None


def build_atoms(
    numbers: Sequence[int],
    coords: Sequence[Sequence[float]],
    *,
    cell=None,
    charge=0,
    spin=1,
):
    import numpy as np
    from ase import Atoms

    # Lazy import to avoid any module import cycles
    try:
        from .models import MAX_ATOMS_PER_REQUEST  # type: ignore
    except Exception:
        MAX_ATOMS_PER_REQUEST = 170  # safe fallback

    if len(numbers) != len(coords):
        raise HTTPException(
            status_code=400,
            detail="Length mismatch atomic_numbers vs coordinates",
        )
    if len(numbers) > int(MAX_ATOMS_PER_REQUEST):
        raise HTTPException(
            status_code=400,
            detail=(
                "Too many atoms: "
                + f"{len(numbers)} > MAX_ATOMS_PER_REQUEST="
                + f"{MAX_ATOMS_PER_REQUEST}"
            ),
        )
    atoms = Atoms(
        numbers=numbers,
        positions=np.array(coords, dtype=float),
        cell=cell,
        pbc=[True, True, True],
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
