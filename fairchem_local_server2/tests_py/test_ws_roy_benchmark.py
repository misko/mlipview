from __future__ import annotations

import asyncio
from pathlib import Path
from typing import List, Tuple

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


def _find_roy_xyz() -> Path:
    """Locate public/molecules/roy.xyz starting from repo root or test dir.

    Searches upward a few levels to be resilient to different cwd setups.
    """
    here = Path(__file__).resolve()
    for up in [here] + list(here.parents)[:6]:
        candidate = up.parent / "public" / "molecules" / "roy.xyz"
        if candidate.exists():
            return candidate
    # Fallback to repo-root relative when running from project root
    c2 = Path("public/molecules/roy.xyz")
    if c2.exists():
        return c2
    raise FileNotFoundError("public/molecules/roy.xyz not found")


def _load_xyz(path: Path) -> Tuple[List[int], List[List[float]]]:
    """Load an .xyz file and return (atomic_numbers, positions).

    Uses ASE if available; otherwise minimal parser for standard XYZ.
    """
    try:
        from ase.io import read as ase_read  # type: ignore

        atoms = ase_read(str(path))  # type: ignore
        return (
            atoms.get_atomic_numbers().tolist(),
            atoms.get_positions().tolist(),
        )
    except Exception:
        text = path.read_text().splitlines()
        if len(text) < 3:
            raise ValueError("Invalid XYZ file: too short")
        nat = int(text[0].strip())
        lines = text[2 : 2 + nat]
        PT = {
            "H": 1,
            "He": 2,
            "Li": 3,
            "Be": 4,
            "B": 5,
            "C": 6,
            "N": 7,
            "O": 8,
            "F": 9,
            "Ne": 10,
            "Na": 11,
            "Mg": 12,
            "Al": 13,
            "Si": 14,
            "P": 15,
            "S": 16,
            "Cl": 17,
            "Ar": 18,
            "K": 19,
            "Ca": 20,
            "Sc": 21,
            "Ti": 22,
            "V": 23,
            "Cr": 24,
            "Mn": 25,
            "Fe": 26,
            "Co": 27,
            "Ni": 28,
            "Cu": 29,
            "Zn": 30,
        }
        Z: List[int] = []
        pos: List[List[float]] = []
        for ln in lines:
            parts = ln.split()
            sym = parts[0]
            x, y, z = map(float, parts[1:4])
            Z.append(PT[sym])
            pos.append([x, y, z])
    return Z, pos


async def _run_ws_md_frames(
    uri: str,
    Z: List[int],
    xyz: List[List[float]],
    n: int,
    calculator: str,
) -> int:
    """
    Connect to websocket, initialize, start MD, ACK, and receive n frames.
    """
    frames = 0
    async with websockets.connect(uri) as ws:
        # Init via USER_INTERACTION (INIT_SYSTEM removed)
        seq = 1
        init = pb.ClientAction()
        init.seq = seq
        init.schema_version = 1
        ui = pb.UserInteractionSparse()
        ui.natoms = len(Z)
        inz = pb.IntDelta()
        inz.indices.extend(list(range(len(Z))))
        inz.values.extend(int(z) for z in Z)
        ui.atomic_numbers.CopyFrom(inz)
        vr = pb.Vec3Delta()
        vr.indices.extend(list(range(len(xyz))))
        for p in xyz:
            vr.coords.extend([float(p[0]), float(p[1]), float(p[2])])
        ui.positions.CopyFrom(vr)
        ui.full_update = True
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # START_SIMULATION (MD, LJ)
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.schema_version = 1
        st = pb.ClientAction.Start()
        st.simulation_type = pb.ClientAction.Start.SimType.MD
        sp = pb.SimulationParams()
        sp.calculator = calculator
        sp.temperature = 0.0
        sp.timestep_fs = 1.0
        sp.friction = 0.02
        st.simulation_params.CopyFrom(sp)
        start.start.CopyFrom(st)
        await ws.send(start.SerializeToString())

        last_seq = 0
        while frames < int(n):
            # Receive binary protobuf result
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            if not isinstance(data, (bytes, bytearray)):
                # Server is protobuf-only; ignore non-bytes
                continue
            res = pb.ServerResult()
            res.ParseFromString(data)
            last_seq = int(res.seq)
            frames += 1
            # Ack to allow server to continue beyond backpressure window
            seq += 1
            ack = pb.ClientAction()
            ack.seq = seq
            ack.ack = last_seq
            # ack-only; no ping payload in this schema
            await ws.send(ack.SerializeToString())

        # STOP_SIMULATION
        seq += 1
        stop = pb.ClientAction()
        stop.seq = seq
        stop.stop.CopyFrom(pb.ClientAction.Stop())
        await ws.send(stop.SerializeToString())

    return frames


@pytest.mark.timeout(60)
def test_ws_roy_md_30_frames(ws_base_url: str):
    xyz_path = _find_roy_xyz()
    Z, xyz = _load_xyz(xyz_path)
    frames = asyncio.run(_run_ws_md_frames(ws_base_url, Z, xyz, 30, "uma"))
    assert frames >= 30
