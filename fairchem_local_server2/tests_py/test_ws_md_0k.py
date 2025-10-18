from __future__ import annotations

import asyncio
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest
import websockets
from fairchem_local_server2 import session_pb2 as pb


def _find_xyz() -> Path:
    """Prefer a small molecule for faster MD; fall back to ROY."""
    here = Path(__file__).resolve()
    candidates = [
        here.parent.parent.parent / "public" / "molecules" / "water.xyz",
        here.parent.parent.parent / "public" / "molecules" / "benzene.xyz",
        here.parent.parent.parent / "public" / "molecules" / "roy.xyz",
        Path("public/molecules/water.xyz"),
        Path("public/molecules/benzene.xyz"),
        Path("public/molecules/roy.xyz"),
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("Could not locate a default molecule XYZ file")


def _load_xyz(path: Path) -> Tuple[List[int], List[List[float]]]:
    try:
        from ase.io import read as ase_read  # type: ignore

        atoms = ase_read(str(path))  # type: ignore
        return (
            atoms.get_atomic_numbers().tolist(),
            atoms.get_positions().tolist(),
        )
    except Exception:
        txt = path.read_text().splitlines()
        nat = int(txt[0].strip())
        lines = txt[2 : 2 + nat]
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
        }
        Z: List[int] = []
        pos: List[List[float]] = []
        for ln in lines:
            s, x, y, z = ln.split()[:4]
            Z.append(PT[s])
            pos.append([float(x), float(y), float(z)])
        return Z, pos


async def _run_md_0k(
    uri: str, Z: List[int], xyz: List[List[float]], frames: int, calculator: str
) -> int:
    """Run MD at 0K for `frames` steps using protobuf-only messages."""
    async with websockets.connect(uri) as ws:
        seq = 1
        # INIT: include zero velocities to avoid MB initialization
        init = pb.ClientAction()
        init.seq = seq
        init.type = pb.ClientAction.Type.INIT_SYSTEM
        init.atomic_numbers.extend(int(z) for z in Z)
        for p in xyz:
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        for _ in Z:
            v = pb.Vec3()
            v.v.extend([0.0, 0.0, 0.0])
            init.velocities.append(v)
        await ws.send(init.SerializeToString())

        # START MD with temperature 0.0
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.type = pb.ClientAction.Type.START_SIMULATION
        start.simulation_type = pb.ClientAction.SimType.MD
        sp = pb.SimulationParams()
        sp.calculator = calculator
        sp.temperature = 0.0
        sp.timestep_fs = 1.0
        sp.friction = 0.02
        start.simulation_params.CopyFrom(sp)
        await ws.send(start.SerializeToString())

        last_seq = 0
        got = 0
        while got < int(frames):
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            assert isinstance(data, (bytes, bytearray))
            res = pb.ServerResult()
            res.ParseFromString(data)
            # basic sanity: lengths match and numbers are finite
            n = len(Z)
            assert len(res.positions) == n
            assert len(res.velocities) == n
            assert len(res.forces) == n
            for vecs in (res.positions, res.velocities, res.forces):
                for v in vecs:
                    assert len(v.v) == 3
                    assert np.isfinite(v.v).all()
            last_seq = int(res.seq)
            got += 1
            # ACK
            seq += 1
            ack = pb.ClientAction()
            ack.seq = seq
            ack.type = pb.ClientAction.Type.PING
            ack.ack = last_seq
            await ws.send(ack.SerializeToString())

        # STOP
        seq += 1
        stop = pb.ClientAction()
        stop.seq = seq
        stop.type = pb.ClientAction.Type.STOP_SIMULATION
        await ws.send(stop.SerializeToString())

        return got


@pytest.mark.timeout(60)
def test_ws_md_runs_at_0k(ws_base_url: str):
    # UMA-only server provided by session fixture
    xyz_path = _find_xyz()
    Z, xyz = _load_xyz(xyz_path)
    nframes = asyncio.run(
        _run_md_0k(ws_base_url, Z, xyz, frames=10, calculator="uma")
    )
    assert nframes >= 10
