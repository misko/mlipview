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
    uri: str,
    Z: List[int],
    xyz: List[List[float]],
    frames: int,
    calculator: str,
) -> int:
    """Run MD at 0K for `frames` steps using protobuf-only messages."""
    async with websockets.connect(uri) as ws:
        seq = 1
        # INIT: include zero velocities to avoid MB initialization
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
        vv = pb.Vec3Delta()
        vv.indices.extend(list(range(len(Z))))
        for _ in Z:
            vv.coords.extend([0.0, 0.0, 0.0])
        ui.velocities.CopyFrom(vv)
        ui.full_update = True
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # START MD with temperature 0.0
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
        got = 0
        while got < int(frames):
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            assert isinstance(data, (bytes, bytearray))
            res = pb.ServerResult()
            res.ParseFromString(data)
            if res.WhichOneof("payload") != "frame":
                continue
            fr = res.frame
            if len(fr.positions) == 0:
                continue
            # basic sanity: lengths match and numbers are finite
            n = len(Z)
            P = np.fromiter(fr.positions, dtype=np.float64)
            V = np.fromiter(fr.velocities, dtype=np.float64)
            F = np.fromiter(fr.forces, dtype=np.float64)
            assert P.size % 3 == 0 and V.size % 3 == 0 and F.size % 3 == 0
            assert P.reshape(-1, 3).shape[0] == n
            assert V.reshape(-1, 3).shape[0] == n
            assert F.reshape(-1, 3).shape[0] == n
            last_seq = int(res.seq)
            got += 1
            # ACK
            seq += 1
            ack = pb.ClientAction()
            ack.seq = seq
            ack.ack = last_seq
            # ack-only; no ping payload in this schema
            await ws.send(ack.SerializeToString())

        # STOP
        seq += 1
        stop = pb.ClientAction()
        stop.seq = seq
        stop.stop.CopyFrom(pb.ClientAction.Stop())
        await ws.send(stop.SerializeToString())

        return got


@pytest.mark.timeout(60)
def test_ws_md_runs_at_0k(ws_base_url: str):
    # UMA-only server provided by session fixture
    xyz_path = _find_xyz()
    Z, xyz = _load_xyz(xyz_path)
    nframes = asyncio.run(_run_md_0k(ws_base_url, Z, xyz, frames=10, calculator="uma"))
    assert nframes >= 10
