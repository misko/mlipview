from __future__ import annotations

import asyncio
from typing import List

import numpy as np
import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _ws_simple_once(
    uri: str, Z: List[int], R: List[List[float]], cell=None
) -> dict:
    async with websockets.connect(uri) as ws:
        seq = 1
        # Init via USER_INTERACTION with atoms + positions (and optional cell)
        init = pb.ClientAction()
        init.seq = seq
        init.type = pb.ClientAction.Type.USER_INTERACTION
        init.atomic_numbers.extend(int(z) for z in Z)
        for p in R:
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        if cell is not None:
            m = pb.Mat3()
            m.m.extend(list(map(float, cell)))
            init.cell.CopyFrom(m)
        await ws.send(init.SerializeToString())

        # Trigger idle compute by sending USER_INTERACTION with positions only
        seq += 1
        req = pb.ClientAction()
        req.seq = seq
        req.type = pb.ClientAction.Type.USER_INTERACTION
        for p in R:
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            req.positions.append(v)
        # Set counters; new schema guarantees these fields
        req.user_interaction_count = 0
        req.sim_step = 0
        await ws.send(req.SerializeToString())

        # Receive frames until a ServerResult carries energy or forces
        for _ in range(10):
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            # Protocol is protobuf-only; enforce binary frames
            assert isinstance(data, (bytes, bytearray))
            res = pb.ServerResult()
            res.ParseFromString(data)
            # Accept any frame that includes at least energy or forces
            has_energy = res.HasField("energy")
            if has_energy or (len(res.forces) > 0):
                energy = getattr(res, "energy", None)
                forces = np.array(
                    [[v.v[0], v.v[1], v.v[2]] for v in res.forces]
                )
                return {"energy": energy, "forces": forces}
    raise AssertionError("No idle compute result received")


@pytest.mark.timeout(60)
def test_ws_idle_compute_returns_finite(ws_base_url: str):
    # simple H2-like system
    Z = [1, 1]
    R = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]
    res = asyncio.run(_ws_simple_once(ws_base_url, Z, R))
    assert res is not None
    # Energy may vary by calculator; just assert finite number if present
    if res.get("energy") is not None:
        assert np.isfinite(float(res["energy"]))
    F = res.get("forces")
    assert F is not None
    assert F.shape == (2, 3)
    assert np.isfinite(F).all()
