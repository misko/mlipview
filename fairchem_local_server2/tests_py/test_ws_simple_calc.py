"""
End-to-end (WS): idle compute returns finite energy/forces.

What this test covers
- Protocol: protobuf over WebSocket with oneof payloads, flat arrays for
    vectors, row‑major Mat3, schema_version=1.
- Flow:
    1) Send USER_INTERACTION to initialize a small system (atoms + positions,
         optional cell).
    2) Send another USER_INTERACTION with positions only to trigger an idle
         compute (no simulation running).
    3) Receive ServerResult frames until a Frame carries energy or forces, then
         validate the payload.

Expectations
- If energy is present, it should be a finite number (calculator-dependent).
- Forces must be present, have shape (N, 3), and be finite.

Why this exists
- Ensures the server’s idle compute path works with the new protobuf schema
    (flat arrays and oneof), without requiring a simulation loop.
"""

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
        init.schema_version = 1
        ui = pb.UserInteractionSparse()
        ui.natoms = len(Z)
        inz = pb.IntDelta()
        inz.indices.extend(list(range(len(Z))))
        inz.values.extend(int(z) for z in Z)
        ui.atomic_numbers.CopyFrom(inz)
        vr = pb.Vec3Delta()
        vr.indices.extend(list(range(len(R))))
        for p in R:
            vr.coords.extend([float(p[0]), float(p[1]), float(p[2])])
        ui.positions.CopyFrom(vr)
        if cell is not None:
            m = (
                pb.UserInteractionSparse.Mat3()
                if hasattr(pb, "UserInteractionSparse")
                else pb.Mat3()
            )
            m.m.extend(list(map(float, cell)))
            ui.cell.CopyFrom(m)
        ui.full_update = True
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # Trigger idle compute by sending USER_INTERACTION with positions only
        seq += 1
        req = pb.ClientAction()
        req.seq = seq
        req.schema_version = 1
        ui2 = pb.UserInteractionSparse()
        vr2 = pb.Vec3Delta()
        vr2.indices.extend(list(range(len(R))))
        for p in R:
            vr2.coords.extend([float(p[0]), float(p[1]), float(p[2])])
        ui2.positions.CopyFrom(vr2)
        ui2.full_update = False
        req.user_interaction.CopyFrom(ui2)
        # Attach optional counters if fields exist
        if hasattr(req, "user_interaction_count"):
            setattr(req, "user_interaction_count", 0)
        if hasattr(req, "sim_step"):
            setattr(req, "sim_step", 0)
        await ws.send(req.SerializeToString())

        # Receive frames until a ServerResult carries energy or forces
        for _ in range(10):
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            assert isinstance(data, (bytes, bytearray))
            res = pb.ServerResult()
            res.ParseFromString(data)
            if res.WhichOneof("payload") != "frame":
                continue
            fr = res.frame
            has_energy_field = hasattr(fr, "energy")
            energy_present = (
                hasattr(fr, "HasField") and fr.HasField("energy")
            ) or getattr(fr, "energy", None) is not None
            got_energy = has_energy_field and energy_present
            if got_energy or (len(fr.forces) > 0):
                energy = getattr(fr, "energy", None)
                f = np.fromiter(fr.forces, dtype=np.float64)
                if f.size % 3 == 0:
                    F = f.reshape(-1, 3)
                else:
                    F = np.empty((0, 3))
                return {"energy": energy, "forces": F}
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
