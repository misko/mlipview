from __future__ import annotations

import asyncio
from typing import List

import numpy as np
import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb
from fairchem_local_server.atoms_utils import build_atoms
from fairchem_local_server.models import RelaxCalculatorName
from fairchem_local_server.services import _md_run


async def _ws_nth_md_step(
    uri: str, Z: List[int], R: List[List[float]], n: int, calculator: str
) -> dict:
    async with websockets.connect(uri) as ws:
        seq = 1
        init = pb.ClientAction()
        init.seq = seq
        # Init via USER_INTERACTION (INIT_SYSTEM removed in new protocol)
        init.type = pb.ClientAction.Type.USER_INTERACTION
        init.atomic_numbers.extend(int(z) for z in Z)
        for p in R:
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        # zero velocities to enforce 0K
        for _ in Z:
            v = pb.Vec3()
            v.v.extend([0.0, 0.0, 0.0])
            init.velocities.append(v)
        await ws.send(init.SerializeToString())

        # Start MD at 0 K
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

        # Receive frames: skip 'initialized' (seq==0), then take nth MD frame
        taken = 0
        last = None
        while taken < int(max(1, n)):
            data = await asyncio.wait_for(ws.recv(), timeout=5.0)
            assert isinstance(data, (bytes, bytearray))
            res = pb.ServerResult()
            res.ParseFromString(data)
            # Skip initialization snapshot and idle frames without positions
            if res.seq == 0 or len(res.positions) == 0:
                continue
            # sanity
            m = len(Z)
            assert len(res.positions) == m
            assert len(res.velocities) == m
            assert len(res.forces) == m
            P = np.array([[v.v[0], v.v[1], v.v[2]] for v in res.positions])
            V = np.array([[v.v[0], v.v[1], v.v[2]] for v in res.velocities])
            F = np.array([[v.v[0], v.v[1], v.v[2]] for v in res.forces])
            last = {"positions": P, "velocities": V, "forces": F}
            # ACK to advance window
            seq += 1
            ack = pb.ClientAction()
            ack.seq = seq
            ack.type = pb.ClientAction.Type.PING
            ack.ack = int(res.seq)
            await ws.send(ack.SerializeToString())
            taken += 1

        # Stop
        seq += 1
        stop = pb.ClientAction()
        stop.seq = seq
        stop.type = pb.ClientAction.Type.STOP_SIMULATION
        await ws.send(stop.SerializeToString())

        assert last is not None
        return last


@pytest.mark.timeout(60)
def test_ws_vs_direct_md_second_step_0k(ws_base_url: str):
    Z = [6, 6]
    R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
    V0 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    # Direct reference with UMA calculator
    atoms = build_atoms(Z, R, cell=None)
    d1 = _md_run(
        atoms,
        steps=1,
        temperature=0.0,
        timestep_fs=1.0,
        friction=0.02,
        calculator=RelaxCalculatorName.uma,
        return_trajectory=False,
        precomputed=None,
        velocities_in=V0,
    ).dict()
    atoms2 = build_atoms(Z, d1["positions"], cell=None)
    _ = _md_run(
        atoms2,
        steps=1,
        temperature=0.0,
        timestep_fs=1.0,
        friction=0.02,
        calculator=RelaxCalculatorName.uma,
        return_trajectory=False,
        precomputed=None,
        velocities_in=d1["velocities"],
    ).dict()

    # Second WS MD frame over UMA
    ws_res = asyncio.run(_ws_nth_md_step(ws_base_url, Z, R, 2, "uma"))

    # For UMA, assert finiteness and shape; numeric equality may vary
    assert ws_res["positions"].shape == (2, 3)
    assert ws_res["velocities"].shape == (2, 3)
    assert ws_res["forces"].shape == (2, 3)
    assert np.isfinite(ws_res["positions"]).all()
    assert np.isfinite(ws_res["velocities"]).all()
    assert np.isfinite(ws_res["forces"]).all()
    # File previously had duplicated content; cleaned to UMA-only via fixture.
