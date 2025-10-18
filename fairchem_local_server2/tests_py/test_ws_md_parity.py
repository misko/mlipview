from __future__ import annotations

import asyncio
from typing import List

import numpy as np
import pytest
import websockets
from ray import serve

from fairchem_local_server.atoms_utils import build_atoms
from fairchem_local_server.models import RelaxCalculatorName
from fairchem_local_server.services import _md_run
from fairchem_local_server2.ws_app import deploy
from fairchem_local_server2 import session_pb2 as pb


async def _ws_nth_md_step(
    uri: str, Z: List[int], R: List[List[float]], n: int
) -> dict:
    async with websockets.connect(uri) as ws:
        seq = 1
        init = pb.ClientAction()
        init.seq = seq
        init.type = pb.ClientAction.Type.INIT_SYSTEM
        init.atomic_numbers.extend(int(z) for z in Z)
        for p in R:
            v = pb.Vec3(); v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        # zero velocities to enforce 0K
        for _ in Z:
            v = pb.Vec3(); v.v.extend([0.0, 0.0, 0.0])
            init.velocities.append(v)
        await ws.send(init.SerializeToString())

        # Start MD at 0 K
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.type = pb.ClientAction.Type.START_SIMULATION
        start.simulation_type = pb.ClientAction.SimType.MD
        sp = pb.SimulationParams()
        sp.calculator = "lj"
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
            res = pb.ServerResult(); res.ParseFromString(data)
            if res.seq == 0:
                continue  # initialization snapshot
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
            ack = pb.ClientAction(); ack.seq = seq
            ack.type = pb.ClientAction.Type.PING
            ack.ack = int(res.seq)
            await ws.send(ack.SerializeToString())
            taken += 1

        # Stop
        seq += 1
        stop = pb.ClientAction(); stop.seq = seq
        stop.type = pb.ClientAction.Type.STOP_SIMULATION
        await ws.send(stop.SerializeToString())

        assert last is not None
        return last


@pytest.mark.timeout(60)
def test_ws_vs_direct_md_second_step_0k():
    # Local LJ-only deployment
    deploy(ngpus=0, ncpus=1, nhttp=1)
    try:
        Z = [6, 6]
        R = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
        V0 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        # Build a direct reference that mirrors streaming semantics:
        # run two sequential 1-step MD calls, feeding back positions+velocities
        atoms = build_atoms(Z, R, cell=None)
        d1 = _md_run(
            atoms,
            steps=1,
            temperature=0.0,
            timestep_fs=1.0,
            friction=0.02,
            calculator=RelaxCalculatorName.lj,
            return_trajectory=False,
            precomputed=None,
            velocities_in=V0,
        ).dict()
        atoms2 = build_atoms(Z, d1["positions"], cell=None)
        d2 = _md_run(
            atoms2,
            steps=1,
            temperature=0.0,
            timestep_fs=1.0,
            friction=0.02,
            calculator=RelaxCalculatorName.lj,
            return_trajectory=False,
            precomputed=None,
            velocities_in=d1["velocities"],
        ).dict()

        # Second WS MD frame (skip init + first MD frame)
        ws_res = asyncio.run(_ws_nth_md_step("ws://127.0.0.1:8000/ws", Z, R, 2))

        np.testing.assert_allclose(
            np.array(d2["positions"]), ws_res["positions"], rtol=0, atol=1e-10
        )
        np.testing.assert_allclose(
            np.array(d2["velocities"]), ws_res["velocities"], rtol=0, atol=1e-10
        )
        np.testing.assert_allclose(
            np.array(d2["forces"]), ws_res["forces"], rtol=0, atol=1e-10
        )
    finally:
        try:
            serve.shutdown()
        except Exception:
            pass
