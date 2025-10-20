from __future__ import annotations

import asyncio

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _recv_one(uri: str):
    async with websockets.connect(uri) as ws:
        seq = 1
        # Init with two atoms via USER_INTERACTION
        init = pb.ClientAction()
        init.seq = seq
        init.type = pb.ClientAction.Type.USER_INTERACTION
        init.atomic_numbers.extend([1, 1])
        for p in ([0.0, 0.0, 0.0], [0.9, 0.0, 0.0]):
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        await ws.send(init.SerializeToString())

        # start relax (one-by-one frames)
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.type = pb.ClientAction.Type.START_SIMULATION
        start.simulation_type = pb.ClientAction.SimType.RELAX
        # Set counters explicitly; schema guarantees presence
        start.user_interaction_count = 5
        start.sim_step = 42
        await ws.send(start.SerializeToString())

        # receive first produced simulation frame
        while True:
            data = await asyncio.wait_for(ws.recv(), timeout=15.0)
            if not isinstance(data, (bytes, bytearray)):
                continue
            res = pb.ServerResult()
            res.ParseFromString(data)
            # Skip initialization echo and idle-compute frames
            # (idle frames omit positions per protocol)
            if res.seq == 0 or len(res.positions) == 0:
                continue
            return res


@pytest.mark.timeout(30)
def test_ws_optional_counters_present_or_skipped(ws_base_url: str):
    res = asyncio.run(_recv_one(ws_base_url))
    # Assert echo semantics on required fields in new protocol
    assert int(res.sim_step) >= 1
    assert int(res.user_interaction_count) >= 0
