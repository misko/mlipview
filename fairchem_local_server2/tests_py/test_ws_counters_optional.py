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
        init.schema_version = 1
        ui = pb.ClientAction.UserInteraction()
        ui.atomic_numbers.extend([1, 1])
        for p in ([0.0, 0.0, 0.0], [0.9, 0.0, 0.0]):
            ui.positions.extend([float(p[0]), float(p[1]), float(p[2])])
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # start relax (one-by-one frames)
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.schema_version = 1
        st = pb.ClientAction.Start()
        st.simulation_type = pb.ClientAction.Start.SimType.RELAX
        # Attach optional counters if fields exist
        if hasattr(start, "user_interaction_count"):
            setattr(start, "user_interaction_count", 5)
        if hasattr(start, "sim_step"):
            setattr(start, "sim_step", 42)
        start.start.CopyFrom(st)
        await ws.send(start.SerializeToString())

        # receive first produced simulation frame
        while True:
            data = await asyncio.wait_for(ws.recv(), timeout=15.0)
            if not isinstance(data, (bytes, bytearray)):
                continue
            res = pb.ServerResult()
            res.ParseFromString(data)
            if res.WhichOneof("payload") != "frame":
                continue
            return res


@pytest.mark.timeout(30)
def test_ws_optional_counters_present_or_skipped(ws_base_url: str):
    res = asyncio.run(_recv_one(ws_base_url))
    # If sim_step exists in compiled schema, expect >= 0 for simulation frames
    # (server may start from 0)
    if hasattr(res, "sim_step"):
        assert int(getattr(res, "sim_step", 0)) >= 0
    # If user_interaction_count exists, expect echo >= 0
    if hasattr(res, "user_interaction_count"):
        assert int(getattr(res, "user_interaction_count", 0)) >= 0
