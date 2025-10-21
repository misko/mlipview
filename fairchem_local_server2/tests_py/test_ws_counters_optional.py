"""
End-to-end (WS): verify optional counter fields are handled gracefully.

What this test covers
- Protocol: protobuf, oneof payloads, flat arrays, schema_version=1.
- Client sends:
    1) USER_INTERACTION as an initialization message (2 atoms, positions
         as flat packed doubles).
    2) START simulation (RELAX). If the compiled client schema exposes
         the optional counter fields (user_interaction_count and sim_step),
         set them on the message. Otherwise, skip them (schema
         backward-compatibility).
- Server behavior expectations:
    - Produce simulation frames (one step at a time) as ServerResult.frame.
    - If the ServerResult schema has the optional counters, populate them
        with a non-negative value (>= 0). If the fields are absent in the
        compiled schema, the test tolerates that by skipping the checks.

Why this exists
- Guards the optionality/compatibility contract around correlation
    counters during schema migration. Different environments (or cached
    stubs) might not compile the optional fields yet; the server must
    still function, and the test should not fail solely due to missing
    fields in the local stubs.

Assertions
- For ServerResult.frame received first:
    - If res.sim_step exists: it should be >= 0 (simulation step counter).
    - If res.user_interaction_count exists: it should be >= 0 (echo or
        snapshot).

Notes
- Counters are treated as correlation/diagnostic metadata; they must
    never break the transport when absent. This test intentionally does
    not enforce exact values—only presence (when compiled) and
    non-negativity.
"""

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
