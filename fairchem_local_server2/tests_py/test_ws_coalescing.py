from __future__ import annotations

import asyncio
import time

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _run_coalescing(uri: str, updates: int = 30) -> tuple[int, int]:
    async with websockets.connect(uri, max_size=16 * 1024 * 1024) as ws:
        seq = 0

        # Minimal init with two atoms
        Z = [1, 1]
        R = [[0.0, 0.0, 0.0], [0.9, 0.0, 0.0]]
        seq += 1
        init = pb.ClientAction()
        init.seq = seq
        init.schema_version = 1
        ui = pb.UserInteractionSparse()
        # natoms and full arrays sent as sparse deltas
        ui.natoms = len(Z)
        z_delta = pb.IntDelta()
        z_delta.indices.extend([0, 1])
        z_delta.values.extend([int(Z[0]), int(Z[1])])
        ui.atomic_numbers.CopyFrom(z_delta)
        r_delta = pb.Vec3Delta()
        r_delta.indices.extend([0, 1])
        r_delta.coords.extend(
            [
                float(R[0][0]),
                float(R[0][1]),
                float(R[0][2]),
                float(R[1][0]),
                float(R[1][1]),
                float(R[1][2]),
            ]
        )
        ui.positions.CopyFrom(r_delta)
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # Drain one idle result (frame or notice), ack
        for _ in range(20):
            try:
                data = await asyncio.wait_for(ws.recv(), timeout=2.0)
            except asyncio.TimeoutError:
                continue
            if not isinstance(data, (bytes, bytearray)):
                continue
            res = pb.ServerResult()
            res.ParseFromString(data)
            which = res.WhichOneof("payload")
            if which == "frame":
                fr = res.frame
                if getattr(fr, "energy", None) is not None and len(fr.positions) == 0:
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.ack = int(res.seq)
                    # No ping payload in current schema; ack-only is valid
                    await ws.send(ack.SerializeToString())
                    break
            elif which == "notice":
                seq += 1
                ack = pb.ClientAction()
                ack.seq = seq
                ack.ack = int(res.seq)
                # No ping payload in current schema; ack-only is valid
                await ws.send(ack.SerializeToString())
                break

        # Send many rapid USER_INTERACTION updates without waiting
        # to trigger server-side coalescing
        for i in range(int(updates)):
            seq += 1
            msg = pb.ClientAction()
            msg.seq = seq
            msg.schema_version = 1
            ui2 = pb.UserInteractionSparse()
            # Alternate small perturbation on H position (index 1)
            dx = 0.01 * ((i % 2) * 2 - 1)
            r_delta2 = pb.Vec3Delta()
            r_delta2.indices.extend([1])
            r_delta2.coords.extend(
                [
                    float(R[1][0] + dx),
                    float(R[1][1]),
                    float(R[1][2]),
                ]
            )
            ui2.positions.CopyFrom(r_delta2)
            msg.user_interaction.CopyFrom(ui2)
            try:
                msg.user_interaction_count = i + 1
            except Exception:
                pass
            await ws.send(msg.SerializeToString())

        # Drain frames for a short window,
        # count idle-only frames and track last UIC echo
        idle_frames = 0
        last_uic = 0
        t0 = time.time()
        while time.time() - t0 < 6.0:
            try:
                data = await asyncio.wait_for(ws.recv(), timeout=0.5)
            except asyncio.TimeoutError:
                # No more frames readily available
                break
            if not isinstance(data, (bytes, bytearray)):
                continue
            res = pb.ServerResult()
            res.ParseFromString(data)
            which = res.WhichOneof("payload")
            if which == "frame":
                fr = res.frame
                if len(fr.positions) == 0:
                    idle_frames += 1
                try:
                    last_uic = int(getattr(res, "user_interaction_count", 0))
                except Exception:
                    last_uic = 0
                # ACK to advance server window
                seq += 1
                ack = pb.ClientAction()
                ack.seq = seq
                ack.ack = int(res.seq)
                # No ping payload in current schema; ack-only is valid
                await ws.send(ack.SerializeToString())
            elif which == "notice":
                try:
                    last_uic = int(getattr(res, "user_interaction_count", 0))
                except Exception:
                    last_uic = 0
                # ACK to advance server window
                seq += 1
                ack = pb.ClientAction()
                ack.seq = seq
                ack.ack = int(res.seq)
                # No ping payload in current schema; ack-only is valid
                await ws.send(ack.SerializeToString())

        return last_uic, idle_frames


@pytest.mark.timeout(60)
def test_ws_coalescing_idle_user_interactions(ws_base_url: str):
    last_uic, idle_frames = asyncio.run(_run_coalescing(ws_base_url, updates=30))
    # Expect to have processed the latest update
    assert last_uic == 30
    # Expect significantly fewer idle frames than updates due to coalescing
    assert idle_frames <= max(5, 30 // 4)
