from __future__ import annotations

import asyncio
from typing import List

import numpy as np
import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _send_ui(
    ws,
    seq: int,
    *,
    natoms: int | None,
    z_deltas: List[tuple[int, int]] | None,
    r_deltas: List[tuple[int, List[float]]] | None,
    v_zero_for: List[int] | None,
) -> int:
    """Helper to send a sparse USER_INTERACTION message.
    Returns next seq.
    """
    msg = pb.ClientAction()
    msg.seq = seq
    msg.schema_version = 1
    ui = pb.UserInteractionSparse()
    if natoms and natoms > 0:
        ui.natoms = int(natoms)
    if z_deltas:
        inz = pb.IntDelta()
        for i, z in z_deltas:
            inz.indices.append(int(i))
            inz.values.append(int(z))
        ui.atomic_numbers.CopyFrom(inz)
    if r_deltas:
        vr = pb.Vec3Delta()
        for i, p in r_deltas:
            vr.indices.append(int(i))
            vr.coords.extend([float(p[0]), float(p[1]), float(p[2])])
        ui.positions.CopyFrom(vr)
    if v_zero_for:
        vv = pb.Vec3Delta()
        for i in v_zero_for:
            vv.indices.append(int(i))
            vv.coords.extend([0.0, 0.0, 0.0])
        ui.velocities.CopyFrom(vv)
    msg.user_interaction.CopyFrom(ui)
    await ws.send(msg.SerializeToString())
    return seq + 1


@pytest.mark.timeout(60)
def test_ws_sparse_preserve_deltas_across_resizes(ws_base_url: str):
    """
    Regression test: applying sparse deltas incrementally (index 0 first, then higher
    indices that trigger internal expand-only sizing) should preserve earlier values and
    perform a valid idle compute (no INVALID_ATOM_Z) once the geometry is complete.
    """

    async def _run():
        async with websockets.connect(ws_base_url) as ws:
            seq = 1
            # Step 1: partial init for index 0 only
            seq = await _send_ui(
                ws,
                seq,
                natoms=3,
                z_deltas=[(0, 1)],
                r_deltas=[(0, [0.0, 0.0, 0.0])],
                v_zero_for=[0],
            )

            # Step 2: complete remaining indices (1: H, 2: O) and positions
            seq = await _send_ui(
                ws,
                seq,
                natoms=None,  # no need to resend natoms
                z_deltas=[(1, 1), (2, 8)],
                r_deltas=[(1, [0.96, 0.0, 0.0]), (2, [-0.24, 0.93, 0.0])],
                v_zero_for=[1, 2],
            )

            # Expect an idle compute frame with forces present and positions sized to 3
            import time

            deadline = time.monotonic() + 15.0
            last_frame: pb.ServerResult | None = None
            while time.monotonic() < deadline:
                data = await asyncio.wait_for(ws.recv(), timeout=5.0)
                res = pb.ServerResult()
                res.ParseFromString(data)
                if res.WhichOneof("payload") != "frame":
                    continue
                fr = res.frame
                P = np.fromiter(fr.positions, dtype=np.float64)
                F = np.fromiter(fr.forces, dtype=np.float64)
                if P.size == 0 or F.size == 0:
                    continue
                # shape checks
                assert P.size % 3 == 0 and F.size % 3 == 0
                assert P.reshape(-1, 3).shape[0] == 3
                assert F.reshape(-1, 3).shape[0] == 3
                # finite energy if provided
                if fr.HasField("energy"):
                    assert np.isfinite(fr.energy)
                last_frame = res
                break

            assert last_frame is not None, "did not receive a valid frame in time"

            # ACK the frame
            ack = pb.ClientAction()
            ack.seq = seq
            ack.ack = int(last_frame.seq)
            await ws.send(ack.SerializeToString())

            # STOP (clean shutdown path)
            stop = pb.ClientAction()
            stop.seq = seq + 1
            stop.stop.CopyFrom(pb.ClientAction.Stop())
            await ws.send(stop.SerializeToString())

    asyncio.run(_run())
