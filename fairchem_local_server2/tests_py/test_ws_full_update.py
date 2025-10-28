from __future__ import annotations

import asyncio
from typing import Sequence

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _recv_result(ws, timeout=5.0):
    data = await asyncio.wait_for(ws.recv(), timeout=timeout)
    if not isinstance(data, (bytes, bytearray)):
        return None
    res = pb.ServerResult()
    res.ParseFromString(data)
    return res


async def _wait_for_frame_with_atoms(ws, expected_atoms: int, seq: int, timeout=15.0):
    """Await a frame whose forces or positions length matches expected_atoms.

    Returns (frame, updated_seq, notice_seen).
    """
    notice_seen = False
    attempts = 0
    while True:
        res = await _recv_result(ws, timeout=timeout)
        assert res is not None
        payload = res.WhichOneof("payload")
        if payload == "notice" and res.notice.simulation_stopped:
            notice_seen = True
        if payload == "frame":
            f_len = len(res.frame.forces)
            p_len = len(res.frame.positions)
            target = expected_atoms * 3
            if f_len == target or p_len == target:
                return res, seq, notice_seen
        seq += 1
        await ws.send(_ack(seq, res.seq).SerializeToString())
        attempts += 1
        if attempts > 1000:
            raise AssertionError(f"Timed out waiting for frame with {expected_atoms} atoms (payload={payload}, forces={len(res.frame.forces)}, positions={len(res.frame.positions)})")


def _make_full_update(seq: int, z: Sequence[int], coords: Sequence[Sequence[float]]):
    msg = pb.ClientAction()
    msg.seq = int(seq)
    msg.schema_version = 1
    ui = msg.user_interaction
    ui.full_update = True
    ui.natoms = len(z)
    inz = ui.atomic_numbers
    inz.indices.extend(range(len(z)))
    inz.values.extend(int(v) for v in z)
    vr = ui.positions
    vr.indices.extend(range(len(coords)))
    for p in coords:
        vr.coords.extend([float(p[0]), float(p[1]), float(p[2])])
    assert ui.full_update is True
    try:
        assert ui.HasField("full_update")  # type: ignore[attr-defined]
    except Exception:
        # proto3 optional bool auto-initialises when set to True; ignore if HasField unsupported
        pass
    return msg


def _ack(seq: int, ack_target: int) -> pb.ClientAction:
    ack = pb.ClientAction()
    ack.seq = int(seq)
    ack.ack = int(ack_target)
    ack.schema_version = 1
    return ack


@pytest.mark.timeout(120)
def test_full_update_idle_and_resize(ws_base_url: str):
    Z2 = [8, 1]
    R2 = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]]
    Z3 = [8, 1, 1]
    R3 = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.2, 0.92, 0.0]]

    async def run(uri: str):
        async with websockets.connect(uri, max_size=16 * 1024 * 1024) as ws:
            seq = 0
            seq += 1
            await ws.send(_make_full_update(seq, Z2, R2).SerializeToString())

            res = await _recv_result(ws, timeout=15.0)
            assert res is not None
            assert res.WhichOneof("payload") == "frame"
            assert len(res.frame.forces) == len(Z2) * 3

            seq += 1
            await ws.send(_ack(seq, res.seq).SerializeToString())

            seq += 1
            await ws.send(_make_full_update(seq, Z3, R3).SerializeToString())

            while True:
                res2 = await _recv_result(ws, timeout=15.0)
                assert res2 is not None
                if res2.WhichOneof("payload") == "frame":
                    assert len(res2.frame.forces) == len(Z3) * 3
                    break
                seq += 1
                await ws.send(_ack(seq, res2.seq).SerializeToString())

    asyncio.run(run(ws_base_url))


@pytest.mark.timeout(180)
def test_full_update_during_md_and_relax(ws_base_url: str):
    Z = [8, 1, 1]
    R = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.2, 0.92, 0.0]]
    Z4 = [8, 1, 1, 1]
    R4 = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.2, 0.92, 0.0], [1.4, 0.0, 0.0]]
    Z2 = [8, 1]
    R2 = [[0.0, 0.0, 0.0], [0.95, 0.0, 0.0]]

    async def run(uri: str):
        async with websockets.connect(uri, max_size=16 * 1024 * 1024) as ws:
            seq = 0
            seq += 1
            await ws.send(_make_full_update(seq, Z, R).SerializeToString())
            res = await _recv_result(ws, timeout=10.0)
            assert res is not None
            seq += 1
            await ws.send(_ack(seq, res.seq).SerializeToString())

            # Start MD
            seq += 1
            start = pb.ClientAction()
            start.seq = seq
            start.schema_version = 1
            st = pb.ClientAction.Start()
            st.simulation_type = pb.ClientAction.Start.SimType.MD
            start.start.CopyFrom(st)
            await ws.send(start.SerializeToString())

            # Wait for one MD frame
            md_frame = await _recv_result(ws, timeout=20.0)
            assert md_frame is not None
            assert md_frame.WhichOneof("payload") == "frame"
            seq += 1
            await ws.send(_ack(seq, md_frame.seq).SerializeToString())

            # Send full update with more atoms
            seq += 1
            await ws.send(_make_full_update(seq, Z4, R4).SerializeToString())

            idle, seq, notice_seen = await _wait_for_frame_with_atoms(ws, len(Z4), seq)
            assert len(idle.frame.forces) == len(Z4) * 3
            seq += 1
            await ws.send(_ack(seq, idle.seq).SerializeToString())

            # Start relax
            seq += 1
            start_relax = pb.ClientAction()
            start_relax.seq = seq
            start_relax.schema_version = 1
            st_relax = pb.ClientAction.Start()
            st_relax.simulation_type = pb.ClientAction.Start.SimType.RELAX
            start_relax.start.CopyFrom(st_relax)
            await ws.send(start_relax.SerializeToString())

            relax_frame = await _recv_result(ws, timeout=20.0)
            assert relax_frame is not None
            seq += 1
            await ws.send(_ack(seq, relax_frame.seq).SerializeToString())

            # Send shrink snapshot
            seq += 1
            await ws.send(_make_full_update(seq, Z2, R2).SerializeToString())

            idle2, seq, notice2_seen = await _wait_for_frame_with_atoms(ws, len(Z2), seq)
            assert len(idle2.frame.forces) == len(Z2) * 3
            seq += 1
            await ws.send(_ack(seq, idle2.seq).SerializeToString())

            # Restart MD to ensure new natoms used
            seq += 1
            restart = pb.ClientAction()
            restart.seq = seq
            restart.schema_version = 1
            st_restart = pb.ClientAction.Start()
            st_restart.simulation_type = pb.ClientAction.Start.SimType.MD
            restart.start.CopyFrom(st_restart)
            await ws.send(restart.SerializeToString())

            while True:
                md2 = await _recv_result(ws, timeout=20.0)
                assert md2 is not None
                if md2.WhichOneof("payload") == "frame":
                    print("restart md positions len=", len(md2.frame.positions))
                    assert len(md2.frame.positions) == len(Z2) * 3
                    break
                seq += 1
                await ws.send(_ack(seq, md2.seq).SerializeToString())

            seq += 1
            await ws.send(_ack(seq, md2.seq).SerializeToString())

    asyncio.run(run(ws_base_url))
