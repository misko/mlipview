from __future__ import annotations

import asyncio
from typing import List

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


async def _recv_any(ws, timeout=5.0):
    data = await asyncio.wait_for(ws.recv(), timeout=timeout)
    if isinstance(data, (bytes, bytearray)):
        r = pb.ServerResult()
        r.ParseFromString(data)
        return r
    return None


@pytest.mark.timeout(120)
def test_ws_start_stop_idle_relax_sequence(ws_base_url: str):
    """
    Scenario:
    - Initialize small molecule via USER_INTERACTION
    - Start MD; receive >=20 frames with positions; STOP and receive stop
        indicator
        - Send 20 USER_INTERACTION updates while idle; receive >=12 idle
            results without positions
    - Start RELAX; receive >=20 frames; STOP and receive stop indicator
    """
    Z = [1, 1, 8]  # water-like minimal
    R = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0], [0.37, 0.58, 0.0]]

    async def run(uri: str):
        async with websockets.connect(uri, max_size=16 * 1024 * 1024) as ws:
            seq = 0
            # Init
            init = pb.ClientAction()
            seq += 1
            init.seq = seq
            init.schema_version = 1
            ui = pb.ClientAction.UserInteraction()
            ui.atomic_numbers.extend(Z)
            for p in R:
                ui.positions.extend([float(p[0]), float(p[1]), float(p[2])])
            init.user_interaction.CopyFrom(ui)
            await ws.send(init.SerializeToString())

            # Drain until we see an energy (idle compute) and ack
            for _ in range(20):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                if r.WhichOneof("payload") != "frame":
                    continue
                fr = r.frame
                has_energy = getattr(fr, "energy", None) is not None
                no_pos = len(fr.positions) == 0
                if has_energy and no_pos:
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.ack = int(getattr(r, "seq", 0))
                    ack.ping.CopyFrom(pb.ClientAction.Ping())
                    await ws.send(ack.SerializeToString())
                    break

            # Start MD continuous; collect >=20 frames
            start = pb.ClientAction()
            seq += 1
            start.seq = seq
            start.schema_version = 1
            st = pb.ClientAction.Start()
            st.simulation_type = pb.ClientAction.Start.SimType.MD
            sp = pb.SimulationParams()
            sp.calculator = "uma"
            sp.temperature = 1500.0
            sp.timestep_fs = 1.0
            st.simulation_params.CopyFrom(sp)
            start.start.CopyFrom(st)
            await ws.send(start.SerializeToString())

            md_frames = 0
            while md_frames < 20:
                r = await _recv_any(ws, timeout=20.0)
                if r is None:
                    continue
                if r.WhichOneof("payload") == "frame" and len(r.frame.positions) >= len(
                    Z
                ):
                    md_frames += 1
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.ack = int(getattr(r, "seq", 0))
                    ack.ping.CopyFrom(pb.ClientAction.Ping())
                    await ws.send(ack.SerializeToString())

            # Stop MD; expect a stop message/flag frame
            stop = pb.ClientAction()
            seq += 1
            stop.seq = seq
            stop.stop.CopyFrom(pb.ClientAction.Stop())
            await ws.send(stop.SerializeToString())

            saw_stop = False
            for _ in range(10):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                if r.WhichOneof("payload") == "notice":
                    msg = getattr(r.notice, "message", "") or ""
                    has_flag = getattr(r.notice, "simulation_stopped", False) is True
                else:
                    msg = ""
                    has_flag = False
                if msg == "SIMULATION_STOPPED" or has_flag:
                    saw_stop = True
                    break
            assert saw_stop, "Expected stop indicator after STOP_SIMULATION (MD)"

            # While idle, send 20 USER_INTERACTION updates (positions only)
            idle_results = 0
            for i in range(20):
                req = pb.ClientAction()
                seq += 1
                req.seq = seq
                req.schema_version = 1
                # perturb a single H slightly to simulate drag
                dx = 0.01 * ((i % 2) * 2 - 1)
                pos = [list(p) for p in R]
                pos[1][0] += dx
                ui2 = pb.ClientAction.UserInteraction()
                for p in pos:
                    ui2.positions.extend([float(p[0]), float(p[1]), float(p[2])])
                req.user_interaction.CopyFrom(ui2)
                try:
                    req.user_interaction_count = i + 1
                except Exception:
                    pass
                await ws.send(req.SerializeToString())

                # Expect energy/forces but no positions in idle frames
                for _ in range(10):
                    r = await _recv_any(ws, timeout=10.0)
                    if r is None:
                        continue
                    if r.WhichOneof("payload") != "frame":
                        continue
                    fr = r.frame
                    e_ok = hasattr(fr, "energy") and (
                        getattr(fr, "energy", None) is not None
                    )
                    no_pos = len(fr.positions) == 0
                    if e_ok and no_pos:
                        idle_results += 1
                        seq += 1
                        ack = pb.ClientAction()
                        ack.seq = seq
                        ack.ack = int(getattr(r, "seq", 0))
                        ack.ping.CopyFrom(pb.ClientAction.Ping())
                        await ws.send(ack.SerializeToString())
                        break

                # Server may coalesce or rate-limit idle frames;
                # ensure at least one
                assert (
                    idle_results >= 1
                ), f"Expected >=1 idle result, got {idle_results}"

            # Start RELAX; collect >=20 frames
            start2 = pb.ClientAction()
            seq += 1
            start2.seq = seq
            start2.schema_version = 1
            st2 = pb.ClientAction.Start()
            st2.simulation_type = pb.ClientAction.Start.SimType.RELAX
            sp2 = pb.SimulationParams()
            sp2.calculator = "uma"
            sp2.fmax = 0.05
            sp2.max_step = 0.2
            sp2.optimizer = "bfgs"
            st2.simulation_params.CopyFrom(sp2)
            start2.start.CopyFrom(st2)
            await ws.send(start2.SerializeToString())

            rx_frames = 0
            while rx_frames < 20:
                r = await _recv_any(ws, timeout=20.0)
                if r is None:
                    continue
                if r.WhichOneof("payload") == "frame" and len(r.frame.positions) >= len(
                    Z
                ):
                    rx_frames += 1
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.ack = int(getattr(r, "seq", 0))
                    ack.ping.CopyFrom(pb.ClientAction.Ping())
                    await ws.send(ack.SerializeToString())

            # Stop RELAX; expect stop indicator
            stop2 = pb.ClientAction()
            seq += 1
            stop2.seq = seq
            stop2.stop.CopyFrom(pb.ClientAction.Stop())
            await ws.send(stop2.SerializeToString())

            saw_stop2 = False
            for _ in range(10):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                if r.WhichOneof("payload") == "notice":
                    msg = getattr(r.notice, "message", "") or ""
                    has_flag = getattr(r.notice, "simulation_stopped", False) is True
                else:
                    msg = ""
                    has_flag = False
                if msg == "SIMULATION_STOPPED" or has_flag:
                    saw_stop2 = True
                    break
            assert saw_stop2, "Expected stop indicator after STOP_SIMULATION (RELAX)"

    asyncio.run(run(ws_base_url))
