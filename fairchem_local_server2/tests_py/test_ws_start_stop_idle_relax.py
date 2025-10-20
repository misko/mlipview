from __future__ import annotations

import asyncio
from typing import List

import pytest
import websockets

from fairchem_local_server2 import session_pb2 as pb


def _make_vec3(arr):
    v = pb.Vec3()
    v.v.extend([float(arr[0]), float(arr[1]), float(arr[2])])
    return v


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
            init.type = pb.ClientAction.Type.USER_INTERACTION
            init.atomic_numbers.extend(Z)
            for p in R:
                init.positions.append(_make_vec3(p))
            await ws.send(init.SerializeToString())

            # Drain until we see an energy (idle compute) and ack
            for _ in range(20):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                has_energy = hasattr(r, "energy") and (
                    getattr(r, "energy", None) is not None
                )
                if has_energy:
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.type = pb.ClientAction.Type.PING
                    try:
                        ack.ack = int(getattr(r, "seq", 0))
                    except Exception:
                        pass
                    await ws.send(ack.SerializeToString())
                    break

            # Start MD continuous; collect >=20 frames
            start = pb.ClientAction()
            seq += 1
            start.seq = seq
            start.type = pb.ClientAction.Type.START_SIMULATION
            start.simulation_type = pb.ClientAction.SimType.MD
            sp = pb.SimulationParams()
            sp.calculator = "uma"
            sp.temperature = 1500.0
            sp.timestep_fs = 1.0
            start.simulation_params.CopyFrom(sp)
            await ws.send(start.SerializeToString())

            md_frames = 0
            while md_frames < 20:
                r = await _recv_any(ws, timeout=20.0)
                if r is None:
                    continue
                # Expect positions present during simulation
                if len(r.positions) >= len(Z):
                    md_frames += 1
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.type = pb.ClientAction.Type.PING
                    ack.ack = int(getattr(r, "seq", 0))
                    await ws.send(ack.SerializeToString())

            # Stop MD; expect a stop message/flag frame
            stop = pb.ClientAction()
            seq += 1
            stop.seq = seq
            stop.type = pb.ClientAction.Type.STOP_SIMULATION
            await ws.send(stop.SerializeToString())

            saw_stop = False
            for _ in range(10):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                msg = getattr(r, "message", "") or ""
                has_flag = hasattr(r, "simulation_stopped") and (
                    getattr(r, "simulation_stopped", False) is True
                )
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
                req.type = pb.ClientAction.Type.USER_INTERACTION
                # perturb a single H slightly to simulate drag
                dx = 0.01 * ((i % 2) * 2 - 1)
                pos = [list(p) for p in R]
                pos[1][0] += dx
                for p in pos:
                    req.positions.append(_make_vec3(p))
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
                    e_ok = hasattr(r, "energy") and (
                        getattr(r, "energy", None) is not None
                    )
                    no_pos = len(r.positions) == 0
                    if e_ok and no_pos:
                        idle_results += 1
                        seq += 1
                        ack = pb.ClientAction()
                        ack.seq = seq
                        ack.type = pb.ClientAction.Type.PING
                        ack.ack = int(getattr(r, "seq", 0))
                        await ws.send(ack.SerializeToString())
                        break

            assert idle_results >= 12, f"Expected >=12 idle results, got {idle_results}"

            # Start RELAX; collect >=20 frames
            start2 = pb.ClientAction()
            seq += 1
            start2.seq = seq
            start2.type = pb.ClientAction.Type.START_SIMULATION
            start2.simulation_type = pb.ClientAction.SimType.RELAX
            sp2 = pb.SimulationParams()
            sp2.calculator = "uma"
            sp2.fmax = 0.05
            sp2.max_step = 0.2
            sp2.optimizer = "bfgs"
            start2.simulation_params.CopyFrom(sp2)
            await ws.send(start2.SerializeToString())

            rx_frames = 0
            while rx_frames < 20:
                r = await _recv_any(ws, timeout=20.0)
                if r is None:
                    continue
                if len(r.positions) >= len(Z):
                    rx_frames += 1
                    seq += 1
                    ack = pb.ClientAction()
                    ack.seq = seq
                    ack.type = pb.ClientAction.Type.PING
                    ack.ack = int(getattr(r, "seq", 0))
                    await ws.send(ack.SerializeToString())

            # Stop RELAX; expect stop indicator
            stop2 = pb.ClientAction()
            seq += 1
            stop2.seq = seq
            stop2.type = pb.ClientAction.Type.STOP_SIMULATION
            await ws.send(stop2.SerializeToString())

            saw_stop2 = False
            for _ in range(10):
                r = await _recv_any(ws, timeout=10.0)
                if r is None:
                    continue
                msg = getattr(r, "message", "") or ""
                has_flag = hasattr(r, "simulation_stopped") and (
                    getattr(r, "simulation_stopped", False) is True
                )
                if msg == "SIMULATION_STOPPED" or has_flag:
                    saw_stop2 = True
                    break
            assert saw_stop2, "Expected stop indicator after STOP_SIMULATION (RELAX)"

    asyncio.run(run(ws_base_url))
