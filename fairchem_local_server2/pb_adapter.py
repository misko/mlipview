from __future__ import annotations

from typing import Optional

from fairchem_local_server2.session_models import (
    ClientAction as JsonAction,
    ServerResult as JsonResult,
    SimulationParams,
)


def _has_pb():
    try:
        from fairchem_local_server2 import session_pb2  # type: ignore

        return session_pb2
    except Exception:
        return None


def decode_action_from_bytes(b: bytes) -> Optional[JsonAction]:
    """Attempt to decode protobuf ClientAction bytes into JSON model.

    Returns None if protobuf module is not available or decoding fails.
    """
    mod = _has_pb()
    if not mod:
        return None
    try:
        msg = mod.ClientAction()
        msg.ParseFromString(b)
        # Map to Pydantic model
        params = None
        if msg.HasField("simulation_params"):
            sp = msg.simulation_params
            params = SimulationParams(
                calculator=sp.calculator or "uma",
                temperature=sp.temperature,
                timestep_fs=sp.timestep_fs,
                friction=sp.friction,
                fmax=sp.fmax,
                max_step=sp.max_step,
                optimizer=(sp.optimizer or "bfgs"),
            )
        pos = [[v.v[0], v.v[1], v.v[2]] for v in msg.positions]
        vel = (
            [[v.v[0], v.v[1], v.v[2]] for v in msg.velocities]
            if len(msg.velocities)
            else None
        )
        cell = None
        if msg.HasField("cell") and len(msg.cell.m) == 9:
            m = msg.cell.m
            cell = [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]]

        action_type_map = {
            mod.ClientAction.INIT_SYSTEM: "init_system",
            mod.ClientAction.UPDATE_POSITIONS: "update_positions",
            mod.ClientAction.START_SIMULATION: "start_simulation",
            mod.ClientAction.STOP_SIMULATION: "stop_simulation",
            mod.ClientAction.PING: "ping",
        }
        sim_type_map = {
            mod.ClientAction.MD: "md",
            mod.ClientAction.RELAX: "relax",
        }

        ja = JsonAction(
            seq=int(msg.seq or 0),
            ack=int(msg.ack) if msg.HasField("ack") else None,
            type=action_type_map.get(msg.type, "ping"),
            atomic_numbers=list(msg.atomic_numbers) or None,
            positions=pos or None,
            velocities=vel,
            cell=cell,
            simulation_type=sim_type_map.get(msg.simulation_type)
            if msg.HasField("simulation_type")
            else None,
            simulation_params=params,
        )
        return ja
    except Exception:
        return None


def encode_result_to_bytes(res: JsonResult) -> Optional[bytes]:
    mod = _has_pb()
    if not mod:
        return None
    try:
        msg = mod.ServerResult()
        msg.seq = int(res.seq)
        if res.client_seq is not None:
            msg.client_seq = int(res.client_seq)
        for p in res.positions:
            v = mod.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            msg.positions.append(v)
        for f in res.forces:
            v = mod.Vec3()
            v.v.extend([float(f[0]), float(f[1]), float(f[2])])
            msg.forces.append(v)
        for vv in res.velocities:
            v = mod.Vec3()
            v.v.extend([float(vv[0]), float(vv[1]), float(vv[2])])
            msg.velocities.append(v)
        if res.cell is not None:
            m = mod.Mat3()
            flat = [
                res.cell[0][0],
                res.cell[0][1],
                res.cell[0][2],
                res.cell[1][0],
                res.cell[1][1],
                res.cell[1][2],
                res.cell[2][0],
                res.cell[2][1],
                res.cell[2][2],
            ]
            m.m.extend([float(x) for x in flat])
            msg.cell.CopyFrom(m)
        if res.message:
            msg.message = res.message
        return msg.SerializeToString()
    except Exception:
        return None
