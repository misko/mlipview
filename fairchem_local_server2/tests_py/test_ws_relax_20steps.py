import asyncio
import json
import os
from typing import List, Tuple

import pytest
import websockets
from ase import Atoms  # type: ignore
from ase.io import read as ase_read  # type: ignore
from ase.optimize import BFGS  # type: ignore

from fairchem_local_server2 import session_pb2 as pb


def _load_water_xyz(path: str) -> tuple[List[int], List[List[float]]]:
    zs: List[int] = []
    pos: List[List[float]] = []
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f.readlines()]
    # Simple XYZ reader: skip first two header lines
    for ln in lines[2:]:
        if not ln:
            continue
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        x, y, z = map(float, parts[1:4])
        znum = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}.get(sym)
        if znum is None:
            raise ValueError(f"Unknown element symbol: {sym}")
        zs.append(int(znum))
        pos.append([x, y, z])
    return zs, pos


async def _recv_server_result(ws):
    # Receive a single binary ServerResult
    msg = await ws.recv()
    if isinstance(msg, (bytes, bytearray)):
        r = pb.ServerResult()
        r.ParseFromString(msg)
        return r
    # Text fallback: allow JSON logs to be interspersed (ignore)
    return None


async def run_relax_20(ws_url: str, xyz_path: str):
    zs, pos = _load_water_xyz(xyz_path)
    energies = []
    async with websockets.connect(ws_url, max_size=16 * 1024 * 1024) as ws:
        client_seq = 0
        # Init system via USER_INTERACTION (atomic_numbers + positions)
        init = pb.ClientAction()
        client_seq += 1
        init.seq = client_seq
        init.schema_version = 1
        ui = pb.ClientAction.UserInteraction()
        ui.atomic_numbers.extend(int(z) for z in zs)
        for p in pos:
            ui.positions.extend([float(p[0]), float(p[1]), float(p[2])])
        init.user_interaction.CopyFrom(ui)
        await ws.send(init.SerializeToString())

        # Wait for idle energy frame (positions omitted per protocol)
        for _ in range(10):
            r = await _recv_server_result(ws)
            if r is None:
                continue
            if r.WhichOneof("payload") == "frame":
                fr = r.frame
                # idle frame has no positions but includes energy
                has_energy = getattr(fr, "energy", None) is not None
                if has_energy and len(fr.positions) == 0:
                    energies.append(float(fr.energy))
                    # ACK idle energy frame to advance server window
                    client_seq += 1
                    ack = pb.ClientAction()
                    ack.seq = client_seq
                    ack.ack = int(getattr(r, "seq", 0))
                    ack.ping.CopyFrom(pb.ClientAction.Ping())
                    await ws.send(ack.SerializeToString())
                    break

        # Perform 20 relax steps: START_SIMULATION RELAX,
        # then read 1 frame and STOP
        for i in range(20):
            start = pb.ClientAction()
            client_seq += 1
            start.seq = client_seq
            start.schema_version = 1
            st = pb.ClientAction.Start()
            st.simulation_type = pb.ClientAction.Start.SimType.RELAX
            sp = pb.SimulationParams()
            sp.calculator = "uma"
            sp.fmax = 0.05
            sp.max_step = 0.2
            sp.optimizer = "bfgs"
            st.simulation_params.CopyFrom(sp)
            start.start.CopyFrom(st)
            await ws.send(start.SerializeToString())

            # Expect one sim frame with positions/forces/energy
            r = None
            for _ in range(50):
                r = await _recv_server_result(ws)
                if r is None:
                    continue
                if r.WhichOneof("payload") == "frame":
                    fr = r.frame
                    if len(fr.positions) >= 3 * len(zs):
                        # capture energy for this relax step
                        if getattr(fr, "energy", None) is not None:
                            energies.append(float(fr.energy))
                        client_seq += 1
                        ack = pb.ClientAction()
                        ack.seq = client_seq
                        ack.ack = int(getattr(r, "seq", 0))
                        ack.ping.CopyFrom(pb.ClientAction.Ping())
                        await ws.send(ack.SerializeToString())
                        break

            # request stop to ensure single-step behavior
            stop = pb.ClientAction()
            client_seq += 1
            stop.seq = client_seq
            stop.stop.CopyFrom(pb.ClientAction.Stop())
            await ws.send(stop.SerializeToString())

    return energies


def _uma_direct_reference(xyz_path: str) -> Tuple[float, float]:
    """Compute UMA reference energies directly (no WS/server) on GPU.

    Uses ASE BFGS optimizer with the same params as our WS relax steps:
    - fmax: 0.05
    - max_step: 0.2
    - steps: 20
    """
    try:
        from fairchem.core import FAIRChemCalculator, pretrained_mlip  # type: ignore
    except Exception as e:  # pragma: no cover
        pytest.skip(f"fairchem not available: {e}")

    try:
        atoms = ase_read(xyz_path)
        if not isinstance(atoms, Atoms):
            raise ValueError("ase_read did not return Atoms")
    except Exception:
        # Fallback to water loader
        zs, pos = _load_water_xyz(xyz_path)
        atoms = Atoms(numbers=zs, positions=pos)

    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
    calc = FAIRChemCalculator(predictor, task_name="omol")
    atoms.calc = calc

    e0 = float(atoms.get_potential_energy())
    opt = BFGS(atoms, maxstep=0.2)
    # Limit to 20 steps with our convergence criterion
    opt.run(fmax=0.05, steps=20)
    eN = float(atoms.get_potential_energy())
    return e0, eN


def _uma_reference_trace(
    xyz_path: str, steps: int = 20
) -> Tuple[List[float], List[List[float]], List[List[float]]]:
    """Compute full UMA energy trace and positions.

    Includes initial + per-step energies, and the initial/final positions.

    Returns (energies, positions_initial, positions_final)

    Matches server-side behavior:
    - Optimizer: ASE BFGS with maxstep=0.2
    - Convergence criterion: fmax=0.05 (but capped at 'steps')
    """
    try:
        from fairchem.core import FAIRChemCalculator, pretrained_mlip  # type: ignore
    except Exception as e:  # pragma: no cover
        pytest.skip(f"fairchem not available: {e}")

    try:
        atoms = ase_read(xyz_path)
        if not isinstance(atoms, Atoms):
            raise ValueError("ase_read did not return Atoms")
    except Exception:
        zs, pos = _load_water_xyz(xyz_path)
        atoms = Atoms(numbers=zs, positions=pos)

    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
    calc = FAIRChemCalculator(predictor, task_name="omol")
    atoms.calc = calc

    trace: List[float] = []
    # Initial energy and positions
    trace.append(float(atoms.get_potential_energy()))
    pos_initial = atoms.get_positions().tolist()
    opt = BFGS(atoms, maxstep=0.2, logfile=None)
    for _ in range(int(steps)):
        opt.step()
        trace.append(float(atoms.get_potential_energy()))
    pos_final = atoms.get_positions().tolist()
    return trace, pos_initial, pos_final


def test_ws_relax_20steps_parity(ws_base_url: str, tmp_path):
    """
    Direct WS test equivalent of the frontend water UMA 20-step relax parity.
    - Initializes water system via USER_INTERACTION
    - Runs 20 single relax steps via START/STOP
    - Asserts first/last energies are within 5e-3 of stored UMA reference
    Prints energies for debugging/logging.
    """
    # Use session-scoped server provided by fixture
    ws_url = ws_base_url
    xyz_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "public",
        "molecules",
        "water.xyz",
    )
    # Compute UMA energies directly (no WS) for 20-step BFGS, and persist
    uma_trace, pos_init, pos_final = _uma_reference_trace(xyz_path, steps=20)
    # Persist UMA trace for web UI matching and inspection
    out_dir = os.path.join(os.getcwd(), "test-results")
    os.makedirs(out_dir, exist_ok=True)
    uma_out = os.path.join(out_dir, "water_uma_bfgs_20_trace.json")
    with open(uma_out, "w", encoding="utf-8") as f:
        json.dump(
            {
                "energies": uma_trace,
                "positions_initial": pos_init,
                "positions_final": pos_final,
            },
            f,
        )

    # Run via WS and collect energies (init + 20)
    ws_trace = asyncio.run(run_relax_20(ws_url, xyz_path))
    print("PY_WS_RELAX20_ENERGIES", json.dumps(ws_trace[:21]))
    # Persist WS trace as well
    ws_out = os.path.join(out_dir, "water_ws_bfgs_20_trace.json")
    with open(ws_out, "w", encoding="utf-8") as f:
        json.dump({"energies": ws_trace}, f)

    assert len(ws_trace) >= 21, "expected at least 21 energies (init + 20)"
    # Compare first and last energy to UMA direct reference
    assert (
        abs(ws_trace[0] - uma_trace[0]) < 5e-3
    ), f"first energy mismatch: {ws_trace[0]} vs {uma_trace[0]}"
    assert (
        abs(ws_trace[20] - uma_trace[20]) < 5e-3
    ), f"last energy mismatch: {ws_trace[20]} vs {uma_trace[20]}"
