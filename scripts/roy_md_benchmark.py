from __future__ import annotations

import argparse
import asyncio
import time
from pathlib import Path
from typing import Tuple, List, Optional
import time as _time

import json
import websockets
from fairchem_local_server2 import session_pb2 as pb

try:
    # Prefer ASE for robust XYZ parsing if available
    from ase.io import read as ase_read  # type: ignore
    from ase.atoms import Atoms  # type: ignore
    ASE_AVAILABLE = True
except Exception:
    ASE_AVAILABLE = False


def _load_xyz(path: Path) -> Tuple[List[int], List[List[float]]]:
    """Load an .xyz file and return (atomic_numbers, positions).

    Uses ASE if available; otherwise falls back to a minimal parser.
    """
    t0 = _time.perf_counter()
    if ASE_AVAILABLE:
        atoms: Atoms = ase_read(str(path))  # type: ignore
        Z = atoms.get_atomic_numbers().tolist()
        pos = atoms.get_positions().tolist()
        dt = _time.perf_counter() - t0
        print(f"[timing] load_xyz (ASE) wall={dt:.4f}s", flush=True)
        return Z, pos

    text = path.read_text().splitlines()
    if len(text) < 3:
        raise ValueError("Invalid XYZ file: too short")
    nat = int(text[0].strip())
    lines = text[2: 2 + nat]
    PT = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
    }
    Z: List[int] = []
    pos: List[List[float]] = []
    for ln in lines:
        parts = ln.split()
        sym = parts[0]
        x, y, z = map(float, parts[1:4])
        Z.append(PT[sym])
        pos.append([x, y, z])
    dt = _time.perf_counter() - t0
    print(f"[timing] load_xyz (fallback) wall={dt:.4f}s", flush=True)
    return Z, pos


async def run_one_session(
    uri: str,
    Z: List[int],
    xyz: List[List[float]],
    duration_s: float,
    calculator: str,
    target_frames: Optional[int] = None,
    timeout_s: Optional[float] = None,
) -> int:
    """Open one WS session, run MD, and count frames.

    Uses protobuf frames by default; falls back to parsing JSON text frames
    if the server sends text.
    """
    frames = 0
    t0 = _time.perf_counter()
    deadline = (
        None
        if target_frames is not None
        else time.time() + float(duration_s)
    )
    last_seq = 0
    server_seq_seen = 0

    def _frame_mode_done() -> bool:
        if target_frames is None:
            return False
        if frames >= int(target_frames):
            return True
        if (
            timeout_s is not None
            and (_time.perf_counter() - t0) > float(timeout_s)
        ):
            return True
        return False

    async with websockets.connect(uri) as ws:
        # INIT_SYSTEM (protobuf)
        seq = 1
        init = pb.ClientAction()
        init.seq = seq
        init.type = pb.ClientAction.Type.INIT_SYSTEM
        init.atomic_numbers.extend([int(z) for z in Z])
        for p in xyz:
            v = pb.Vec3()
            v.v.extend([float(p[0]), float(p[1]), float(p[2])])
            init.positions.append(v)
        await ws.send(init.SerializeToString())

        # START_SIMULATION (protobuf)
        seq += 1
        start = pb.ClientAction()
        start.seq = seq
        start.type = pb.ClientAction.Type.START_SIMULATION
        start.simulation_type = pb.ClientAction.SimType.MD
        sp = pb.SimulationParams()
        sp.calculator = calculator
        sp.temperature = 298.0
        sp.timestep_fs = 1.0
        sp.friction = 0.02
        start.simulation_params.CopyFrom(sp)
        await ws.send(start.SerializeToString())

        # Main receive loop
        while True:
            if target_frames is None:
                if time.time() >= deadline:  # type: ignore[arg-type]
                    break
            else:
                if _frame_mode_done():
                    break
            try:
                msg = await asyncio.wait_for(ws.recv(), timeout=0.5)
            except asyncio.TimeoutError:
                # send periodic protobuf ack ping even on idle
                seq += 1
                ack = pb.ClientAction()
                ack.seq = seq
                ack.type = pb.ClientAction.Type.PING
                ack.ack = last_seq
                await ws.send(ack.SerializeToString())
                continue

            # report running FPS upon each received frame
            elapsed = _time.perf_counter() - t0
            if elapsed > 0:
                print(
                    (
                        f"[roy-fps] frames={frames} "
                        f"fps={frames/elapsed:.2f}"
                    ),
                    flush=True,
                )

            # Try binary protobuf first
            if isinstance(msg, (bytes, bytearray)):
                try:
                    res = pb.ServerResult()
                    res.ParseFromString(msg)
                    server_seq_seen = max(server_seq_seen, int(res.seq or 0))
                except Exception:
                    server_seq_seen += 1
            else:
                # Try JSON text
                try:
                    data = json.loads(msg)
                    server_seq_seen = max(
                        server_seq_seen, int(data.get("seq") or 0)
                    )
                except Exception:
                    server_seq_seen += 1

            frames += 1
            last_seq = server_seq_seen

            # ack so server backpressure window advances (protobuf ack)
            seq += 1
            ack = pb.ClientAction()
            ack.seq = seq
            ack.type = pb.ClientAction.Type.PING
            ack.ack = last_seq
            await ws.send(ack.SerializeToString())

        # STOP_SIMULATION
        seq += 1
        stop = pb.ClientAction()
        stop.seq = seq
        stop.type = pb.ClientAction.Type.STOP_SIMULATION
        await ws.send(stop.SerializeToString())

    dt = max(1e-9, _time.perf_counter() - t0)
    print(
        (
            f"[timing] run_one_session duration={duration_s}s "
            f"frames={frames} wall={dt:.4f}s"
        ),
        flush=True,
    )
    return frames


def main():
    ap = argparse.ArgumentParser(
        description="ROY MD steps/sec benchmark over WebSocket"
    )
    ap.add_argument(
        "--base",
        default="ws://127.0.0.1:8000",
        help="Base ws URL (no trailing /ws)",
    )
    ap.add_argument(
        "--file",
        default=str(Path("public/molecules/roy.xyz")),
        help="Path to ROY XYZ file",
    )
    ap.add_argument(
        "--duration",
        type=float,
        default=10.0,
        help="Benchmark duration in seconds",
    )
    ap.add_argument(
        "--duration-ms",
        type=float,
        default=None,
        help="Benchmark duration in milliseconds (overrides --duration)",
    )
    ap.add_argument(
        "--frames",
        type=int,
        default=None,
        help="Target frames to receive (overrides duration mode)",
    )
    ap.add_argument(
        "--timeout",
        type=float,
        default=5.0,
        help="Timeout seconds when using --frames mode",
    )
    ap.add_argument(
        "--calculator",
        choices=["lj", "uma"],
        default="lj",
        help="Calculator to request",
    )
    args = ap.parse_args()

    xyz_path = Path(args.file)
    if not xyz_path.exists():
        raise SystemExit(f"XYZ file not found: {xyz_path}")
    Z, xyz = _load_xyz(xyz_path)

    uri = args.base.rstrip("/") + "/ws"
    # Select duration from seconds or milliseconds flag
    duration_s = args.duration
    if args.duration_ms is not None:
        duration_s = float(args.duration_ms) / 1000.0

    t0 = _time.perf_counter()
    frames = asyncio.run(
        run_one_session(
            uri,
            Z,
            xyz,
            duration_s,
            args.calculator,
            target_frames=args.frames,
            timeout_s=args.timeout if args.frames is not None else None,
        )
    )
    dt = _time.perf_counter() - t0
    print(f"[timing] benchmark total wall={dt:.4f}s", flush=True)

    sps = frames / float(duration_s) if duration_s > 0 else 0.0
    print("ROY benchmark results:")
    print(
        (
            f"atoms= {len(Z)}  duration_s= {duration_s}  "
            f"calculator= {args.calculator}"
        )
    )
    print(f"frames= {frames}  steps_per_sec= {sps:.2f}")


if __name__ == "__main__":
    main()
