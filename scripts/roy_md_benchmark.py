from __future__ import annotations

import argparse
import asyncio
import time
from pathlib import Path
from typing import Tuple, List

import json
import websockets

try:
    # Prefer ASE for robust XYZ parsing if available
    from ase.io import read as ase_read  # type: ignore
    from ase.atoms import Atoms  # type: ignore
    ASE_AVAILABLE = True
except Exception:
    ASE_AVAILABLE = False


def _load_xyz(path: Path) -> Tuple[List[int], List[List[float]]]:
    """Load an .xyz file and return (atomic_numbers, positions).

    Uses ASE if available; otherwise falls back to a minimal parser that
    expects standard XYZ format: first line natoms, second line comment,
    following lines 'Elem x y z'.
    """
    if ASE_AVAILABLE:
        atoms: Atoms = ase_read(str(path))  # type: ignore
        Z = atoms.get_atomic_numbers().tolist()
        pos = atoms.get_positions().tolist()
        return Z, pos

    # Minimal fallback parser
    text = path.read_text().splitlines()
    if len(text) < 3:
        raise ValueError("Invalid XYZ file: too short")
    try:
        nat = int(text[0].strip())
    except Exception as e:
        raise ValueError(f"Invalid XYZ first line: {e}")
    lines = text[2 : 2 + nat]
    if len(lines) < nat:
        raise ValueError("Invalid XYZ: not enough coordinate lines")

    # Periodic table mapping
    PT = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
        'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
        'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    }

    Z: List[int] = []
    pos: List[List[float]] = []
    for ln in lines:
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"Bad XYZ line: {ln}")
        sym = parts[0]
        try:
            x, y, z = map(float, parts[1:4])
        except Exception as e:
            raise ValueError(f"Bad XYZ coords: {e}")
        if sym not in PT:
            raise ValueError(f"Unknown element symbol: {sym}")
        Z.append(PT[sym])
        pos.append([x, y, z])
    return Z, pos


async def run_one_session(uri: str, Z: List[int], xyz: List[List[float]], duration_s: float, calculator: str) -> int:
    """Open one WS session, run MD for duration_s, and count frames received.

    Returns the number of result frames (1 result ~= 1 MD step in our server).
    """
    frames = 0
    async with websockets.connect(uri) as ws:
        # init
        seq = 1
        init = {
            "seq": seq,
            "type": "init_system",
            "atomic_numbers": Z,
            "positions": xyz,
            "cell": None,
        }
        await ws.send(json.dumps(init))
        # start MD
        seq += 1
        start = {
            "seq": seq,
            "type": "start_simulation",
            "simulation_type": "md",
            "simulation_params": {
                "calculator": calculator,
                "temperature": 298.0,
                "timestep_fs": 1.0,
                "friction": 0.02,
            },
        }
        await ws.send(json.dumps(start))

        deadline = time.time() + float(duration_s)
        last_seq = 0
        while time.time() < deadline:
            try:
                msg = await asyncio.wait_for(ws.recv(), timeout=0.5)
            except asyncio.TimeoutError:
                # send periodic ack ping even on idle
                seq += 1
                await ws.send(json.dumps({"seq": seq, "type": "ping", "ack": last_seq}))
                continue
            frames += 1
            # try to parse JSON to extract seq; if binary (protobuf), just count
            try:
                data = json.loads(msg)
                last_seq = max(last_seq, int(data.get("seq") or 0))
            except Exception:
                # binary frame
                pass
            # ack so server backpressure window advances
            seq += 1
            await ws.send(json.dumps({"seq": seq, "type": "ping", "ack": last_seq}))
    return frames


def main():
    ap = argparse.ArgumentParser(description="ROY MD steps/sec benchmark over WebSocket")
    ap.add_argument("--base", default="ws://127.0.0.1:8000", help="Base ws URL (no trailing /ws)")
    ap.add_argument("--file", default=str(Path("public/molecules/roy.xyz")), help="Path to ROY XYZ file")
    ap.add_argument("--duration", type=float, default=10.0, help="Benchmark duration in seconds")
    ap.add_argument("--calculator", choices=["lj", "uma"], default="lj", help="Calculator to request")
    args = ap.parse_args()

    xyz_path = Path(args.file)
    if not xyz_path.exists():
        raise SystemExit(f"XYZ file not found: {xyz_path}")
    Z, xyz = _load_xyz(xyz_path)

    uri = args.base.rstrip("/") + "/ws"
    frames = asyncio.run(run_one_session(uri, Z, xyz, args.duration, args.calculator))

    sps = frames / float(args.duration) if args.duration > 0 else 0.0
    print("ROY benchmark results:")
    print(f"atoms= {len(Z)}  duration_s= {args.duration}  calculator= {args.calculator}")
    print(f"frames= {frames}  steps_per_sec= {sps:.2f}")


if __name__ == "__main__":
    main()
