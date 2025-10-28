from __future__ import annotations

import argparse
import asyncio
import json
import time
from pathlib import Path
from typing import List

import websockets

# Reuse XYZ loader from benchmark if available, else include a minimal one
try:
    from scripts.roy_md_benchmark import _load_xyz  # type: ignore
except Exception:
    from ase.io import read as ase_read  # type: ignore

    def _load_xyz(path: Path):  # type: ignore
        atoms = ase_read(str(path))
        Z = atoms.get_atomic_numbers().tolist()
        pos = atoms.get_positions().tolist()
        return Z, pos


async def receive_n_frames(uri: str, Z: List[int], xyz: List[List[float]], n: int, calculator: str, timeout_s: float = 20.0) -> int:
    """Open a WS session and return once n frames have been received, or timeout.

    Returns the number of frames received (>= n indicates success).
    """
    async with websockets.connect(uri) as ws:
        seq = 1
        init = {
            "seq": seq,
            "type": "init_system",
            "atomic_numbers": Z,
            "positions": xyz,
            "cell": None,
        }
        await ws.send(json.dumps(init))
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

        server_seq_seen = 0
        frames = 0
        deadline = time.time() + float(timeout_s)
        while time.time() < deadline and frames < n:
            try:
                msg = await asyncio.wait_for(ws.recv(), timeout=0.8)
            except asyncio.TimeoutError:
                # ping keepalive + advance ack
                seq += 1
                await ws.send(json.dumps({"seq": seq, "type": "ping", "ack": server_seq_seen}))
                continue
            frames += 1
            # track server seq for ACK even if binary frames
            try:
                data = json.loads(msg)
                server_seq_seen = max(server_seq_seen, int(data.get("seq") or 0))
            except Exception:
                server_seq_seen += 1
            # send ack for backpressure window
            seq += 1
            await ws.send(json.dumps({"seq": seq, "type": "ping", "ack": server_seq_seen}))
        return frames


def main():
    ap = argparse.ArgumentParser(description="ROY WS 20-frames smoke test")
    ap.add_argument("--base", default="ws://127.0.0.1:8000", help="Base ws URL (no trailing /ws)")
    ap.add_argument("--file", default=str(Path("public/molecules/roy.xyz")), help="Path to ROY XYZ file")
    ap.add_argument("--calculator", choices=["lj", "uma"], default="lj")
    ap.add_argument("--frames", type=int, default=20, help="Frames to wait for before success")
    ap.add_argument("--timeout", type=float, default=20.0, help="Timeout seconds")
    args = ap.parse_args()

    Z, xyz = _load_xyz(Path(args.file))
    uri = args.base.rstrip("/") + "/ws"

    frames = asyncio.run(receive_n_frames(uri, Z, xyz, args.frames, args.calculator, args.timeout))
    ok = frames >= args.frames
    print(f"received={frames} target={args.frames} ok={ok}")
    raise SystemExit(0 if ok else 2)


if __name__ == "__main__":
    main()
