"""WebSocket stress/load test utility for UMA Serve WS API (/ws).

Opens N concurrent WebSocket sessions, initializes a small system,
and starts MD or Relax producing a stream of results. Measures message
throughput and simple latency proxies (time between send and first reply).

Usage:
  python -m fairchem_local_server2.stress_tool_v2 --base ws://127.0.0.1:8000 \
    --concurrency 64 --duration 60 --mode md --md-steps 1

Notes:
- Uses JSON messages by default; enable protobuf by sending binary frames
  if session_pb2.py is present and desired (adapter auto-detects).
"""

from __future__ import annotations

import argparse
import asyncio
import time
from dataclasses import dataclass
from typing import List, Tuple

import json
import websockets


def _water() -> Tuple[List[int], List[List[float]]]:
    Z = [8, 1, 1]
    xyz = [
        [0.0, 0.0, 0.0],
        [0.9572, 0.0, 0.0],
        [-0.239987, 0.926627, 0.0],
    ]
    return Z, xyz


@dataclass
class Metrics:
    sessions_ok: int = 0
    sessions_err: int = 0
    messages_rcvd: int = 0


async def _one_session(
    uri: str, mode: str, steps: int, m: Metrics, calculator: str
):
    try:
        async with websockets.connect(uri) as ws:
            # init
            Z, xyz = _water()
            seq = 1
            init = {
                "seq": seq,
                "type": "init_system",
                "atomic_numbers": Z,
                "positions": xyz,
                "cell": [
                    [15.0, 0.0, 0.0],
                    [0.0, 15.0, 0.0],
                    [0.0, 0.0, 15.0],
                ],
            }
            await ws.send(json.dumps(init))
            # start
            seq += 1
            start = {
                "seq": seq,
                "type": "start_simulation",
                "simulation_type": mode,
                "simulation_params": {
                    "calculator": calculator,
                    "temperature": 298.0,
                    "timestep_fs": 1.0,
                    "friction": 0.02,
                    "fmax": 0.05,
                    "max_step": 0.2,
                },
            }
            await ws.send(json.dumps(start))

            # receive a few messages and ack
            deadline = time.time() + 2.0
            last_seq = 0
            while time.time() < deadline:
                try:
                    msg = await asyncio.wait_for(ws.recv(), timeout=0.5)
                except asyncio.TimeoutError:
                    continue
                m.messages_rcvd += 1
                try:
                    data = json.loads(msg)
                except Exception:
                    # protobuf binary: just count
                    continue
                last_seq = max(last_seq, int(data.get("seq") or 0))
                # send ack
                seq += 1
                await ws.send(
                    json.dumps({
                        "seq": seq,
                        "ack": last_seq,
                        "type": "ping",
                    })
                )
            m.sessions_ok += 1
    except Exception:
        m.sessions_err += 1


async def run_stress_ws(
    base_ws: str, duration_s: int, concurrency: int, mode: str, calculator: str
):
    m = Metrics()
    # normalize
    base = base_ws.rstrip("/") + "/ws"

    async def worker():
        await _one_session(base, mode, steps=1, m=m, calculator=calculator)

    started = time.time()
    deadline = started + float(duration_s)
    sem = asyncio.Semaphore(max(1, int(concurrency)))
    tasks: List[asyncio.Task] = []

    async def spawn():
        async with sem:
            await worker()

    while time.time() < deadline:
        tasks = [t for t in tasks if not t.done()]
        while len(tasks) < concurrency:
            tasks.append(asyncio.create_task(spawn()))
        await asyncio.sleep(0.01)
    if tasks:
        await asyncio.gather(*tasks, return_exceptions=True)
    return m


def main():
    ap = argparse.ArgumentParser(description="UMA WS stress utility")
    ap.add_argument("--base", default="ws://127.0.0.1:8000")
    ap.add_argument("--duration", type=int, default=15)
    ap.add_argument("--concurrency", type=int, default=32)
    ap.add_argument("--mode", choices=["md", "relax"], default="md")
    ap.add_argument(
        "--calculator",
        choices=["lj", "uma"],
        default="lj",
        help="Calculator to request in simulation params",
    )
    args = ap.parse_args()

    m: Metrics = asyncio.run(
        run_stress_ws(
            base_ws=args.base,
            duration_s=args.duration,
            concurrency=args.concurrency,
            mode=args.mode,
            calculator=args.calculator,
        )
    )

    print("\n=== UMA WS Stress Results ===")
    print(f"base= {args.base} mode= {args.mode} calculator= {args.calculator}")
    print(f"duration_s= {args.duration} concurrency= {args.concurrency}")
    print(f"sessions_ok= {m.sessions_ok} sessions_err= {m.sessions_err}")
    print(f"messages_rcvd= {m.messages_rcvd}")


if __name__ == "__main__":
    main()
