from __future__ import annotations

import asyncio
import json
import sys

import aiohttp
import websockets


async def main(base: str = "http://127.0.0.1:8000"):
    # reset stats
    async with aiohttp.ClientSession() as sess:
        try:
            async with sess.post(base.rstrip("/") + "/uma/reset") as r:
                _ = await r.json()
        except Exception as e:
            print("reset error:", e)
            # continue anyway
    # drive one UMA session
    ws_url = base.replace("http", "ws").rstrip("/") + "/ws"
    try:
        async with websockets.connect(ws_url) as ws:
            Z = [8, 1, 1]
            xyz = [
                [0.0, 0.0, 0.0],
                [0.9572, 0.0, 0.0],
                [-0.239987, 0.926627, 0.0],
            ]
            await ws.send(
                json.dumps(
                    {
                        "seq": 1,
                        "type": "init_system",
                        "atomic_numbers": Z,
                        "positions": xyz,
                        "cell": [
                            [15.0, 0.0, 0.0],
                            [0.0, 15.0, 0.0],
                            [0.0, 0.0, 15.0],
                        ],
                    }
                )
            )
            await ws.send(
                json.dumps(
                    {
                        "seq": 2,
                        "type": "start_simulation",
                        "simulation_type": "md",
                        "simulation_params": {
                            "calculator": "uma",
                            "temperature": 298.0,
                            "timestep_fs": 1.0,
                            "friction": 0.02,
                        },
                    }
                )
            )
            # receive a couple messages
            for _ in range(2):
                msg = await ws.recv()
                if isinstance(msg, (bytes, bytearray)):
                    print("got binary frame", len(msg))
                else:
                    print("got text frame", str(msg)[:120])
                await ws.send(
                    json.dumps({"seq": 3, "type": "ping", "ack": 10})
                )
    except Exception as e:
        print("ws error:", e)
        return 1
    # read stats
    async with aiohttp.ClientSession() as sess:
        try:
            async with sess.get(base.rstrip("/") + "/uma/stats") as r:
                js = await r.json()
                print("UMA stats:", js)
        except Exception as e:
            print("stats error:", e)
            return 1
    return 0


if __name__ == "__main__":
    base = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8000"
    rc = asyncio.run(main(base))
    sys.exit(rc)
