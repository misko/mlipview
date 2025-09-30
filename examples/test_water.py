"""Example script to validate local fairchem server with a simple
water molecule.

Usage:
    source mlipview/bin/activate  # if not already active
    python examples/test_water.py

Requirements:
    - fairchem_local_server/server.py running (uvicorn)
    - Network accessible endpoint (default http://127.0.0.1:8000)

This script kept as an example; not used by the web frontend nor
automated tests.
"""

import http.client
import json

HOST = "127.0.0.1"
PORT = 8000
PATH = "/compute"

payload = {
    "atoms": [
        {"symbol": "O", "x": 0.0, "y": 0.0, "z": 0.0},
        {"symbol": "H", "x": 0.96, "y": 0.0, "z": 0.0},
        {"symbol": "H", "x": -0.24, "y": 0.93, "z": 0.0},
    ]
}

body = json.dumps(payload)
headers = {"Content-Type": "application/json"}

conn = http.client.HTTPConnection(HOST, PORT, timeout=10)
try:
    conn.request("POST", PATH, body=body, headers=headers)
    resp = conn.getresponse()
    data = resp.read().decode()
    print("Status:", resp.status, resp.reason)
    print("Raw Response:", data)
    try:
        j = json.loads(data)
        print("Parsed keys:", list(j.keys()))
    except Exception as e:
        print("Failed to parse JSON:", e)
finally:
    conn.close()
