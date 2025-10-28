# Examples

This directory contains optional example scripts that are not used by the web UI or automated tests. They are handy for quick smoketests when validating environments outside the full websocket stack.

## `test_water.py`
Simple connectivity test against the legacy JSON REST endpoint (`fairchem_local_server/server.py`). It posts a small water molecule payload to `/compute` and prints the response payload. The main MLIP viewer uses the protobuf WebSocket server (`fairchem_local_server2/ws_app.py`), so this script is only needed when sanity-checking the older REST interface.

Usage:
```
source mlipview_venv/bin/activate
python examples/test_water.py
```
Ensure the REST server is running first (default: `uvicorn fairchem_local_server.server:app --reload`). For the WebSocket timeline experience, refer to the instructions in `README.md`.
