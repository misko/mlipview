# Examples

This directory contains optional example scripts not used by the web UI or Jest tests.

## `test_water.py`
Simple connectivity test against a running `fairchem_local_server/server.py` instance. Posts a small water molecule payload to the `/compute` endpoint and prints the response.

Usage:
```
source mlipview/bin/activate
python examples/test_water.py
```
Ensure the server is running first (default: `uvicorn fairchem_local_server.server:app --reload`).
