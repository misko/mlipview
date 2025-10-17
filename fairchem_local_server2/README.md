# fairchem_local_server2 (WebSocket server)

Single-connection WebSocket API for UMA+ASE simulations.

## Endpoints
- `GET /serve/health`: health snapshot
- `WS /ws`: single WebSocket per session

## Message protocol
Default: JSON objects matching Pydantic models in `session_models.py`.
Optional: Protobuf frames if client sends binary messages using `session.proto`.
The server auto-detects: binary → protobuf, text → JSON.

ClientAction fields:
- `seq` (int), `ack` (int | null), `type` ("init_system" | "update_positions" | "start_simulation" | "stop_simulation" | "ping")
- `atomic_numbers` (N), `positions` (N×3), `velocities` (N×3, optional), `cell` (3×3, optional)
- `simulation_type` ("md" | "relax"), `simulation_params`

ServerResult fields:
- `seq` (int), `client_seq` (int), `positions` (N×3), `velocities` (N×3), `forces` (N×3), `cell` (3×3), `message` (optional)

Backpressure: server will not compute more than 10 replies beyond last client ack.

## Protobuf codegen
- Schema: `fairchem_local_server2/session.proto`
- Generate Python modules (example):

```
python -m pip install protobuf
python -m grpc_tools.protoc \
  -I fairchem_local_server2 \
  --python_out=fairchem_local_server2 \
  fairchem_local_server2/session.proto
```

This creates `fairchem_local_server2/session_pb2.py`. The server will use it automatically if available.

## Deploy
Use the entrypoint:

```
python -m fairchem_local_server2.serve_ws_app --ngpus 1 --ncpus 8 --nhttp 1
```

- UMA predictor replicas: `--ngpus` (one per GPU)
- ASE worker pool size: `--ncpus`
- WS ingress replicas: `--nhttp` (default 1)

## Test
Run the parity test (uses Lennard-Jones; no GPU required):

```
pytest -q fairchem_local_server2/tests_py/test_md_parity.py
```
