# Project Overview

## Overview

This project provides a lightweight viewer + backend integration for running remote ML and classical (LJ) single-point calculations and fixed-step relaxations through HTTP endpoints. All client-side Lennard-Jones and BFGS optimizer code has been removed; relaxations are now performed exclusively by the backend.

Key FastAPI endpoints (served by `fairchem_local_server/server.py`):
* `POST /simple_calculate` – single-point energy / forces (and optional stress) for a given geometry.
* `POST /relax` – fixed number of BFGS steps returning final positions, forces, stress, and energies.

Both endpoints accept a `calculator` field (`uma` | `lj`). UMA uses a cached FairChem MLIP; `lj` uses ASE's built-in Lennard-Jones calculator.

## Quick Start

Backend (Python / FastAPI):

```bash
uvicorn fairchem_local_server.server:app --host 127.0.0.1 --port 8000
```

Frontend / tests (Node):

```bash
npm install
npm test
```

## Endpoint Schemas (Simplified)

`/simple_calculate` request body:
```json
{
  "atomic_numbers": [8,1,1],
  "coordinates": [[0,0,0],[0.96,0,0],[ -0.24,0.93,0]],
  "properties": ["energy","forces"],
  "calculator": "uma"
}
```
Response:
```json
{ "results": { "energy": -76.123, "forces": [[...],[...],[...]] } }
```

`/relax` request body:
```json
{
  "atomic_numbers": [8,1,1],
  "coordinates": [[0,0,0],[0.96,0,0],[ -0.24,0.93,0]],
  "steps": 10,
  "calculator": "uma"
}
```
Response (abridged):
```json
{
  "initial_energy": -75.9,
  "final_energy": -76.12,
  "positions": [[...],[...],[...]],
  "forces": [[...],[...],[...]],
  "stress": null,
  "steps_completed": 10,
  "calculator": "uma"
}
```

## Design Notes

* UMA calculator objects are cached per composition to amortize graph construction.
* Relaxations run a fixed number of BFGS steps; early convergence (fmax) is logged but not yet used to stop.
* Switching calculators is a constant-time enum branch; adding more backends would extend the enum + selection logic.
* Client code fetches forces/energies remotely; no physics or optimizer logic remains in the browser bundle.

## Development

Run type & lint checks (if configured) and the Jest tests:
```bash
npm test
```

### Backend ML API (Python)

Relaxation, forces and energies are provided by FastAPI app `fairchem_local_server/server.py`.

Start it manually for local hacking (if not relying on Jest global setup):

```bash
source mlipview_venv/bin/activate  # or set MLIPVIEW_PYTHON to another env
uvicorn fairchem_local_server.server:app --host 0.0.0.0 --port 8000
```

Key endpoints:
* `POST /simple_calculate` (single point) – body: `{ atomic_numbers, coordinates, properties:["energy","forces"], calculator:"uma"|"lj" }`
* `POST /relax` (fixed BFGS steps) – body adds `steps` and optional `fmax` (currently informational)

### Integration Test (Water Relax Parity)

Jest global setup now launches BOTH servers automatically:
1. Node static server (port 4000)
2. Python uvicorn ML server (port 8000)

Override Python interpreter with `MLIPVIEW_PYTHON=/path/to/python`.

Run just the parity test:
```bash
npm test -- -t "water relax step energy parity"
```
It prints initial, relaxed and reference energies; ensures relaxed energy is not higher and reference matches the initial single-point.

### Removal of Legacy Client Optimizer

File `public/lj_bfgs.js` (in-browser LJ potential + BFGS) was removed. All relaxations must go through `/relax`.

Regenerate a relaxation trace directly via Python for debugging:
```bash
python fairchem_local_server/relax_water_http.py
```

### Playwright E2E Tests

UI tests in `tests-e2e/` are run with Playwright (`npx playwright test`). Global setup starts Node + Python servers only if their ports are free. Two added relax sanity tests:

* `waterDirectRelax.spec.js` – calls `viewerApi.relaxStep()` directly and logs before/after energy.
* `uiRelaxParity.spec.js` – simulates a user click on the Relax Step button.

Run just those:
```bash
npx playwright test tests-e2e/waterDirectRelax.spec.js tests-e2e/uiRelaxParity.spec.js
```

Use a custom Python interpreter:
```bash
MLIPVIEW_PYTHON=/path/to/python npx playwright test
```

## Future Work

* Add early-stop convergence to `/relax` (respect `fmax`).
* Streaming or chunked relaxation progress (Server-Sent Events / WebSocket).
* Additional calculators (e.g., EAM, GAP) behind the same abstraction.
* Optional stress-aware relaxations with cell degrees of freedom.

## License

See repository root for license information.
