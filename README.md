# Project Overview

## Overview

This project provides a lightweight viewer + backend integration for running remote ML and classical (LJ) single-point calculations and fixed-step relaxations through HTTP endpoints. All client-side Lennard-Jones and BFGS optimizer code has been removed; relaxations are now performed exclusively by the backend.

Key FastAPI endpoints (served by `fairchem_local_server/server.py`):
* `POST /simple_calculate` – single-point energy / forces (and optional stress) for a given geometry.
* `POST /relax` – fixed number of BFGS steps returning final positions, forces, stress, and energies.
* `POST /md` – advance NVT MD a fixed number of steps (VelocityVerlet + per-step velocity rescale thermostat).

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
`/md` request body (advance 5 MD steps @298K):
```json
{
  "atomic_numbers": [8,1,1],
  "coordinates": [[0,0,0],[0.96,0,0],[-0.24,0.93,0]],
  "steps": 5,
  "temperature": 298.0,
  "timestep_fs": 1.0,
  "calculator": "uma"
}
```
Response (abridged):
```json
{
  "initial_energy": -75.90,
  "final_energy": -75.85,
  "positions": [[...],[...],[...]],
  "velocities": [[...],[...],[...]],
  "forces": [[...],[...],[...]],
  "steps_completed": 5,
  "temperature": 301.2,
  "calculator": "uma"
}
```
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
* MD uses a lightweight deterministic NVT scheme (VelocityVerlet + velocity rescale each step) chosen for interactive stability and minimal state.
* Switching calculators is a constant-time enum branch; adding more backends would extend the enum + selection logic.
* Client code fetches forces/energies remotely; no physics or optimizer logic remains in the browser bundle.

## Development
### API Endpoints (Post-Cleanup)

The viewer now exclusively uses Ray Serve namespaced endpoints:

Canonical:
  * `/serve/health`
  * `/serve/simple`  (single-point energy/forces)
  * `/serve/relax`   (BFGS relaxation steps)
  * `/serve/md`      (MD integration steps)

Previously supported root aliases (`/simple_calculate`, `/relax`, `/md`, `/calculate`, `/health`) have been removed.
If you had scripts using legacy paths, update them to the `/serve/*` equivalents.


Run type & lint checks (if configured) and the Jest tests:
```bash
npm test
```

### Dev Proxy /md

The frontend calls all backend endpoints (`/simple_calculate`, `/relax`, `/md`, `/health`) using relative paths so that in development Vite proxies them to the FastAPI process. If you receive a 404 for `/md` in the browser while other endpoints work, confirm:

1. `vite.config.js` contains a proxy entry for `/md` (it should look like `'/md': { target: backendTarget, changeOrigin: true, secure: false }`).
2. The FastAPI backend is actually running (check `curl http://127.0.0.1:8000/health`).
3. You're not hardcoding an absolute different origin in `window.__MLIPVIEW_SERVER` (leave it unset for proxy).
4. If using a custom hostname (e.g. `https://kalman:5173`), the name is listed in `allowedHosts` in `vite.config.js`.

Override backend target with env vars before launching Vite:
```bash
export FASTAPI_HOST=kalman
export FASTAPI_PORT=8000
npx vite
```
HTTPS dev certs (optional): place `localhost-key.pem` and `localhost-cert.pem` in repo root or set `NO_VITE_HTTPS=1` to force HTTP.

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
* `POST /md` (NVT MD) – body adds `steps`, `temperature` (K), optional `timestep_fs` (default 1.0)

### Frontend MD Controls

HUD buttons:
* `MD Step` – performs one 1 fs MD step via `/md`.
* `MD Run` – continuous MD loop (enable via `window.__MLIP_FEATURES.MD_LOOP = true`).
* `Relax Run` – continuous relaxation loop (auto-enabled by default; disable via `viewerApi.enableFeatureFlag('RELAX_LOOP', false)`).

Loops are auto-enabled on load. You can disable or re-enable them dynamically via the runtime API:
```js
// After viewerApi is created
// Disable
viewerApi.enableFeatureFlag('RELAX_LOOP', false);
viewerApi.enableFeatureFlag('MD_LOOP', false);
// Re-enable
viewerApi.enableFeatureFlag('RELAX_LOOP', true);
viewerApi.enableFeatureFlag('MD_LOOP', true);
```

Continuous run implementation details (Relax Run / MD Run):
* Makes at most ONE backend step request every 30 ms under normal operation (request pacing, configurable).
* On server/network/parse failures the loop engages exponential backoff starting at 200 ms doubling each failure up to 5 s, then resumes normal pacing after a successful step.
* Aborts automatically after 10 consecutive step errors (protects against persistent 500 responses).
* `startRelaxContinuous({ maxSteps })` returns `{ converged:boolean, steps:number }`.
* `startMDContinuous({ steps })` returns `{ completed:boolean, steps:number }` (steps actually executed; may be lower if aborted by error streak).
* Runtime pacing config stored at `window.__MLIP_CONFIG.minStepIntervalMs` (default 30).
  * Adjust with `viewerApi.setMinStepInterval(50)` (clamped to >=1ms).
* You can prematurely halt a running loop via `viewerApi.stopSimulation()` which clears the internal `running.kind` flag.

Testing:
* Added `tests/water_relax_run_parity.spec.js` – compares a 50-step browser continuous relax run to a reference aggregate `/relax` call (ASE backend) within a tight energy tolerance.
* Added `tests/water_md_run_stability.spec.js` – performs a 200-step continuous MD run (pacing/backoff active) ensuring energies remain finite and sufficient steps complete.

Console temperature override:
```js
window.__MD_TEMP = 500; // target 500 K for subsequent mdStep / runs
viewerApi.mdStep();
```
The backend returns instantaneous kinetic temperature; exposed as `viewerApi.state.dynamics.temperature`.
Safety: MD aborts server-side if any per-step displacement > 5 Å or coordinates become non-finite.

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
