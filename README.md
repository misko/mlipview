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

#### Temperature Slider (Desktop HUD)
The desktop HUD provides a temperature slider controlling the target thermostat temperature for MD steps and runs:

* Range: 0–3000 K
* Resolution: 30 discrete steps (internally mapped; includes an exact 298 K step)
* Visible tick labels: 298 K, 400 K, 3000 K (positioned proportionally by value; 0 K omitted to avoid overlap)
* Selected target stored at `window.__MLIP_TARGET_TEMPERATURE` (default 298) and mirrored to `viewerApi.state.dynamics.targetTemperature`.
* UI buttons (MD step / MD run) pass this target value; programmatic calls can still supply an explicit `temperature` parameter overriding the slider.

The instantaneous kinetic temperature reported by the backend remains available at `viewerApi.state.dynamics.temperature` and may differ from the target during early thermostat equilibration.

Live adjustment: During a continuous MD run started with `MD: (run)`, moving the slider now changes the target temperature of subsequent `/serve/md` step requests immediately. The loop re-samples `window.__MLIP_TARGET_TEMPERATURE` each iteration instead of closing over the initial value. Single-step MD (`MD: (step)`) still uses the temperature at the moment the button is pressed. Values are clamped to the slider domain [0, 2000] K as a safety measure.

Instantaneous Temperature Display: A HUD element labeled `T:` shows the instantaneous kinetic temperature returned by the backend after each MD step (single or continuous). This value may fluctuate around (and deviate temporarily from) the target slider temperature—especially early in a run or following abrupt geometry edits—because it is derived from current atomic velocities. The target (thermostat) temperature is the slider value; the instantaneous value is what the system actually exhibits at that timestep.

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

## 2025-10 UI Updates

Recent interface refinements (October 2025):
* Simplified HUD: removed Stop, Reset Energy, and Recompute Bonds buttons.
* Relax controls condensed to: `Relax: (step) (run)` where the run button toggles to `stop` while a continuous relaxation is active.
* MD controls now mirror Relax: `MD: (step) (run)` with identical run↔stop toggle behavior.
* Scene background (desktop & XR base) set to white for cleaner embedding and screenshots.
* XR dropdown initial label changed from `XR:none` to `SelectVR` for clearer action affordance.
* Atom dragging now respects a maximum radius: an atom cannot be moved farther than 5× its initial radial distance from the origin at drag start (prevents losing atoms far off screen).

These changes are non-breaking for programmatic APIs; only the HUD markup and initial XR label string changed.

## License

See repository root for license information.

### Precomputed Results Injection (Server Optimization)

Clients that already possess MLIP (or other calculator) results for the *initial* geometry can now avoid the first backend compute by supplying a `precomputed` object in `/relax` or `/md` requests.

Example:
```json
{
  "atomic_numbers": [8,1,1],
  "coordinates": [[0,0,0],[0.9575,0,0],[-0.2399872,0.92662721,0]],
  "steps": 1,
  "calculator": "lj",
  "precomputed": {
    "energy": -14.12345,
    "forces": [[0,0,0],[0,0,0],[0,0,0]],
    "stress": [0.1,0.2,0.3,0.01,0.02,0.03]
  }
}
```

Details:
* All fields optional; supply any subset of `energy`, `forces`, `stress`.
* `stress` may be Voigt length-6 `[xx, yy, zz, yz, xz, xy]` or flattened 3x3 length-9 row‑major (auto-converted to Voigt order).
* Provided `energy` becomes `initial_energy` (skips first calculator evaluation). If only forces are provided an energy compute still occurs when needed.
* Response includes `precomputed_applied` listing which items were injected.
* Subsequent relax/MD steps always recompute using the attached calculator.
* Non-finite or shape-mismatched inputs return HTTP 400.

Response excerpt:
```json
{
  "initial_energy": -14.12345,
  "final_energy": -14.13001,
  "forces": [[...],[...],[...]],
  "steps_completed": 1,
  "calculator": "lj",
  "precomputed_applied": ["energy","forces","stress"]
}
```

Use this when orchestrating single-step loops client-side to cut latency and redundant model evaluations.

### Frontend Automatic Precomputed Injection
The viewer now automatically attaches a `precomputed` block (energy, forces, stress) to `/serve/relax` and `/serve/md` requests when a fresh force cache is available for the current geometry. This avoids an initial backend calculator call for single-step operations. The cache is invalidated on any user geometry edit (drag, bond rotation, etc.). No user action is required—this is transparent and logged via `[API][relax][precomputed-attach]` or `[API][md][precomputed-attach]` debug messages when API debug is enabled.
