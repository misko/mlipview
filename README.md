# Minimal LJ + BFGS Parity

This repository has been reduced to the essentials needed to reproduce an ASE Lennard-Jones BFGS relaxation for a water molecule in the browser (Node/Jest environment here). The goal is parity with the provided Python/ASE ground truth:

```
Initial energy:  2.802523
Final energy:   -2.983548  (BFGS, fmax 0.05)
```

## Implementation

File: `public/lj_bfgs.js`

Contains:
* `ljEnergyForces(positions)` — Lennard-Jones energy + forces with ASE-like defaults (epsilon=1, sigma=1, rc=3*sigma, energy shift so U(rc)=0).
* `bfgsOptimize({...})` — Dense BFGS with Armijo backtracking and max per-component step control.
* `optimizeWaterExample()` convenience helper.

## Parity Test

File: `tests/lj_bfgs_parity.spec.js` asserts:
* Initial energy within 5e-3 of 2.802523
* Relaxed energy within 5e-3 of -2.983548 with fmax < 0.05

Run it:

```
npm test --silent -- lj_bfgs_parity.spec.js
```

## Usage Example

```js
import { ljEnergyForces, bfgsOptimize } from './public/lj_bfgs.js';
const d = 0.9575; const t = Math.PI/180*104.51;
const positions = [ [d,0,0], [d*Math.cos(t), d*Math.sin(t), 0], [0,0,0] ];
console.log('Initial', ljEnergyForces(positions).energy);
const res = await bfgsOptimize({ positions, fmax:0.05, maxSteps:500, compute: async p => ljEnergyForces(p) });
console.log('Final', res.energy, 'steps', res.steps);
```

## Notes

* No neighbor list is needed for this tiny system.
* Energy shifting matches ASE's default (non-smooth) cutoff handling.
* Forces are analytical; dense BFGS is fine for very small atom counts.

## Future Extensions (Out of Scope Now)

* Add element-specific epsilon/sigma with Lorentz-Berthelot mixing.
* Optional smooth cutoff (ro, rc) parity with ASE `smooth=True`.
* Serialization of optimization trace for visualization.

---
Everything else (viewer, VR, remote force providers) was removed per request to focus solely on BFGS + LJ parity.

### Extending

Add a new backend by editing `public/forcefield/registry.js`:

```
FACTORIES.mybackend = (opts) => createMyBackendForceField(opts);
```

Then load with `?ff=mybackend`.

### Raw vs Scene Compute

`compute()` returns Babylon `Vector3` forces for immediate visualization.
`computeRaw({Z, xyz})` is backend-neutral numeric data (no Babylon dependency) and can be used for batch or headless analysis.

## Local UMA FAIR-Chem Server (Simplified)

When developing locally over `http://` the FAIR-Chem adapter talks to `http://localhost:8000` (default UMA). If you access the viewer itself over HTTPS, use an external reverse proxy or tunnel (e.g. Caddy, Traefik, nginx, ngrok) to terminate TLS and forward to the Node server + UMA backend. The previous dual HTTP+HTTPS helper script has been removed to reduce maintenance.

Mixed-content avoidance is now handled by:
1. Same-origin proxy endpoints in `server.js` (`/simple_calculate`, `/calculate`).
2. Automatic upgrade attempt inside `public/forcefield/fairchem.js` (logs a single warning on fallback).

If you need HTTPS for VR quickly:
```
caddy reverse-proxy --from https://localhost:8443 --to :3000
```
(Accept the local trust prompt / add a dev certificate.)

Production-grade deployment should still use a hardened reverse proxy with proper certificates.


## Test Pyramid Strategy

The test suite is layered for speed and fidelity:

1. Core (Jest / Node) – logic & state (`tests`, `mlipviewer2/tests`). Uses lightweight Babylon mocks.
2. Browser (Jest / jsdom) – DOM integration / highlight & canvas behavior (`tests-browser`). Adds minimal canvas + pointer event polyfills.
3. Engine (NullEngine harness) – optional Babylon scene graph checks without real GL (`tests/engine`).
4. Future (Playwright) – real headless browser for WebGL rendering regression (bonds present, no stray primitives). Not yet added, but structure is prepared.

Run all:
```
npm test
```
Only browser project:
```
npx jest --selectProjects browser
```
Only core project:
```
npx jest --selectProjects core
```
Engine-only sample:
```
npx jest tests/engine/nullEngine.spec.js --selectProjects core
```

Adding a browser test: drop a spec in `tests-browser/*.spec.js` and it will run under `jsdom` automatically.

## Static Build (Vite)

You can produce a fully static production bundle (HTML + ES modules + assets) using Vite.

Scripts added:

```
npm run build    # Generates dist/ using Vite (root = public/)
npm run preview  # Serves the built dist/ locally on port 5174
```

Build output goes to `dist/` and can be deployed behind any static file host or CDN (Netlify, GitHub Pages, nginx, S3 + CloudFront, etc.). Requires Node 20.19+ (see `.nvmrc`).

Notes:
1. External Babylon CDN scripts referenced in `public/index.html` are preserved as-is (they are not bundled). This keeps the bundle lean and leverages browser caching. For full offline packaging you can vendor these scripts by removing the `<script src=...>` tags and importing modules instead.
2. Environment-specific settings (e.g. remote FAIR-Chem endpoints) should be passed via query params or a small runtime config file; no `.env` integration has been added yet.
3. If you only need a static viewer (no Node server API routes), you can deploy just the contents of `dist/`. The express server remains useful for local dev and any dynamic endpoints (molecule listing, force proxy, etc.).
4. Source maps are enabled for easier debugging of production issues (`dist/**/*.js.map`).

Minimal end-to-end check:

```
npm install
npm run build
npm run preview
# open http://localhost:5174
```

If you see the viewer and can load the default molecule, the static build is functioning.

## FairChem HTTP Reference Generation

To generate and update the FairChem 20-step BFGS trace for water:

1. Start the local FastAPI server:
	```bash
	uvicorn fairchem_local_server.server:app --host 127.0.0.1 --port 8000
	```
2. Run the relax script:
	```bash
	python fairchem_local_server/relax_water_http.py
	```
3. This writes `public/reference/water_fairchem_bfgs_trace.json`.
4. Run parity tests:
	```bash
	npm test -- fairchem_bfgs_parity.spec.js
	npx playwright test tests-e2e/waterFairChemRelaxation.spec.js
	```

The JS + browser parity tests skip gracefully if the server is unreachable or the reference file is empty.

