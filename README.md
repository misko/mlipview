## MLIPView: Desktop Molecular Viewer + Remote MD/Relax

This repository provides a lightweight, fast molecular viewer wired to a Python backend for ML/classical forces, relaxations, and MD. The front end renders atoms/bonds/forces with Babylon.js and drives short “step” APIs on the backend. The desktop flow prioritizes responsiveness with robust pacing, backoff, and safe partial-apply of results when the geometry changes mid-flight.

Highlights:
- Thin‑instance rendering for atoms/bonds/forces (Babylon.js) with selection highlights
- Event-driven architecture separating domain state, rendering, UI, and picking
- Remote API for forces/relax/MD via FastAPI (Ray Serve-style endpoints)
- Resilient loops with exponential backoff and request pacing
- Comprehensive desktop UI including selection, PBC cell controls, MD temperature/friction, uploads, and SMILES


## Quick start

Prereqs:
- Node.js >= 20.19
- Python 3.10+ with FastAPI/uvicorn for the backend (see `fairchem_local_server/`)

Install and run:
- Backend (Python):
  - Activate your env (optional): `source mlipview_venv/bin/activate`
  - Start FastAPI: `uvicorn fairchem_local_server.server:app --host 127.0.0.1 --port 8000`
- Frontend (Node):
  - `npm install`
  - Dev server + static: `npm run dev` (serves the web app on port 4000 by default)

Open http://localhost:4000 and load a molecule (left panel). The desktop HUD includes Relax and MD controls; see MD flow below.

Notes:
- HTTPS dev certs: place `localhost-key.pem` and `localhost-cert.pem` in repo root to enable HTTPS on 4443.
- Tests: `npm test` (Jest + jsdom). Playwright E2E available in `tests-e2e/`.


## Backend API (served by FastAPI)

Canonical endpoints (Ray Serve style):
- `POST /serve/health`
- `POST /serve/simple`  – single-point energy/forces (and optional stress)
- `POST /serve/relax`   – fixed-step relax; returns final positions/forces/stress
- `POST /serve/md`      – NVT MD for N steps; returns updated positions/velocities/forces/temperature

The frontend calls relative URLs; in dev, requests are routed by the Node server. To explicitly target a different origin, set `window.__MLIPVIEW_SERVER` in the console before loading.


## Desktop MD flow: request → parse → render → next request

Where in the code (front end):
- Entry/runtime: `public/index.js` (initNewViewer, startMDContinuous/startRelaxContinuous)
- API paths: `public/api_endpoints.js`
- State and events: `public/domain/moleculeState.js`, `public/domain/eventBus.js`
- Rendering: `public/render/scene.js`, `public/render/moleculeView.js`
- Picking/interaction: `public/core/pickingService.js`, `public/domain/manipulationService.js`, `public/domain/selectionService.js`
- UI (desktop panel): `public/ui/desktopPanel.js` with helpers in `public/ui/*.js`

End-to-end MD step (happy path):
1) Trigger: User presses “MD: run” or “MD: step” in the desktop HUD (wired in `public/index.html` via `buildDesktopPanel()` and button handlers). Programmatic calls use `viewerApi.mdStep()` or `viewerApi.startMDContinuous()` from `public/index.js`.
2) Request formation: `public/index.js` collects the current elements, positions, and optional cell from `moleculeState`. For continuous MD, it uses steps=1 and the current target temperature from the temperature slider. If a fresh force cache exists for this exact geometry, a `precomputed` block is attached to save a redundant backend compute.
3) POST /serve/md: Sent to the backend. Request pacing guarantees at most one request every N ms (default 30ms, configurable). If the previous request is still in flight, the loop naturally waits.
4) Parse: On a successful response, `public/index.js` updates `molState.positions`, `molState.dynamics.forces`, velocity/temperature fields, and caches forces/energy for potential `precomputed` usage on the next step.
5) Render: `molState.markPositionsChanged()` emits `positionsChanged`; `moleculeView.js` receives events and updates thin instances (atoms/bonds/forces) and selection highlight meshes. Any bond-length-dependent visuals update automatically. The temperature display and energy plot are also updated.
6) Next request: For “run”, the loop schedules the next /serve/md step honoring the min-step interval. For “step”, it stops after one iteration.

Resilience and safety:
- Exponential backoff after failures (starts at ~200ms, doubles up to ~5s), then resumes normal pacing after success
- Auto-abort after a streak of failures to prevent log spam
- Partial-apply guards: each response carries a “version” stamped at request time; if the user edited atoms (drag/rotate) while the request was in flight, the viewer discards stale parts and only applies safe fields (e.g., temperature), avoiding geometry teleportation
- Force cache is invalidated on user edits

Relax flow is identical, but calls `POST /serve/relax` and applies returned positions/forces in batches of steps.


## Supported desktop interactions and code locations

Selection & picking
- Click atom: selects atom, updates Selection panel and highlights (index/color) – `public/core/pickingService.js`, `public/domain/selectionService.js`, `public/render/moleculeView.js`
- Click bond: toggles bond selection and orientation metadata – `public/selection-model.js`, `public/domain/selectionService.js`

Atom manipulation
- Drag atom: click-drag moves an atom on a screen-aligned plane; camera detaches while dragging. On release, bonds are recomputed if needed – `public/core/pickingService.js`, `public/domain/manipulationService.js`
- Rotate selected bond: keyboard/UI/gesture triggers pass delta angle to `manipulationService.rotateBond(delta)` – `public/domain/manipulationService.js`

Rendering & overlays
- Thin-instance atoms/bonds/forces; selection highlight meshes parented to masters – `public/render/moleculeView.js`
- Cell lines and ghost atoms/bonds (PBC visualization) – `public/render/moleculeView.js`, `public/util/pbc.js`

Periodic boundary conditions
- Toggle cell visibility, edit parameters, wrap positions – `public/ui/cellToggle.js`, `public/util/pbc.js`, `public/domain/moleculeState.js`

Forces, MD & relax controls
- Toggle force vectors, MD temperature & friction sliders, run/step buttons – `public/ui/*.js`, `public/ui/desktopPanel.js`, wired into `public/index.js`

Loaders & inputs
- Load default or example molecules, parse XYZ, fetch SMILES via PubChem – `public/util/moleculeLoader.js`, `public/util/xyzLoader.js`, `public/util/smilesLoader.js`

Touch & mobile shims
- Touch gestures mapped to camera orbit/zoom and picking – `public/ui/touchControls.js`

VR (optional)
- Session setup and controller interactions (laser pick, rotate/drag, two-hand scaling) – `public/vr/setup.js`


## Running, building, testing

Local servers
- Node static/dev server: `npm run dev` (HTTP 4000; HTTPS 4443 when local certs exist)
- Python FastAPI backend: `uvicorn fairchem_local_server.server:app --host 127.0.0.1 --port 8000`

Scripts (package.json)
- `npm start` – start Node server
- `npm run dev` – development mode
- `npm run build` – Vite build
- `npm run preview` – Vite preview
- `npm test` – run Jest tests (unit + jsdom based)

Troubleshooting
- 404s for /serve/*: ensure the FastAPI server is running and reachable; confirm you didn’t override the base URL unexpectedly
- Render looks empty: check browser console for errors; try loading an example via the left panel; confirm `positionsChanged` events fire when editing
- MD loop not advancing: verify min-step interval and backoff in console logs; ensure the backend returns finite positions/forces


## Where to look in the code (jump list)

- Entry & orchestration: `public/index.js`, `public/index.html`
- Endpoints and providers: `public/api_endpoints.js`, `public/fairchem_provider.js`
- Domain state & events: `public/domain/moleculeState.js`, `public/domain/eventBus.js`
- Bonds and recompute: `public/domain/bondService.js`, `public/bond_render.js`
- Picking & manipulation: `public/core/pickingService.js`, `public/domain/manipulationService.js`, `public/domain/selectionService.js`
- Rendering: `public/render/scene.js`, `public/render/moleculeView.js`
- UI desktop HUD: `public/ui/desktopPanel.js` (+ toggles/sliders in `public/ui/`)
- Loaders & utils: `public/util/*`
- Physics helpers: `public/physics/sim-model.js` (structures), placeholders in `public/physics/*`
- VR/XR (optional): `public/vr/setup.js`


## License

See repository root for license information.
