# MLIPView

MLIPView is UMA Fairchem’s reference workstation for interactively steering machine-learned interatomic potential (MLIP) simulations. A Babylon.js viewer renders atoms, bonds, energies, and VR overlays in real time while a Ray Serve–backed FastAPI service streams UMA predictions, handles user edits, and records a 500-frame timeline for rewind and playback. The project targets desktop, VR, and AR displays so computational chemists can inspect trajectories, adjust structures, and validate results from a single tool.

Key capabilities:
- Real-time idle, MD, and relax streaming with per-atom gating during drags and bond rotations.
- Timeline mode with scrubbable history, retained energy markers, and read-only camera controls for investigating prior frames; playback now resumes live streaming with bounded frame replay so simulations rejoin the backend seamlessly.
- Session snapshots with JSON export/import so timelines, energy traces, and counters can be restored exactly across runs (backed by the unified `SessionStateManager` interface).
- VR/AR support tailored for Meta Quest, including controller pick/rotate and HUD overlays.
- Rich debugging surface (`viewerApi`) plus exhaustive Jest, Playwright, and pytest coverage to protect behaviour.

## Repository Guides

- `backend_design.md` – current backend architecture, protobuf schema, and session lifecycle.
- `frontend_design.md` – viewer orchestration, WebSocket client behaviour, and interaction model.
- `testing.md` – catalogue of Jest, Playwright, and pytest coverage.
- `test_hooks.md` – convenient hooks and globals for inspecting WebSocket traffic during tests.
- `QUEST_DEBUG_CONSOLE.md` – capturing Meta Quest headset logs for diagnosing VR issues.
- `VR_SETUP_README.md` – configuring VR runtimes and controller mappings for local testing.

## Prerequisites

- Python 3.12 (CUDA 12.4 capable GPU recommended for UMA; set `MLIPVIEW_FORCE_CPU=1` to stay on CPU).
- Node.js ≥ 20.19.0 and npm.
- `protoc` 30.x (used by `npm run gen:proto:js`).

## Installation

```bash
git clone <repo-url> mlipview
cd mlipview

# Python environment
python -m venv mlipview_venv
source mlipview_venv/bin/activate            # Windows: mlipview_venv\Scripts\activate
pip install -U pip
pip install -e .                             # installs backend + test dependencies from pyproject.toml

# Node / frontend toolchain
npm install
```

Regenerate protobuf stubs any time `fairchem_local_server2/session.proto` changes:

```bash
npm run gen:proto:js
```

## Running the Backend

The backend uses Ray Serve to host both the WebSocket ingress and UMA predictor.

```bash
source mlipview_venv/bin/activate
python -m fairchem_local_server2.serve_ws_app \
  --ngpus 1          # use 0 to stay on CPU
```

Key environment variables:

- `MLIPVIEW_FORCE_CPU=1` – force CPU execution even if CUDA is available.
- `UMA_DEVICE=cuda|cpu` – override predictor device selection.
- `MLIPVIEW_RESUME_DEBUG=1` – enable verbose console logging of ACK/seq bookkeeping during timeline resume debugging.

The server exposes WebSocket `/ws` plus JSON health endpoints under `http://127.0.0.1:8000/serve`.
The WebSocket is protobuf-only; the viewer auto-generates the required stubs via `npm run gen:proto:js`.

## Running the Frontend

```bash
npm run dev          # Vite dev server on http://localhost:5173 (proxying backend ws/http)
```

To point the viewer at a different backend origin, set `window.__MLIPVIEW_SERVER = 'http://hostname:port'` before calling `initNewViewer`.

Production build & preview:

```bash
npm run build
npm run preview      # serves dist/ at http://localhost:5174
```

## Testing

Backend (pytest):

```bash
./mlipview_venv/bin/python -m pytest -q
```

Frontend unit / integration (Jest):

```bash
npm test -- --runInBand
# or faster loop without spinning backend mocks:
npm run test:fast
```

Playwright E2E suites:

```bash
npm run test:e2e          # full battery
npm run test:e2e:smoke    # smoke subset
```

Combined convenience targets are available in `package.json` (`test:py`, `test:integration`).

## Deployment Notes

- The repository ships sample nginx config (`nginx_routing.conf`) for proxying `/serve` to the Python API and static assets to the Vite build.
- `server-app.js` can be used to serve the built frontend (`dist/`) behind Express if nginx is not available.
- For VR usage, follow the headset-specific setup in `VR_SETUP_README.md` and consult `QUEST_DEBUG_CONSOLE.md` for debugging tools.
- Timeline resume issues: run the viewer with `?debug=1` and set `MLIPVIEW_RESUME_DEBUG=1` on the backend to mirror ACK/seq flow; the Playwright regression `ws-session-playback-resume.spec.js` demonstrates a healthy JSON load → playback → live resume cycle.

## Troubleshooting & Support

- Force UMA to stay on CPU: `MLIPVIEW_FORCE_CPU=1` when launching the backend.
- Reconnect issues or `WAITING_FOR_ACK` notices: see the counter/backpressure tests listed in `testing.md`.
- Need to inspect outbound WebSocket payloads during development: wire up the helpers in `test_hooks.md` or enable `window.__MLIPVIEW_DEBUG_WS`.

For deeper architectural context or extension points, start with `backend_design.md` and `frontend_design.md`. Let the maintainers know if you spot documentation gaps—everything above is kept current with the main branch.***
