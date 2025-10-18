## Design and Architecture

This document describes the front‑end architecture, desktop MD/Relax workflows, and a per‑file map of the `public/` folder. It complements the top‑level README and focuses on how user interactions propagate through state, rendering, and the backend API.


## Architecture overview

The viewer follows a modular pattern:
- Domain (state): authoritative molecule data and toggles, event bus for change propagation
- Render: Babylon.js scene setup and thin‑instance meshes for atoms/bonds/forces/ghosts; selection highlight meshes
- Core interactions: picking and drag orchestration; bond rotation; selection model
- UI: desktop panel with MD/Relax controls, PBC cell controls, toggles and sliders; touch shims
- Utilities: molecule loading/parsing, periodic boundary math, constants and element info
- Backend wiring: endpoint selection, optional local force providers
- Optional VR/XR: session setup, controller lasers, drag, bond rotation via joystick, two‑hand scaling

Dataflow (desktop):
1) Input (mouse/touch/keys) → Picking & manipulation → Domain state updates (positions, selection, toggles)
2) Domain emits events via a simple event bus
3) Render layer subscribes and updates thin instances and highlight meshes
4) UI listens to selection/positions/metrics to update labels, plots, and button states
5) API loops (MD/Relax) form requests from current domain state, parse results, and re‑emit updates

Key invariants:
- The domain state is the single source of truth for positions, bonds, selection, and simulation metadata.
- All render updates happen via events; there is no direct mutation of Babylon instances from UI.
- In-flight API results are version‑checked to avoid applying stale frames after user edits.


## Desktop MD/Relax flow (details)

Where: `public/index.js` orchestrates step/run functions and installs handlers during viewer initialization.

Request formation:
- Gather elements (Z or symbols), positions, and optional periodic cell vectors from `moleculeState`.
- MD: steps=1 (per iteration when running), temperature from the UI slider; Relax: steps=N (configurable batch size) or steps=1 in a loop.
- If a force/energy cache matches the current geometry, attach a `precomputed` block to skip redundant backend work.
- POST to `public/api_endpoints.js` paths (`/serve/simple`, `/serve/relax`, `/serve/md`).

Response parsing and application:
- On success, update forces, energy, positions (and velocities for MD). Use `markPositionsChanged()` which wraps PBC when enabled and emits `positionsChanged`.
- Emit `forcesChanged` and `bondsChanged` as needed; recompute bonds after significant geometry changes.
- Update energy plot and instantaneous temperature label (UI hooks).

Loop pacing and robustness:
- Minimum inter‑request interval enforced (default ~30 ms). Continuous runs schedule the next iteration after parsing and render event emission.
- Errors trigger exponential backoff (start ~200 ms up to ~5 s) before retrying; after a failure streak, the loop aborts with a user‑visible message.
- Edits during flight update an interaction version; stale responses are partially applied or dropped.


## Rendering strategy

Thin instances for performance:
- Atoms: instanced spheres with per‑instance color/scale derived from element info
- Bonds: instanced cylinders with per‑instance alpha to shape appearance for crossings/opacity effects
- Force vectors: optional instanced arrows aligned with per‑atom forces
- Ghost atoms/bonds and cell lines for PBC visualization

Highlights:
- Dedicated highlight meshes for atom and bond selection are created once and re‑positioned; they are parented to molecule masters so they inherit global transforms (including VR root rotation/scale).

Lighting:
- Scene sets up hemispheric fill and a camera‑aligned directional “headlight”. A helper (`public/render/lighting.js`) can attach consistent lighting synced to camera on every frame.


## Per‑file map (public/)

Entry and orchestration
- `index.html` – Boots the viewer, wires panel buttons, and installs latency wrappers for debugging
- `index.js` – The main entry: creates the scene/view/state, sets endpoints, runs relax/MD loops, maintains force cache, and exports the runtime API
- `api_endpoints.js` – Centralized endpoint map for `/serve/simple`, `/serve/relax`, `/serve/md`, `/serve/health`
- `fairchem_provider.js` – Simple HTTP provider to call a backend FairChem force/energy service (with retry)

Domain (state and services)
- `domain/moleculeState.js` – Core state (elements, positions, bonds, selection, toggles, dynamics, cell) and event bus bridges; wrap on position updates
- `domain/bondService.js` – Recomputes bonds (min‑image and dedup), sets per‑edge opacity and emits `bondsChanged`
- `domain/selectionService.js` – Click handlers and selection state transitions
- `domain/manipulationService.js` – Atom drag intersectors and bond rotation (across strong bonds only); emits `positionsChanged`
- `domain/eventBus.js` – Lightweight emitter (on/once/emit/clear)

Rendering
- `render/scene.js` – Babylon engine, scene, camera, pointer observables and lighting setup
- `render/moleculeView.js` – Creates thin‑instance masters for atoms/bonds/forces/ghosts, selection highlight meshes, and resolves picks; subscribes to state events
- `render/lighting.js` – Optional camera‑synced ambient + directional lighting helper

Core interactions
- `core/pickingService.js` – Atom/bond picking, selection, drag orchestration (camera detach/freeze while dragging); integrates with manipulation service and energy/plot hooks
- `selection-model.js` – Selection model helpers and bond orientation coloring

UI (desktop)
- `ui/desktopPanel.js` – Left panel HUD: live metrics, selection summary, MD/Relax buttons, force/cell toggles, temperature/friction sliders, molecule selector, SMILES, XYZ upload, reset
- `ui/forcesToggle.js`, `ui/cellToggle.js`, `ui/temperatureSlider.js`, `ui/frictionSlider.js` – Widget helpers
- `ui/touchControls.js` – Touch gesture mapping to camera and picking
- `ui/moleculeSelect.js` – Example selector and URL helpers
- `ui/errorBanner.js` – Transient error banner

Utilities
- `util/moleculeLoader.js` – Parses URL params, loads molecules (examples/SMILES/XYZ), enables PBC, and triggers bond recompute
- `util/xyzLoader.js`, `util/smilesLoader.js` – Parsers and remote fetch for SMILES via PubChem
- `util/pbc.js` – Reciprocal lattice wrapping, orthorhombic cell estimation, and cell helpers
- `util/constants.js` – Pacing and MD constants
- `elements.js`, `data/periodicTable.js` – Element info (colors, radii, scales)
- `bond_render.js` – Bond inference and opacity shaping (non‑PBC version)

Physics (helpers & placeholders)
- `physics/sim-model.js` – Simulation arrays, masses, stress/kinetic helpers
- `physics/forcefield.js` – Local toy forcefield (LJ + harmonic bonds) used by optional local provider
- `physics/force-provider.js` – Local provider façade wrapping `forcefield.js`
- `physics/optimizer.js` – Standalone optimizer (SD/CG) for local experiments
- `physics/integrators.js`, `physics/relaxation-runner.js` – Placeholders (reserved for future)

VR/XR (optional)
- `vr/setup.js` – Session bootstrap, controller lasers, pick/drag, joystick bond rotation, two‑hand scaling, and HUD bars


## User interactions (desktop) and code paths

Atom select: click on sphere → `core/pickingService.js` → `selectionService.clickAtom` → event → highlight in `render/moleculeView.js`

Bond select: click on cylinder → `core/pickingService.js` → `selectionService.clickBond` → orientation metadata in `selection-model.js` → highlight update

Atom drag: mousedown → detach camera → intersector in `manipulationService.beginDrag` → live `positionsChanged` → bonds optionally recomputed → mouseup reattach camera

Bond rotate: UI or keyboard gesture → `manipulationService.rotateBond(deltaAngle)` → updates a branch of atoms around bond (strong bonds only)

Force vectors: toggle in HUD → `moleculeState.toggleForceVectorsVisibility` → view rebuilds instanced arrows

PBC controls: toggle and cell edits in HUD → `moleculeState.toggleCellVisibilityEnhanced` / `util/pbc.js` updates → wrap on position changes; ghost atoms/bonds updated in view

MD temperature/friction: sliders update values exposed on `moleculeState.dynamics` and sampled when building /serve/md requests


## Eventing and versioning

Change events:
- `positionsChanged`, `bondsChanged`, `selectionChanged`, `forcesChanged` emitted via `moleculeState.bus`

Version guards in `index.js`:
- `structureVersion` and `interactionVersion` are incremented on edits; outgoing requests carry these values and results are applied only if versions still match. This prevents stale or conflicting updates (e.g., when dragging during a long request).

Force/energy cache:
- Latest force field result for the current geometry is cached and attached as `precomputed` to the next `/serve/relax` or `/serve/md` call. Cache is invalidated on edits.


## Extensibility notes

- New calculators: extend the Python backend and keep the viewer unchanged; it posts to `/serve/*` endpoints using relative URLs by default.
- Custom renderables: add new thin‑instance masters in `moleculeView.js` and subscribe to domain events.
- New UI toggles/sliders: add a small helper in `ui/` and wire it in `ui/desktopPanel.js`.
- VR/XR: `vr/setup.js` cleanly composes with the same state/view so highlights and transforms remain consistent.


## Known placeholders and future work

- `public/physics/integrators.js`, `public/physics/relaxation-runner.js`, `public/interaction/dragCore.js`, `public/interaction/pickingFacade.js`, `public/ui/rotationControllerCore.js`, `public/ui/energy-plot.js` are scaffolds/placeholders. They can be filled or removed as the implementation stabilizes.
- Consider unifying element data across UI and render via `data/periodicTable.js` exclusively.
- Optional: migrate lighting setup to `render/lighting.js` everywhere and tune intensities per scene.


## Appendix: Desktop flow checklist (engineer’s view)

Inputs/outputs
- Input: elements (Z/symbols), positions (Å), optional cell (Å vectors), temperature/friction (MD), step counts
- Output: updated positions, forces, energies, stress (optional), velocities/temperature (MD)

Edge cases
- Non‑finite backend responses: MD/Relax loop aborts and shows an error banner
- Edits during in‑flight response: partial-apply guards prevent teleports; force cache invalidated
- Very large or tiny molecules: camera auto‑framing and selection highlight scales ensure visibility

Success criteria
- Continuous MD/Relax maintains near‑interactive pacing with stable visuals and predictable selection behavior
- User edits reflect immediately; subsequent steps use the new geometry

