# Running the Test Suites

```
npm test               # Jest unit/browser suites under tests/ and tests-browser/
npm run test:integration  # Integration target configured via jest.integration.config.js
npm run test:py        # Python parity checks under fairchem_local_server2/tests_py
npm run test:e2e       # Playwright WebSocket/DOM end-to-end coverage (serial execution)
```

# Diagnostic Flags
- `window.__MLIP_DEBUG_STRETCH = true` (or append `?bondStretchDebug=1` to the viewer URL) prints bond-opacity summaries from both the bond service and the dual-mesh renderer. Use this when validating stretched-bond translucency in live mode; disable it once finished to keep the console quiet.

# Ported Test Catalogue

This catalogue tracks only the suites that are fully ported to the protobuf/WebSocket stack. Entries are grouped hierarchically so related coverage is easy to locate. For each test, the intent captures what the suite safeguards, and the implementation notes summarise how the current code exercises that behaviour.

## Playwright WebSocket Suites (`tests-e2e/ws-*.spec.js`)

- **Backpressure & Counter Stability**
  - `ws-backpressure-ack.spec.js`: Adds a runtime shim that drops binary ACK frames, waits for the viewer to expose `WAITING_FOR_ACK`, then releases a single ACK and asserts the pipeline resumes with monotonically increasing seq counters.
  - `ws-sim-counters.spec.js`: Runs a continuous MD session, scrapes console telemetry, and confirms `simStep` and ACK identifiers strictly increase while the socket stays open.
  - `ws-stream-single-socket.spec.js`: Chains multiple MD/relax runs without reloading, observing the shared `WebSocket.readyState` and seq progression to ensure connection reuse without churn.
  - `ws-safe-sphere.spec.js`: Toggles the safe-sphere overlay twice during a stream and inspects outbound ClientActions to verify counters stay stable and no duplicate toggles fire.

- **Idle Initialisation & Baselines**
  - `ws-initial-idle-energy.spec.js`: Blocks auto-MD, captures the first idle `forces+energy` frame, and asserts the energy plot baseline is populated exactly once.
  - `ws-idle-positions-omitted.spec.js`: Requests an idle compute where the backend omits `positions`, then verifies atom transforms remain steady while energy ticks still appear.
  - `ws-mock-load.spec.js`: Loads the mock dataset through the system panel and inspects the next USER_INTERACTION payload to confirm `natoms`, `atomic_numbers`, and `cell` match the mock definition.

- **Cell Stress & Parameter Updates**
  - `ws-cell-stress.spec.js`: Starts a relax run against UMA, collects streamed frames, and ensures stress tensors feed through to the HUD alongside force-arrow updates.
  - `ws-live-md-param-updates.spec.js`: Adjusts MD temperature/friction mid-run, monitors outbound ClientActions, and confirms the following server frames echo the new parameters.

- **Counters, Interaction Echo, and UI Control**
  - `ws-counters-echo.spec.js`: Drags an atom, pauses, and begins MD while logging the HUD counters; validates server echoes track the viewer’s `userInteractionVersion` and total-interaction counts.
  - `ws-start-stop-idle-relax.spec.js`: Drives idle, MD, then relax via UI buttons, verifying each phase emits the correct verb and the button state/-LEDs sync with the run status.
  - `ws-ui-buttons-start-stop-idle-relax.spec.js`: Smoke-tests toolbar controls for stop/start across idle, MD, and relax, ensuring the underlying viewer flags (`isMdRunning`, `isRelaxRunning`) toggle in step.
- **Timeline Replay & Read-only Inspection**
  - `ws-timeline-controls.spec.js`: Validates play/pause/live button policy, ensuring mode transitions gate the correct controls.
  - `ws-timeline-interaction-lock.spec.js`: Enters timeline mode, asserts manipulations are rejected while paused, then resumes live streaming.
  - `ws-timeline-visibility.spec.js`: Checks hover reveal of the dock, overlay styling, and the live-resume path.
  - `ws-timeline-replay.spec.js`: Scrubs to historical frames, verifies signatures remain stable, and confirms playback returns to live streaming.
  - `ws-timeline-camera.spec.js`: Scrubs history, dispatches pointer and wheel gestures, and confirms camera rotation/zoom stay responsive while geometry edits remain blocked.
  - `ws-timeline-energy-marker.spec.js`: Enters timeline playback and asserts the energy plot displays (and clears) the playback marker.
  - `ws-timeline-slider-select.spec.js`: Clicks a single slider position and waits for timeline mode to activate on the requested offset (guards against the historical double-click requirement).
  - `ws-session-save-load.spec.js`: Captures a JSON snapshot, mutates geometry, loads the snapshot, verifies reset-to-last-load rehydrates the state/counters, and asserts MD runs resume after exiting timeline playback.
  - `ws-timeline-mesh-mode.spec.js`: Drives a historical selection, applies an opacity mask, and asserts atom/bond instances migrate to soft masters during the fade then return to their baseline solid/soft split when the mask clears.
- `ws-session-playback-resume.spec.js`: Runs MD on ROY, captures a JSON snapshot, reloads it, scrubs five frames back, then plays forward until live mode resumes—asserting at least 18 and at most 170 new MD frames arrive before the viewer reattaches to the live stream (protects against both stalls and runaway playback).

- **Ghost & Transparency Modes**
  - `ws-ghost-periodic.spec.js`: Loads a periodic XYZ, toggles ghost cells, and confirms ghost atoms/bonds live on the soft masters during live streaming and timeline playback (no regression to opaque pools).

## Session Snapshot Fixtures (manual QA)

- `fixtures/sessionSnapshots/idle-water.json`: Idle water geometry with a single idle energy sample to validate baseline restores without timeline playback.
- `fixtures/sessionSnapshots/md-benzene.json`: Benzene captured mid-MD run with resume metadata and timeline frames; load via the System panel session controls to benchmark file size and automatic MD resume behaviour.

## Jest Suites (`tests/x-*.spec.js`)

- **Cell & Periodic Boundary Controls**
  - `x-acetic-acid-ghosts.spec.js`: Loads a periodic XYZ, enables cell + ghosts, triggers `requestSimpleCalculateNow`, and confirms USER_INTERACTION carries a 3×3 cell plus ghost instance counts rise.
  - `x-api-cell-payload.spec.js`: Toggles the periodic UI switch with a WS hook installed and verifies the captured payload contains cell vectors, origin, and correct atom counts.
  - `x-cell-optimizer.spec.js`: Documents (via `test.failing`) the desired convergence behaviour for the variable-cell optimiser once reintroduced.
  - `x-cell-toggle.dom.spec.js`: Builds the HUD, injects a status node, and clicks the cell toggle twice to assert `showCell`, `showGhostCells`, and the button text rotate through ON/OFF labels.
  - `x-pbc-ghost-master-disabled.dom.spec.js`: Ensures ghost master meshes stay disabled while no thin instances exist, only enabling once ghosts are generated.
  - `x-pbc-long-bond-crossing.dom.spec.js`: Rebuilds bonds for a sample long-bond structure and checks thin-instance counts match expectations without double counting.
  - `x-pbc-long-bond-crossing-x.dom.spec.js`: Variant covering X-direction wrapping, validating midpoint diagnostics line up with wrapped coordinates.
  - `x-periodic-monoclinic-ui.dom.spec.js`: Renders the periodic section, seeds a monoclinic cell into `viewer.state`, and confirms lattice lengths/angles display correctly.
  - `x-periodic-omol25.dom.spec.js`: Ensures all OMOL25 elements render with their highlight class inside the mini periodic table while non-eligible elements remain unstyled.
  - `x-ui-pbc-bond-break-no-artifact.dom.spec.js`: Breaks a bond under PBC and asserts the scene graph contains no stray meshes at the origin.

- **Molecule Loading & Rendering**
  - `x-molecules.spec.js`: Iterates the molecule registry used on the landing page and asserts each entry exposes the required metadata fields.
  - `x-moleculeState.spec.js`: Validates `createMoleculeState` for non-periodic inputs by inspecting positions, bonds, and default selection state.
  - `x-moleculeState.pbc.spec.js`: Same factory but with a 3×3 cell, confirming periodic flags, lengths, and angles are normalised.
  - `x-molecule-reload.spec.js`: Loads a molecule twice through the viewer API and ensures geometry versions increment without duplicating atoms or ghosts.
  - `x-molecule-switch-highlight.spec.js`: Switches molecules after selecting an atom and verifies highlights and selection reset.
  - `x-no-stray-origin-primitives.spec.js`: Runs through common toggles (ghosts, forces) and asserts no mesh remains at `(0,0,0)` unless representing real geometry.
  - `x-visual-integration.spec.js`: Snapshot-style check that atom and bond thin-instance buffers match expected counts/material IDs for benzene.
  - `x-lighting-consistency.spec.js`: Inspects material colour coefficients on masters to ensure they match the tuned palette.
  - `x-red-sphere-switch.spec.js`: Toggles highlight modes and confirms the atom highlight mesh remains cyan (no regression to red).
  - `moleculeView.meshModes.spec.js`: Builds thin-instance groups with the dual-master renderer and verifies baseline assignments plus migrations keep buffers and counts in sync.

- **Parsers, Importers, and Prefills**
  - `x-base64Utf8.spec.js`: Exercises Base64⇄UTF-8 helpers and ensures XYZ comments yield temperature/cell metadata.
  - `x-xyz-loader.spec.js`: Parses XYZ samples, asserting `atoms`, `temperature`, and `cell` fields populate correctly.
  - `x-xyz-extended.spec.js`: Covers extended XYZ (`Properties=`) parsing, verifying property arrays align per atom.
  - `x-smilesLoader.unit.spec.js`: Normalises valid SMILES strings and rejects unsupported syntax.
  - `x-smiles-and-upload.dom.spec.js`: Simulates SMILES submission and XYZ upload in the System panel, confirming parsed output lands in `viewer.state.pendingLoad`.
  - `x-methyl-temperature-prefill.dom.spec.js`: Loads the methyl example and checks the MD slider/target temperature default to 500 K.
  - `x-temperature-slider-prefill.dom.spec.js`: Loads an XYZ with `temp=300K` and asserts UI/global targets reflect 300 K.
  - `x-temperature-tick-position.dom.spec.js`: Renders the temperature gauge and confirms tick labels/positions match design.

- **Highlight & Ghost Rendering**
  - `x-highlight.spec.js`: Selects an atom and validates the highlight mesh translation mirrors the atom position.
  - `x-highlight-visibility.spec.js`: Ensures the highlight mesh starts hidden and becomes visible only when selection exists.
  - `x-highlight-stray-cylinder.spec.js`: Selects and clears a bond, verifying the bond highlight mesh hides and buffers zero out.
  - `moleculeView.ghostBonds.spec.js`: Calls `createMoleculeView`, runs the bond service, and asserts `rebuildBonds()` now refreshes ghost thin instances automatically—no manual `rebuildGhosts()` calls required.
  - `x-ghost-bonds.spec.js`: Enables ghosts, rebuilds, and asserts ghost bond matrices match expected neighbour counts.
  - `x-ghost-crossing-bonds.spec.js`: Similar coverage focused on cross-cell bonds to prevent duplicate instances.
  - `x-ghost-selection.spec.js`: Picks a ghost and confirms the selection resolves to the primary atom with highlight repositioned.

- **Selection & Interaction Mechanics**
  - `x-selection-clear.spec.js`: Clicks empty space after selecting an atom to ensure the selection resets to `kind=null` and highlights hide.
  - `x-selection-null-guard.spec.js`: Manually nulls the selection state, calls selection helpers, and verifies no crashes while state normalises.
  - `x-selection-panel.dom.spec.js`: Renders the selection panel, simulates atom/bond selections, and checks DOM text and action buttons update.
  - `x-selection-periodic-click.dom.spec.js`: Fires synthetic clicks on periodic table entries and asserts the selection service queues the element and updates UI.
  - `x-selection-persistence.spec.js`: Mutates positions (should keep selection) then alters topology (should clear selection) to confirm the guard logic.
  - `x-selection-service.spec.js`: Exercises `clickAtom`/`clickBond` ordering, ensuring mutual exclusion and rotation-group metadata propagation.
  - `x-pointer-observable-bug.spec.js`: Re-initialises the viewer twice against a stubbed scene and asserts only one Babylon pointer observer attaches.
  - `x-joystick-delta.spec.js`: Runs joystick input samples through the bond rotation helper to confirm deltas throttle and scale correctly.

- **Bond Rotation & Picking**
  - `x-bondOpacityIsolation.spec.js`: Rotates a bond and confirms opacity isolation only affects the active pair.
  - `x-bondRotateEnergyCount.spec.js`: Steps through discrete rotations and checks the energy plot registers the correct number of ticks.
  - `x-bondRotationSoftBond.spec.js`: Marks a bond as soft, attempts rotation, and asserts atom coordinates remain unchanged.
  - `x-bondRotationUI.spec.js`: Drives the rotation UI controls and ensures both the viewer manipulation calls and DOM labels update.
  - `x-bondService.spec.js`: Mutates the cell and runs bond service refresh to verify bond lists/caches regenerate.
  - `x-picking.spec.js`: Constructs synthetic pick hits and ensures atom/bond resolvers return accurate metadata.
  - `x-picking-integration.spec.js`: Dispatches pointer events through the picking service, validating selection changes, drag state, and camera detaching.
  - `x-rotate-recompute.spec.js`: Executes bond rotation commands and checks recompute helpers bump geometry versions and atom indices as expected.
  - `x-rotation-inversion.spec.js`: Applies opposing rotations and verifies the resulting quaternion correctly inverts to prevent direction flips.

- **Simulation Controls & Auto-MD**
  - `x-auto-md-panel-toggle-reset.spec.js`: Toggles auto-MD via the panel, injects a completion frame, and asserts the toggle resets with updated status text.
  - `x-auto-md-run-reset.spec.js`: Starts auto-MD, triggers a viewer reset, and verifies start/stop calls plus timer cleanup.
  - `x-auto-md-starts-roy.spec.js`: Loads the ROY example with auto-MD enabled and confirms the first outbound ClientAction is `START_SIMULATION` (MD).
  - `x-md-temperature-slider.spec.js`: Adjusts the temperature slider, runs `viewer.mdStep`, and checks the stubbed `requestSingleStep` receives the updated temperature.
  - `x-md-temperature-sync.dom.spec.js`: Verifies slider labels start at 1500 K, update after user input, and propagate to MD step parameters.
  - `x-md-initial-temperature.ws.spec.js`: Loads an XYZ containing temperature metadata and asserts the initial USER_INTERACTION + MD params respect the embedded value.
  - `x-md-rps-label.ws.spec.js`: Starts/stops MD twice and checks the HUD RPS label reflects the latest frame rate.
  - `x-methyl-temperature-prefill.dom.spec.js`: (Shared with importer coverage) specifically ensures the example pins the 500 K target.

- **MD & Relaxation Streaming**
  - `x-md-drag-exclusion-vr.spec.js`: Simulates VR grip selection and injects MD frames to confirm held atoms stay fixed while others move.
  - `x-md-drag-position-lock.dom.spec.js`: Starts MD, performs a desktop drag, feeds frames, and ensures the dragged atom ignores server updates until release.
  - `x-md-forces-update.spec.js`: Injects sequential MD frames with varying forces and verifies shaft/head matrices change accordingly.
  - `x-relax-force-cache.spec.js`: Runs relax single-step, inspects `viewer.state.forceCache`, and confirms idle recompute bumps the version only after geometry change.
  - `x-relax-forces-update.spec.js`: Feeds relax frames with different forces and checks force arrow thin-instance buffers refresh.
  - `x-relax-run-post-init-enable.spec.js`: Toggles the relax run feature flag, reloads, and confirms persistence across viewer resets.
  - `x-relax-single-step-network.spec.js`: Calls the relax helper, ensuring exactly one `START_SIMULATION` fires and the promise resolves on the first injected frame.

- **Force & Energy Visualisation**
  - `x-force-arrowheads.spec.js`: Uses the Babylon stub, populates forces, emits `forcesChanged`, and checks thin-instance matrices for shafts/heads.
  - `x-force-render.spec.js`: Confirms the colour buffer retains the calibrated red emissive values (`0.95/0.05/0.05/1`).
  - `x-force-provider.spec.js`: Mocks a provider, performs recomputes, and verifies the cache version changes only when geometry does.
  - `x-force-update-perturbation.spec.js`: Perturbs atom positions, triggers force recompute, and asserts the magnitude histogram recalculates.
  - `x-forces-toggle.dom.spec.js`: Toggles force visibility via the HUD, checking mesh visibility and preference persistence.
  - `x-energy-reset.dom.spec.js`: Seeds the energy plot, presses reset, and confirms datasets clear and controls hide until new data arrives.

- **UI Panels, Toggles, and Resets**
  - `x-desktop-panel.dom.spec.js`: Builds the desktop panel and ensures expected controls (MD slider, friction slider, molecule select, cell toggle) exist.
  - `x-desktop-bond-pick.dom.spec.js`: Stubs picking, dispatches a pointer event, asserts bond selection updates, and exercises rotation buttons.
  - `x-desktop-rotation-input.dom.spec.js`: Modifies numeric rotation inputs, switches tabs, and confirms values persist.
  - `x-ui-toggles.dom.spec.js`: Iterates the panel toggles, validating `aria-checked` and viewer flags align with button text.
  - `x-run-completion-toggles.dom.spec.js`: Starts MD/relax runs, injects completion frames, and checks toggles revert to OFF with updated labels.
  - `x-reset-button.dom.spec.js`: Clicks the reset button, ensures it disables while awaiting the `resetToInitialPositions` promise, then re-enables.
  - `x-reset-invalidation.dom.spec.js`: Seeds caches, triggers reset, asserts caches/HUD overlays return to defaults, and verifies baseline geometry is restored after add/remove edits.

- **User Interaction & Counter Handling**
  - `x-user-interaction-during-md.dom.spec.js`: Starts MD, injects stale USER_INTERACTION frames (ignored) followed by in-order frames that are applied.
  - `x-user-interaction-when-idle.dom.spec.js`: Injects idle frames through the WS stub and confirms geometry updates immediately.
  - `x-user-interaction-version-isolation.spec.js`: Bumps local counters, feeds lower-count frames (ignored), then matches counts to verify acceptance.
  - `x-interaction-invalidation.spec.js`: Triggers interaction invalidation and assures stale frames with mismatched counters are dropped.
  - `x-ws-double-connect.dom.spec.js`: Inits the viewer twice with a stubbed WS client and asserts only one underlying connection is established.

- **Dragging & Atom Manipulation**
  - `x-atom-drag-camera.spec.js`: Simulates drag start, checks manipulation state, and ensures camera inertial offsets zero out while dragging.
  - `x-atom-drag-screen-plane.spec.js`: Runs the drag solver with sample pointer data and verifies resulting positions stay on the computed plane.
  - `x-drag-recompute.spec.js`: Performs drag start/end, confirms bond recompute occurs once, and geometry version increments.
  - `x-cameraSuppression.spec.js`: Ensures the ArcRotate camera remains frozen during active drags by checking inertial offsets and position locks frame-by-frame.

- **VR / XR & Related Assets**
  - `x-http-vr-asset.spec.js`: Mocks `fetch` to fail and ensures the VR asset loader handles the error path cleanly.
  - `x-vr-setup.spec.js`: Invokes the VR setup helper with Babylon stubs, checking controller/helper nodes register correctly.
  - `x-vr-module-availability.spec.js`: Stubs `navigator.xr` absence and asserts capability flags fall back gracefully.
  - `x-vr-controller-handshake.spec.js`: Simulates controller attach/detach events and verifies viewer latches update.
  - `x-vr-atom-drag-initiation.spec.js`: Emits mock controller picks to ensure drag begins and MD pauses.
  - `x-vr-drag-stability.spec.js`: Feeds a pose sequence into the drag solver and checks position deltas stay within tolerance.
  - `x-vr-dynamic-import.spec.js`: Exercises the dynamic import path for VR modules, ensuring the promise resolves and caches results.
  - `x-vr-highlight-alignment.spec.js`: Updates head pose while a highlight is active and confirms the highlight mesh repositions correctly.
  - `x-vr-hud-energy-position.spec.js`: Moves the camera and validates the HUD energy panel follows with the expected offset.
  - `x-vr-laser-thickness.spec.js`: Adjusts controller laser configuration and asserts mesh widths respect clamp values.
  - `x-vr-late-bond-rotation.spec.js`: Queues rotation callbacks post-latch and ensures the rotation group updates.
  - `x-vr-post-release-latch.spec.js`: Releases a VR drag, advances timers, and verifies atoms remain latched until timeout expires.
  - `x-vr-reuse-parity.spec.js`: Re-enters VR mode twice, ensuring controllers/meshes are reused without duplication.
  - `x-vr-trigger.spec.js`: Fires synthetic trigger events and checks viewer callbacks execute with normalised payloads.
  - `x-xr-controls-core.spec.js`: Validates the XR controls model exposes default state and emits observer callbacks on mutation.
  - `x-xr-dropdown.dom.spec.js`: Renders the XR dropdown, selects each option, and ensures `viewer.xr.setMode` receives the right value.
  - `x-xr-energy-hud.dom.spec.js`: Builds the XR energy HUD and confirms binding to the energy data source updates labels.
  - `x-xr-hud.dom.spec.js`: Constructs the XR HUD scaffold, verifying button groups and event wiring.
  - `x-xr-hud-buttons.spec.js`: Clicks each HUD button and checks the associated viewer handlers fire.
  - `x-xr-hud-energy-depth.spec.js`: Adjusts the depth slider and asserts stored settings plus style transforms reflect the new depth.
  - `x-xr-hud-energy-plot-picking.spec.js`: Synthesises pointer events over the HUD chart and ensures handlers return the correct step index.
  - `x-xr-hud-energy-scale.spec.js`: Moves the energy scale slider and confirms the chart rerenders with the new amplitude.
  - `x-xr-hud-position.spec.js`: Applies custom offsets and checks HUD anchor matrices translate accordingly.
  - `x-xr-hud-texture-clamp.spec.js`: Inspects HUD materials to confirm texture clamp modes use `CLAMP_ADDRESSMODE`.

- **Miscellaneous Utilities**
  - `x-ar-auto-normalize.spec.js`: Runs augmented-reality auto-normalisation on various vectors and asserts values clamp within configured bounds while preserving direction.
  - `x-spherical-radial-modes.spec.js`: Exercises the spherical drag helper to verify radial displacement magnitudes and directions.
  - `frameBuffer.serialize.spec.js`: Ensures timeline frame buffers export/import without losing metadata or coordinate fidelity.
  - `sessionStateManager.spec.js`: Stubs the viewer/WS stack to verify snapshot capture and restoration seed counters, timeline state, and full-update syncs.
  - `controlMessageEngine.spec.js`: Asserts control-message ranges resolve by priority, preserve labels/notes, and return the correct per-frame speed/callout/opacity actions.
  - `timelinePlaybackController.spec.js`: Exercises the standalone playback scheduler, validating fps cadence, overrides, snapshot helpers, and the new setter utilities for auto-play/loop configuration.

> _Cross-reference_: Tests listed only once above cover all ported `x-*` suites. If you add a new suite, update this hierarchy so related coverage stays discoverable.
