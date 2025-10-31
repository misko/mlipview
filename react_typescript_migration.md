# React + TypeScript Migration Plan

## 1. Goals & Success Criteria
- Introduce TypeScript for all browser code so molecule/WS/timeline logic gains static guarantees, while preserving current Vite build speed.
- Adopt React for UI orchestration (panels, overlays, timeline, VR HUD shells) without regressing Babylon scene performance.
- Keep backend APIs (protobuf WS, timeline JSON) unchanged during the migration; schema updates remain a separate concern.
- Maintain green Jest, Playwright, and pytest suites at the end of every phase so we can ship incremental improvements rather than a single “big bang”.

## 2. Current Frontend Baseline (2025-01)
- Vite (ESM) builds plain JS modules under `public/`; Babylon scene orchestration lives in `public/index.js`.
- UI is hand-wired DOM manipulation + custom event bus. State surfaces (`SessionStateManager`, `createMoleculeState`) are plain objects with implicit contracts.
- Tests lean on Jest + Playwright in JS. Proto stubs are generated via Buf (`@bufbuild/protobuf`) into JS.
- Global viewer API (`window.viewerApi`) exposes hooks that tests and VR tooling depend on.

These traits favour an incremental migration: Vite already understands TS/JS mixtures, and the Babylon layer can stay outside React until we are ready to wrap it.

## 3. Tooling Decisions (recommended)
- Use Vite’s built-in TypeScript pipeline (`tsconfig.json`, `vite.config.js` updates). Enable `allowJs` + `checkJs` initially so existing modules can live alongside TS.
- Adopt ESLint + `@typescript-eslint` with the new config once the TS foundation is in place.
- Keep Jest as the unit runner via `ts-jest` or Vite’s `vitest` (evaluate during Stage 1). Playwright already supports TS specs.
- Continue generating protobuf stubs with Buf, but switch output to TypeScript (`--es_opt=target=ts`) once core runtime compiles under TS.
- Source directories:
  - `src/react/` for new React code (hooks, providers, UI components).
  - `src/viewer/` (or reuse `public/`) for Babylon + domain services, migrated gradually to `.ts`.
- Preserve the existing `public/index.js` entry during the bridge period; introduce a new React entry (`src/main.tsx`) and run both until parity is proven.

## 4. Why TypeScript + React
**Pros**
- Static typing for molecule state, protobuf payloads, and WS counters reduces regressions when touching complex flows (drag gating, timeline resume).
- React brings declarative UI state, making it easier to reason about panel/timeline/VR overlays, and to integrate future desktop ↔ VR layouts.
- Ecosystem access: React Query, Zustand/Recoil, and Suspense patterns can replace ad-hoc singletons (`reconnectState`, viewer caches) over time.
- Tests gain richer tooling (React Testing Library, component-level snapshots).

**Cons / Costs**
- Initial productivity dip while retrofitting type information onto large JS modules.
- React render loops must not interfere with Babylon’s render cadence; extra care needed around refs and effect timing.
- Some Playwright specs poke DOM structure directly; migration will require adapter shims or updated selectors.
- Proto codegen + TypeScript may increase build time without careful caching.

## 5. Migration Stages

### Stage 0 – Preparation (1 sprint)
- Add `tsconfig.json` with `"allowJs": true`, `"checkJs": false`, `"jsx": "react-jsx"`.
- Update Vite config to accept `.ts`/`.tsx`.
- Configure `npm` scripts: `typecheck`, `lint`, `test:unit` (Jest+ts-jest or Vitest), `test:e2e`.
- Install React 18, ReactDOM, `@types/react`, `@types/react-dom`, `typescript`, `@typescript-eslint/*`.
- Add a CI lane for `npm run typecheck` (runs `tsc --noEmit`).
- Tests: run existing Jest/Playwright/pytest baseline to prove parity before touching runtime code.

### Stage 1 – Type Safety Foundation (2–3 sprints)
- Convert shared utilities and newly added modules to `.ts` with minimal changes (helpers in `public/util`, constants, periodic math).
- Introduce JSDoc typedefs for high-churn files (`public/index.js`, `moleculeView.js`) so TS consumers can import their types while they stay JS.
- Define central TypeScript interfaces for:
  - Molecule state (`Atom`, `Bond`, `MeshMode`).
  - WebSocket client contract (`FairchemWsClient`).
  - Timeline frame structures.
- Update `SessionStateManager` to `.ts` and export these types; fix schema migration bug (case `schemaVersion===6` never executes) while under type checking.
- Tests: enable `ts-jest` (or Vitest) for the converted modules. Add `tsc --noEmit` to CI gating.

### Stage 2 – React Shell & Providers (2 sprints)
- Create `src/react/App.tsx` that mounts:
  - Babylon canvas container.
  - Desktop panels (wrapped in React components).
  - Timeline dock.
- Build context providers for viewer state (`ViewerStateContext`), WebSocket status, timeline playback, and configuration.
- Wrap existing imperative modules by exposing React hooks (`useViewerApi`, `useTimeline`, `useWsClient`). Internally they delegate to the existing singletons until those modules become TypeScript.
- During the bridge period, render the existing DOM-driven UI via a compatibility layer so React and the legacy code can share the Babylon instance.
- Add integration tests with React Testing Library for new components (`TimelinePanel`, `ErrorBanner`).
- Smoke Playwright tests against the React entry flag (e.g., `?ui=react`) to ensure end-to-end parity.

### Stage 3 – Component Migration (4–5 sprints, iterative)
- Incrementally port UI modules to React components:
  1. Desktop control panel (temperature slider, MD/Relax toggles).
  2. Timeline dock (buttons, slider, playback indicators).
  3. VR HUD overlay shells (maintain Babylon GUI integration via refs).
  4. Error banner, reconnect prompts, library picker.
- Replace ad-hoc event wiring with React state/hooks. Use `useReducer` or Zustand for shared state currently coordinated by the viewer orchestrator and SessionStateManager.
- Convert high-value domain modules to `.ts`: `fairchem_ws_client`, `controlMessageEngine`, `timelinePlaybackController`, `frameBuffer`.
- For each module migration, update Jest specs to consume TypeScript types; convert specs themselves to `.ts` where it simplifies mocks.
- Maintain a feature toggle to switch between legacy and React UI so Playwright can run both during transition (`?legacyUI=1`).

### Stage 4 – Remove Legacy Entry (2 sprints)
- Once all UI panels and hooks run under React, retire `public/index.js` entry in favor of `src/main.tsx`.
- Move Babylon orchestrator into a dedicated hook (`useBabylonViewer`) that mounts once and communicates through contexts.
- Delete unused DOM helper modules, replace with React wrappers.
- Finalize TypeScript conversion of remaining files; disable `allowJs`.
- Update docs (`frontend_design.md`, testing guides) to describe the React architecture.

### Stage 5 – Hardening & Expansion (ongoing)
- Enable strict type checking (`"strict": true`), progressively fixing any remaining `any`.
- Add Storybook (optional) for VR/desktop panels to ease design iterations.
- Explore state-management refinements (e.g., React Query for WS data, Zustand for viewer state) once base migration is stable.
- Monitor bundle size and performance; leverage React lazy loading for heavy editor panels.

## 6. Testing & Release Strategy
- **Per stage** run full Jest + Playwright + pytest. Require at least `ws-*` E2E suites before merging UI migrations.
- Add smoke tests for both legacy and React entries until cutover (Playwright config can parametrize the UI flag).
- Hook `npm run typecheck` and `npm run lint` into pre-merge CI.
- For backend parity, keep protobuf schema untouched; run pytest after each proto codegen/tooling update.

## 7. Risks & Mitigations
- **React rendering interfering with Babylon**: keep Babylon scene management inside a stable ref; run `scene.render()` via the existing render loop, not React re-renders. Use `useLayoutEffect` carefully.
- **Typed protobuf generation churn**: cache generated files in repo; add tests to diff JS vs TS outputs to avoid accidental behaviour changes.
- **Playwright selector drift**: create data-testid attributes in React components mirroring current DOM structure so existing specs change minimally.
- **Team ramp-up**: schedule pair sessions on TS/React patterns; document conversions in `testing.md` or a new `frontend_migration.md`.

## 8. Recommended Milestones
1. **M0** – Tooling landed (`tsconfig`, React deps, CI typecheck). All tests green.
2. **M1** – Core state modules (`SessionStateManager`, WS client) compiling under TS with types exported.
3. **M2** – React shell optionally renderable (`npm run dev:react`), Playwright smoke passing.
4. **M3** – Feature flag default switched to React UI; legacy entry kept behind `?legacy=1`.
5. **M4** – Legacy removed, strict TS on, documentation updated.

Each milestone should ship independently (tag or release branch) with changelog notes and validated regression logs.

## 9. Open Questions for Planning Session
- Should we adopt Vitest in lieu of Jest when TypeScript is enabled, or keep Jest to avoid widespread config changes?
- Do we want to generate TypeScript types from protobuf on the backend as well (for shared contracts), or keep Python codegen separate?
- Is there appetite for introducing a lightweight state manager (Zustand/Redux Toolkit) during Stage 3, or should we defer until after cutover?

Collect answers during Step 5 review so implementation tickets can be sized accurately.

## 10. Current Progress & Decisions (2025-10-31)
- Step 1 review completed: documentation and implementation snapshots updated to reflect schema v6 and the single-session-manager architecture.
- Legacy snapshot migrations and the unused `stateStore` singleton were removed; `SessionStateManager` now hard-rejects non-v6 schemas.
- Targeted Jest suites plus the full `npm run test:fast` pass succeeded after the cleanup; backend parity (`npm run test:py`) also green (pydantic v2 migration deferred).
- Stage 2 in progress: added optional React shell (`?ui=react`) via `src/main.tsx`/`react/App.tsx`, which reuses `public/bootstrapLegacyViewer.js` to host the legacy viewer logic while we build out new providers and components.
- Introduced `ViewerProvider`, `WsProvider`, and `TimelineProvider` contexts plus a lightweight `ReadyBanner` so the React surface can observe viewer + websocket + timeline lifecycle today without disrupting legacy DOM controls; migration of additional providers/components remains TODO.
- Legacy bootstrap now returns cleanup handles for DOM listeners and intervals so the React shell can mount/unmount safely while the DOM UI remains authoritative.
- Pydantic validator refactor postponed until after core frontend work; warnings remain acceptable for now.
- React shell now consumes live websocket/timeline events via new legacy bridge hooks (`viewerApi.ws.addStateListener`, `viewerApi.timeline.subscribe`); providers maintain connection status, countdown metadata, and timeline activity inside React state.
- Added a React `ReconnectBanner` that replaces the DOM-driven banner when `?ui=react` is active (legacy `showReconnectBannerUi` now no-ops in React mode). Banner wires to the shared websocket context and exposes manual reconnect.
- Introduced React Testing Library coverage (`src/react/__tests__/ReconnectBanner.test.tsx`) and extended Jest/Babel config (`@babel/preset-react`, jsdom matchers). `npm run test:fast` passes with the new suite enabled.
- Playwright fixtures now accept a `uiMode` toggle; config defines `legacy` and `react` projects so the same e2e specs can be exercised against both shells. React runs succeed for `ws-atom-edit` and the updated timeline slider/replay specs (warnings about Babylon UMD bundling remain benign and logged for follow-up).
- Added `tests-e2e/utils/timeline.js` to centralise React timeline interactions (buffer waits, offset selection, loop/FPS helpers) so specs no longer poke hidden legacy controls. `ws-timeline-slider-select` and `ws-timeline-replay` were ported to the helper and validate under both UI modes.
- Babylon global now comes from a dedicated shim (`src/babylonLegacy.ts`) that loads the scoped `@babylonjs/core`, loaders, and GUI bundles; the legacy UMD scripts were dropped from `public/index.html`.
- `bootstrapLegacyViewer` accepts `useReactUi` and skips building the DOM control panel + polling timers in React mode, preventing duplicate handlers while the new React `StructureControls` / `SimulationControls` take over those responsibilities.
- React `StructureControls` and `SimulationControls` fully cover cell/ghost toggles, force visibility, and MD/relax start-stop actions; associated DOM/Jest suites (`tests/x-selection-panel.dom.spec.js`, `tests/x-desktop-bond-pick.dom.spec.js`, `tests/x-auto-md-panel-toggle-reset.spec.js`) now pass against the React components.
- `ws-timeline-controls` now runs cleanly against the React dock after syncing `timelineState.playing` with `timelinePlayback.isPlaying()` ahead of every event emission (`public/index.js`). Playwright React runs complete without the prior pause/poll flake.
- Added a `StatusHud` React overlay driven by the shared contexts so websocket connectivity and timeline mode/offset are visible in React mode. Unit tests cover its rendering states and `npm run test:fast` stays green.
- First React-native simulation controls landed (`SimulationControls.tsx`) to mirror the Relax/MD start-stop actions without touching the legacy panel. Buttons use the viewer context for mutation and the timeline context for run-state so React UI can drive simulations. Jest coverage added and the full `npm run test:fast` suite remains green.
- Additional React Playwright runs cover `ws-start-stop-idle-relax`, `ws-live-md-param-updates`, `ws-initial-idle-energy`, `ws-ws-reconnect`, `ws-session-playback-resume`, and `ws-session-save-load`, further validating simulation control surfaces (start/stop) and session persistence flows.
- Session load/save UI now lives in React (`SessionControls.tsx`) with data-testids matching legacy buttons; it wraps `viewerApi.session` for file uploads, downloads, and library manifests. Jest coverage added and legacy hidden controls remain for compatibility.
- Timeline dock (play/pause/live + scrub slider) is now implemented in React (`TimelineDock.tsx`) powered by the timeline context; the legacy HTML timeline is skipped automatically when the React shell is active. Jest coverage ensures the controls delegate to the timeline API, and Playwright selectors continue to work through matching `data-testid`s.
- React `SimulationControls` now drives MD temperature and developer friction sliders, mirroring the legacy inputs while keeping `__MLIP_TARGET_TEMPERATURE` and `__MLIP_CONFIG.mdFriction` synchronised with the viewer API. The RTL suite asserts the new handlers and dev gating.
- React `SessionControls` absorbed the molecule preset picker, SMILES navigation, and XYZ upload flow (with validation + banner errors) so the legacy panel can retire. Tests cover navigation calls, manifest/JSON uploads, and invalid SMILES handling.
- React `XRControls` replaces the DOM XR widget, polling `vr.debugInfo()` for session state and delegating to `switchXR`. Jest coverage ensures mode changes/polling work while `bootstrapLegacyViewer` simply skips the legacy widget in React mode.
- React `LiveMetricsPanel` now reads metrics via the viewer bus (energy, temperature, RPS, max force/stress) and renders the sparkline/marker energy plot that used to live in the legacy panel, including range/sample readouts.
- React `RenderControls` exposes solid/soft atom & bond toggles plus a reset action, wiring through `viewerApi.view.applyMeshModes`/`resetToBaselineModes`; unit tests cover the new API surface.
- Added a `timelineBridge` runtime adapter that feeds the React provider with canonical status snapshots (`playing`, `active`, loop/fps, last events) and exposes async helpers (`ensurePaused`, `live`) so the React dock no longer polls legacy flags.
- Timeline-focused Playwright suites (`ws-timeline-controls`, `ws-timeline-replay`, `ws-timeline-slider-select`, plus helpers in `tests-e2e/utils/timeline.js`) now consume the bridge status and `ensurePaused` helpers rather than polling legacy flags, stabilising React-mode runs.
- Extended the remaining timeline-react Playwright specs (`ws-timeline-camera`, `ws-timeline-visibility`, `ws-timeline-mesh-mode`, `ws-ghost-periodic`) to use the shared helpers and `computeTimelineOffset`, and taught the Playwright global setup/teardown to respect `MLIPVIEW_SKIP_SERVERS=1` by reusing the externally launched backend while still starting/stopping the Vite preview.
- React Playwright runs now share the same CUDA backend when the skip flag is set; `ws-session-save-load` passes end-to-end, but the longer timeline suites still hit the current 60 s harness timeout and `ws-library-sn2-autoplay` stalls waiting for the timeline to surface `mode === 'playing'`. Raw logs live under `tmp/test-logs/` for debugging.
- Advanced timeline controls (`TimelineAdvancedControls.tsx`) now surface loop toggles, range editing, playback FPS overrides, and per-frame metadata; timeline context exposes loop/speed/frame information via new helpers. Unit tests cover the bridge to the playback controller, and the session loader gained robust fallbacks for environments without `File.text()`.
- WebSocket restores now seed both `nextSeq` and `lastSeq`, ensuring the first post-load `USER_INTERACTION` reuses the stored sequence number (`public/fairchem_ws_client.js`, `public/core/sessionStateManager.ts`, `tests/fairchem_ws_client.seedSequencing.spec.js`).
- Snapshot bonds serialize in columnar form (`atomIndexA/B`, offsets, flags, `weight`, `imageDelta`); loaders normalise older array-shaped data and `json_library.md` documents the schema.
- Error banner wiring now routes through the React shell (`ErrorBannerProvider` + `ErrorBanner`), while `public/ui/errorBanner.js` delegates to the provider when active so legacy DOM remains a fallback.
- `SessionStateManager` has been migrated to TypeScript (`public/core/sessionStateManager.ts`), with helper utilities for columnar bond serialisation/deserialisation. Imports/tests target the new module, leaving room for future type refinement.
- Live metrics banner (`LiveMetricsPanel.tsx`) polls `viewerApi.getMetrics()` + viewer state to show energy/force/stress/temperature; covered by fake-timer Jest tests.
- React `SelectionDetails` replaces the legacy selection panel: shares a new `public/data/elementCatalog.js` for names/masses/colors, mirrors atom/bond metadata (positions, weights, vdW, bond rotation) and is covered by RTL tests (`SelectionDetails.test.tsx`).
- Simulation controls now include the Show Forces toggle in React: `SimulationControls.tsx` observes `state.showForces`, calls `setForceVectorsEnabled`, and `SimulationControls.test.tsx` verifies UI/WS parity.
- New `StructureControls.tsx` handles unit-cell and ghost-cell visibility via the React shell, keeps them in sync with `moleculeState` events, and supersedes the legacy periodic panel wiring when `?ui=react` is active; covered by `StructureControls.test.tsx`.
- Playwright global setup respects `MLIPVIEW_SKIP_SERVERS=1`, letting us reuse an already running backend/preview when iterating locally (set `BASE_URL` if it isn’t the default 5174).
- React Playwright `ws-atom-edit` currently stalls in the idle geometry loop (Chromium pegs CPU with no further console output). Need to tighten the test’s ready check once the React shell is active.
- Continue porting the remaining Playwright suites (legacy render toggles, session/timeline flows) to React selectors, then remove the hidden DOM shims once parity is proven.
- Playwright specs still target legacy DOM selectors for render toggles/energy plot; add React-mode variants once the new components are exercised, then delete the hidden legacy buttons (`legacyControlsHidden`).
- `ws-live-md-param-updates.spec.js` now detects React UI mode and drives the new temperature slider / metrics HUD via `data-testid` selectors so the test runs against both panels.
- Next focus: finish migrating the remaining Playwright suites (render toggles, session replay smoke) to React-aware helpers, then remove the hidden DOM shims once those tests stay green.

## 12. Design Proposals (2025-11-03)

## 12. Current Issues & Next Steps (2025-11-03)
- Playwright React timeline suites (`ws-timeline-camera`, `ws-timeline-visibility`, `ws-timeline-mesh-mode`, `ws-ghost-periodic`) now use the shared helpers but still exceed the hard 60 s harness timeout. Decide whether to extend the cap or trim the scenarios; the captured traces live in `tmp/test-logs/ws-*.react.log`.
- `ws-library-sn2-autoplay` consistently times out waiting for `timeline.getStatus().mode === 'playing'` even though frames stream in. Need to trace the bridge status flow and surface the playback mode before the wait expires (log snapshot in `tmp/test-logs/ws-library-sn2-autoplay.react.log`).
- Hidden legacy controls remain in `public/ui/desktopPanel.js` until the React render toggles and metrics tests finish porting; once the React suites are green, remove the shims and retire the DOM-only specs.
- Document the `MLIPVIEW_SKIP_SERVERS=1` workflow in `testing.md` so local developers understand how to share the CUDA backend while running React Playwright projects; update CI notes if we adopt a longer timeout.
- After the outstanding suites stabilise, flip the default Playwright project to React mode and keep the legacy project only for regression coverage until the DOM panel is formally removed.

### React Shell Blueprint
1. **Bootstrapping & Routing**
   - Default path keeps loading `public/bootstrapLegacyViewer.js`. Opt-in React mode uses `?ui=react` to dynamically import `/src/main.tsx` which mounts `<App />` via `createRoot`.
   - Both paths share the same `#app` container so we can toggle between implementations during migration and for Playwright coverage.

2. **App Layers**
   - `App` renders the Babylon canvas + status nodes, then calls a *viewer bridge hook* that encapsulates:
     - Initialisation (calling the legacy bootstrap today).
     - Exposing the viewer API through React context providers once available.
     - Cleanup (`viewerApi.shutdown`, interval disposal, DOM handler restore) when the component unmounts or reloads.
   - Initial contexts to introduce:
     1. `ViewerContext` – wraps the legacy `viewerApi`, state snapshots, and mutators.
     2. `WsContext` – surfaces connection state, counters, and helper actions (`startSimulation`, `stopSimulation`, etc.).
     3. `TimelineContext` – exposes frame buffer metadata, playback state, and control-message APIs.
   - Each context should include a hook (e.g., `useViewer`, `useWsClient`, `useTimeline`) that describes the contract React components will consume later.

3. **Incremental React Components**
   - Phase A: Render read-only overlays (status banner, reconnect toast, test HUD) using React components fed by the contexts. These components live alongside the legacy DOM to validate data flow without taking over controls yet.
   - Phase B: Rebuild control surfaces (desktop panel sections, timeline dock) as React components, wired to the contexts and removing their legacy DOM counterparts.
   - Phase C: Wrap VR overlays/portals via refs to the Babylon GUI, keeping the render loop unaffected.

4. **Legacy Bridge Evolution**
   - During early Stage 2 the bridge simply mounts the legacy bootstrap; as React replacements land, progressively move logic out of `bootstrapLegacyViewer` into typed hooks:
     - Extract viewer initialisation into `useBabylonViewer` hook.
     - Promote session management helpers into a React provider, letting contexts own snapshot/load/reset.
     - Once all UI functionality lives in React components, retire the bridge and switch the default `index.html` boot path to React.

5. **Testing Strategy**
   - Add React Testing Library coverage for provider contracts and new components (e.g., context hooks, timeline panel).
   - Extend Playwright to run both `?ui=legacy` (default) and `?ui=react` suites while in transition; ensure shared fixtures use data-testids so selectors stay stable.

## 11. Review Concerns (2025-11-01)
- **Snapshot schema hard downgrade** *(resolved 2025-11-02)* – `public/core/sessionStateManager.ts:266` now rejects snapshots with a `schemaVersion` greater than 6, preventing silent data loss when future formats appear.
- **Jest exit strategy masking async leaks** *(mitigated 2025-11-02)* – `jest.config.js` keeps `forceExit`/`detectOpenHandles` by default but honours `MLIPVIEW_JEST_FORCEEXIT=0` so local runs can opt out while CI remains protected.
- **Babylon dependency duplication** *(resolved 2025-11-01)* – migrated to the scoped `@babylonjs/*` bundles and load them via an explicit `babylonLegacy` shim so we no longer ship the legacy UMD copies alongside the ESM build.
- **Legacy VR / cell tests still failing** – `test_to_port:39-44` documents `tests/vrPostReleaseLatch.spec.js` and `tests-browser/aceticAcidGhostsAndApi.dom.spec.js` as failing because they still depend on the deprecated REST wiring. These leave gaps in the VR joystick latch regression and periodic ghost parity, so we should decide whether to port or replace them during Stage 2.
- **Legacy XR dropdown tests to port** – `tests/x-xr-dropdown.dom.spec.js` and related suites still target `#xrModeSelect`. Once the React XR controls stabilise we should clone those assertions into RTL/Playwright coverage and retire the DOM harness.

## 13. Current Status (2025-01-09)
- React shell now initialises the viewer through `useBabylonViewer` (`src/react/runtime/useBabylonViewer.ts`) and wires the session manager via the new `SessionStateProvider` (`src/react/context/SessionStateContext.tsx`), while legacy bootstrap remains available for the non-React entry.
- Simulation controls call the typed viewer API (`viewerApi.setTargetTemperature`, `viewerApi.setMDFriction`) and honour the timeline lock. Global temperature/friction flags remain in sync for legacy code paths, but React no longer mutates them directly (`src/react/components/SimulationControls.tsx`).
- Legacy REST-era tests listed in `test_to_port` have either been replaced by their `tests/x-*.spec.js` WS equivalents or formally retired with equivalent WS coverage noted.
- Timeline overlay/interaction lock is owned by React via `TimelineOverlay.tsx`; the legacy DOM overlay is still instantiated when `window.__MLIPVIEW_REACT_UI_ACTIVE` is false so non-React builds retain previous behaviour.
