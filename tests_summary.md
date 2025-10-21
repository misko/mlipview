# UI Tests Summary

Generated: 2025-10-20T15:13:57.044Z

| File | Area | Test summary | Design alignment | Currently passing |
| --- | --- | --- | --- | --- |
| tests-browser/aceticAcidCellLoad.dom.spec.js | pbc/cell | acetic acid cell loads and PBC enabled | Neutral (non-protocol UI) | Pass |
| tests-browser/aceticAcidGhostsAndApi.dom.spec.js | ui/core | acetic acid: ghosts render and API includes cell | Legacy REST usage (needs migration) | Fail |
| tests-browser/aceticAcidPeriodicToggle.dom.spec.js | ui/core | Periodic UI toggle reflects PBC On after XYZ load with cell | Neutral (non-protocol UI) | Pass |
| tests-browser/energyPlot.dom.spec.js | energy/plot | energy plot DOM integration | Neutral (non-protocol UI) | Pass |
| tests-browser/highlightSphere.dom.spec.js | ui/core | highlight sphere hidden initially | Neutral (non-protocol UI) | Pass |
| tests-browser/loadAppliesTemperatureAndCell.dom.spec.js | pbc/cell | applyParsedToViewer temperature and cell handling | Neutral (non-protocol UI) | Pass |
| tests-browser/mdInitialRequestsTemperature.dom.spec.js | md/simulation | MD initial requests use XYZ temperature when present | Neutral (non-protocol UI) | Fail |
| tests-browser/mdInitialTemperatureRequests.dom.spec.js | md/simulation | MD initial requests honor XYZ temperature | Neutral (non-protocol UI) | Fail |
| tests-browser/methylTemperaturePrefill.dom.spec.js | ui/core | methyl radicals example initializes temperature to 500K | Neutral (non-protocol UI) | Pass |
| tests-browser/temperatureSliderPrefillFromXYZ.dom.spec.js | xyz/parser | temperatureSliderPrefillFromXYZ.dom.spec.js | Neutral (non-protocol UI) | Pass |
| tests-browser/vrBondJoystickRotation.dom.spec.js | vr/xr | VR bond joystick rotation | Neutral (non-protocol UI) | Pass |
| tests-browser/ws_counters_and_drag.spec.js | websocket/protocol | WS userInteraction counter gating and drag exclusion | Legacy REST usage (needs migration) | Pass |
| tests-browser/xyzCommentParsing.dom.spec.js | xyz/parser | XYZ comment parsing for temperature and cell | Neutral (non-protocol UI) | Pass |
| tests-e2e/autoMd0-baseline-energy.ws.spec.js | websocket/protocol | autoMD=0 baseline simple_calculate over real WS | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/autoMd0-drag-energy.ws.spec.js | websocket/protocol | autoMD=0 baseline and drag energy updates (real WS) | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/autoMdRoy.integration.spec.js | md/simulation | [integration] auto MD on ROY over real WS | Neutral (non-protocol UI) | Unknown |
| tests-e2e/benzene-selection.spec.js | selection | Benzene selection / highlight cleanliness | Neutral (non-protocol UI) | Unknown |
| tests-e2e/drag-throttle-emit.spec.js | ui/core | drag-throttle-emit.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/energy-drag.spec.js | energy/plot | Energy plot updates on atom drag (full page) | Neutral (non-protocol UI) | Unknown |
| tests-e2e/energyApiOnlyTicks.integration.spec.js | energy/plot | [integration] API-only energy ticks (live WS) | Neutral (non-protocol UI) | Unknown |
| tests-e2e/forceCache.integration.spec.js | forces/render/cache | [integration] force cache versioning (live WS) | Neutral (non-protocol UI) | Unknown |
| tests-e2e/idle-drag-energy-ticks-accept-stale.spec.js | energy/plot | idle-drag-energy-ticks-accept-stale.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/idle-drag-energy-updates.spec.js | energy/plot | idle-drag-energy-updates.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/md-drag-uic-echo.spec.js | md/simulation | md-drag-uic-echo.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/mdVelocityContinuity.integration.spec.js | md/simulation | [integration] md velocity continuity (live WS) | Neutral (non-protocol UI) | Unknown |
| tests-e2e/relax-drag-latch.spec.js | relaxation | relax-drag-latch.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/relaxForceCache.integration.spec.js | relaxation | [integration] relax + force cache (live WS) | Neutral (non-protocol UI) | Unknown |
| tests-e2e/royMultiStepRelaxParity.spec.js | relaxation | royMultiStepRelaxParity.spec.js | Mixed (migrate to WS-only) | Unknown |
| tests-e2e/smoke-page-load-drag.spec.js | ui/core | smoke-page-load-drag.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/smoke-page-load-md.spec.js | md/simulation | smoke-page-load-md.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/smoke-page-load.spec.js | ui/core | smoke-page-load.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/uiInitialForcesParity.spec.js | forces/render/cache | uiInitialForcesParity.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/uiRelaxParity.spec.js | relaxation | uiRelaxParity.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/waterDirectRelax.spec.js | relaxation | waterDirectRelax.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/waterFairChemRelaxation.spec.js | relaxation | waterFairChemRelaxation.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/waterMultiStepRelaxParity.spec.js | relaxation | waterMultiStepRelaxParity.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/waterRelaxation.spec.js | relaxation | waterRelaxation.spec.js | Neutral (non-protocol UI) | Unknown |
| tests-e2e/ws-backpressure-ack.spec.js | websocket/protocol | WS protocol: backpressure and ACK | WS aligned (protobuf) | Unknown |
| tests-e2e/ws-cell-stress.spec.js | websocket/protocol | WS protocol: cell and optional stress | WS aligned (protobuf) | Unknown |
| tests-e2e/ws-counters-echo.spec.js | websocket/protocol | ws-counters-echo.spec.js | WS aligned (protobuf) | Unknown |
| tests-e2e/ws-idle-positions-omitted.spec.js | websocket/protocol | ws-idle-positions-omitted.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/ws-mock-load.spec.js | websocket/protocol | WS mock page load | Neutral (non-protocol UI) | Unknown |
| tests-e2e/ws-sim-counters.spec.js | websocket/protocol | WS protocol: simulation counters | WS aligned (protobuf) | Unknown |
| tests-e2e/ws-start-stop-idle-relax.spec.js | websocket/protocol | ws-start-stop-idle-relax.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/ws-stream-single-socket.spec.js | websocket/protocol | WS streaming single-socket | Neutral (non-protocol UI) | Unknown |
| tests-e2e/ws-ui-buttons-start-stop-idle-relax.spec.js | websocket/protocol | ws-ui-buttons-start-stop-idle-relax.spec.js | Legacy REST usage (needs migration) | Unknown |
| tests-e2e/xr-dropdown.spec.js | vr/xr | XR dropdown presence | Neutral (non-protocol UI) | Unknown |
| tests/apiCellPayload.spec.js | pbc/cell | API includes cell when PBC enabled | WS aligned (protobuf) | Pass |
| tests/arAutoNormalize.spec.js | ui/core | ar auto-normalization | Neutral (non-protocol UI) | Pass |
| tests/atomDragCamera.spec.js | md/simulation | atom drag does not move camera | Neutral (non-protocol UI) | Pass |
| tests/atomDragDebug.spec.js | md/simulation | atomDragDebug placeholder | Neutral (non-protocol UI) | Pass |
| tests/atomDragScreenPlane.spec.js | md/simulation | camera-aligned drag plane | Neutral (non-protocol UI) | Pass |
| tests/autoMdDisabledBaseline.dom.spec.js | md/simulation | autoMD=0 baseline energy | WS aligned (protobuf) | Pass |
| tests/autoMdDisabledInteractionEnergies.dom.spec.js | md/simulation | autoMD=0 interaction energies | Mixed (migrate to WS-only) | Fail |
| tests/autoMdPanelToggleReset.dom.spec.js | md/simulation | Auto MD Panel Toggle Reset | Neutral (non-protocol UI) | Pass |
| tests/autoMdRunReset.spec.js | md/simulation | autoMdRunReset.spec.js | Neutral (non-protocol UI) | Pass |
| tests/autoMdStartsRoy.ws.spec.js | websocket/protocol | auto MD on page load (ROY) | WS aligned (protobuf) | Pass |
| tests/base64Utf8.spec.js | ui/core | UTF-8 base64 helpers and XYZ parsing | Neutral (non-protocol UI) | Pass |
| tests/benzeneBondRotationLengths.spec.js | bond/rotate | Benzene bond rotation keeps reasonable bond lengths (with cell + ghosts) | Neutral (non-protocol UI) | Pass |
| tests/bondHighlight.spec.js | bond/rotate | bond highlight selection | Neutral (non-protocol UI) | Pass |
| tests/bondHighlightRotation.spec.js | bond/rotate | bond highlight rotation inheritance | Neutral (non-protocol UI) | Pass |
| tests/bondOpacityIsolation.spec.js | bond/rotate | bond opacity isolation when displacing a single carbon | Neutral (non-protocol UI) | Pass |
| tests/bondRotateEnergyCount.spec.js | energy/plot | single bond rotation energy step count (API-only policy) | Neutral (non-protocol UI) | Pass |
| tests/bondRotationSoftBond.spec.js | bond/rotate | bond rotation exclusion of soft bonds | Neutral (non-protocol UI) | Pass |
| tests/bondRotationUI.spec.js | bond/rotate | bond rotation manipulation | Neutral (non-protocol UI) | Pass |
| tests/bondService.spec.js | bond/rotate | bondService periodic | Neutral (non-protocol UI) | Pass |
| tests/bond_recompute_toggle.spec.js | bond/rotate | bond recompute on manual position changes | Legacy REST usage (needs migration) | Pass |
| tests/browserPick.spec.js | websocket/protocol | browserPick.spec.js | Neutral (non-protocol UI) | Pass |
| tests/cacheReuse.spec.js | ui/core | cacheReuse.spec.js | Legacy REST usage (needs migration) | skipped |
| tests/cameraSuppression.spec.js | ui/core | camera suppression during atom drag | Neutral (non-protocol UI) | Pass |
| tests/cellOptimizer.spec.js | pbc/cell | cell optimizer | Neutral (non-protocol UI) | Pass |
| tests/cellToggle.dom.spec.js | pbc/cell | Cell toggle (unified) | Neutral (non-protocol UI) | Pass |
| tests/debugLatencyInterRequest.dom.spec.js | ui/core | debugLatencyInterRequest.dom.spec.js | Neutral (non-protocol UI) | Pass |
| tests/desktopArcRotateMouseOrbit.dom.spec.js | ui/core | desktopArcRotateMouseOrbit.dom.spec.js | Neutral (non-protocol UI) | skipped |
| tests/desktopBondPick.dom.spec.js | bond/rotate | desktop bond selection via picking | Neutral (non-protocol UI) | Pass |
| tests/desktopPanel.dom.spec.js | ui/core | desktop left panel UI | Neutral (non-protocol UI) | Pass |
| tests/desktopRotationInputPreserved.dom.spec.js | ui/core | desktop rotation inputs remain intact | Neutral (non-protocol UI) | Pass |
| tests/dragActiveExclusion.spec.js | ui/core | drag exclusion from incoming sim updates | Neutral (non-protocol UI) | Fail |
| tests/dragEndResumesSimApply.dom.spec.js | ui/core | drag end resumes sim application | Neutral (non-protocol UI) | Fail |
| tests/dragRecompute.spec.js | ui/core | atom drag triggers bond recompute on end | Neutral (non-protocol UI) | Pass |
| tests/dynamics.spec.js | ui/core | dynamics.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/energy-plot.spec.js | energy/plot | energy-plot.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/energyApiOnlyTicks.spec.js | energy/plot | API-only energy ticks | WS aligned (protobuf) | Pass |
| tests/energyInteractions.spec.js | energy/plot | energy plot tick policy (API-only) | Neutral (non-protocol UI) | Pass |
| tests/energyPlot.spec.js | energy/plot | energy plot integration (API-only ticks) | Neutral (non-protocol UI) | Pass |
| tests/energyReset.dom.spec.js | energy/plot | Energy plot clears on reset | Neutral (non-protocol UI) | Fail |
| tests/engine/nullEngine.spec.js | ui/core | NullEngine harness basic | Neutral (non-protocol UI) | Pass |
| tests/fairchem_bfgs_parity.spec.js | ui/core | fairchem water BFGS parity via server /relax | Legacy REST usage (needs migration) | Pass |
| tests/forceArrowheadRender.spec.js | forces/render/cache | force arrowheads | Neutral (non-protocol UI) | Pass |
| tests/forceCache.spec.js | forces/render/cache | force cache: version increments only on geometry changes (WS idle compute) | WS aligned (protobuf) | Pass |
| tests/forceProvider.spec.js | forces/render/cache | force provider abstraction | Neutral (non-protocol UI) | Pass |
| tests/forceRender.spec.js | forces/render/cache | force rendering | Neutral (non-protocol UI) | Pass |
| tests/forceUpdatePerturbation.spec.js | forces/render/cache | force update after perturbation | Neutral (non-protocol UI) | Fail |
| tests/forcesAsyncRefresh.spec.js | forces/render/cache | forcesAsyncRefresh placeholder | Neutral (non-protocol UI) | Pass |
| tests/forcesToggle.dom.spec.js | forces/render/cache | forces toggle button | Neutral (non-protocol UI) | Pass |
| tests/ghostBonds.spec.js | bond/rotate | ghost bonds render when cell and ghosts enabled | Neutral (non-protocol UI) | Pass |
| tests/ghostCrossingBonds.spec.js | bond/rotate | cross-image ghost bonds are generated | Neutral (non-protocol UI) | Pass |
| tests/ghostSelection.spec.js | selection | ghost atoms/bonds are not pickable | Neutral (non-protocol UI) | Pass |
| tests/highlight.spec.js | ui/core | selection highlight | Neutral (non-protocol UI) | Pass |
| tests/highlightStrayCylinder.spec.js | ui/core | highlight stray bond cylinder regression | Neutral (non-protocol UI) | Pass |
| tests/highlightVisibility.spec.js | ui/core | highlight visibility | Neutral (non-protocol UI) | Pass |
| tests/httpAppFactory.spec.js | ui/core | server-app factory | Neutral (non-protocol UI) | Pass |
| tests/httpVrAsset.spec.js | vr/xr | HTTP VR asset | Neutral (non-protocol UI) | Fail |
| tests/interactionInvalidation.spec.js | ui/core | interaction invalidation | Neutral (non-protocol UI) | Fail |
| tests/joystickDelta.spec.js | ui/core | joystickDelta | Neutral (non-protocol UI) | Pass |
| tests/lightingConsistency.spec.js | ui/core | lighting consistency | Neutral (non-protocol UI) | Pass |
| tests/line-plot-core.spec.js | ui/core | line-plot-core | Neutral (non-protocol UI) | Pass |
| tests/lj_bfgs_parity.spec.js | ui/core | LJ + BFGS parity via server /relax | Legacy REST usage (needs migration) | Pass |
| tests/manipulation.spec.js | ui/core | manipulationService | Neutral (non-protocol UI) | Pass |
| tests/mdDragExclusionDesktop.spec.js | md/simulation | desktop drag excludes MD updates for held atom | Legacy REST usage (needs migration) | Fail |
| tests/mdDragExclusionVR.spec.js | md/simulation | VR drag excludes MD updates for held atom | Legacy REST usage (needs migration) | Fail |
| tests/mdDragPositionLock.dom.spec.js | md/simulation | Dragged atom position lock (client-side semantics) | WS aligned (protobuf) | Fail |
| tests/mdForcesUpdate.spec.js | md/simulation | md forces visualization | WS aligned (protobuf) | Fail |
| tests/mdFrictionSlider.spec.js | md/simulation | MD friction slider | Neutral (non-protocol UI) | Fail |
| tests/mdInstantTemperatureDisplay.spec.js | md/simulation | Instantaneous MD temperature HUD | WS aligned (protobuf) | Pass |
| tests/mdRpsLabel.ws.spec.js | websocket/protocol | MD RPS label (WS) | WS aligned (protobuf) | Pass |
| tests/mdTemperatureDynamicRun.spec.js | md/simulation | MD run dynamic temperature | Neutral (non-protocol UI) | Pass |
| tests/mdTemperatureSlider.spec.js | md/simulation | MD temperature slider (WS) | WS aligned (protobuf) | Fail |
| tests/mdTemperatureSync.dom.spec.js | md/simulation | Temperature slider sync | WS aligned (protobuf) | Pass |
| tests/mdVelocityContinuity.dom.spec.js | md/simulation | md velocity continuity (jsdom) | WS aligned (protobuf) | Fail |
| tests/mdVelocityContinuity.spec.js | md/simulation | md velocity continuity | WS aligned (protobuf) | Pass |
| tests/mobileCameraDetachOnTouch.dom.spec.js | mobile | mobile: camera detach/reattach during touch gestures | Neutral (non-protocol UI) | Pass |
| tests/mobileGestures.dom.spec.js | mobile | mobile gestures | Neutral (non-protocol UI) | Pass |
| tests/mobileMenu.dom.spec.js | mobile | mobile top bar tabs | Neutral (non-protocol UI) | Pass |
| tests/mobileResponsivePanel.dom.spec.js | mobile | responsive: desktop panel vs mobile top bar | Neutral (non-protocol UI) | Pass |
| tests/mobileRotateNoZoom.babylonInput.dom.spec.js | mobile | mobile: custom touch controls suppress default pointer zoom | Neutral (non-protocol UI) | Pass |
| tests/mobileRotateNoZoom.dom.spec.js | mobile | mobile: rotate does not zoom | Neutral (non-protocol UI) | Pass |
| tests/mobileSingleTouchOrbit.dom.spec.js | mobile | mobile: single-finger orbit rotation | Neutral (non-protocol UI) | Pass |
| tests/mobileSingleTouchZoomRepro.dom.spec.js | mobile | mobile: single-finger should not zoom (bug repro) | Neutral (non-protocol UI) | Pass |
| tests/mobileTouchAtomDrag.dom.spec.js | md/simulation | mobileTouchAtomDrag.dom.spec.js | Neutral (non-protocol UI) | skipped |
| tests/mobileViewerWiring.dom.spec.js | mobile | mobileViewerWiring.dom.spec.js | Neutral (non-protocol UI) | skipped |
| tests/moleculeReload.spec.js | molecule/load/state | molecule reload selection flow | Neutral (non-protocol UI) | Pass |
| tests/moleculeState.pbc.spec.js | pbc/cell | moleculeState PBC integration | Neutral (non-protocol UI) | Pass |
| tests/moleculeState.spec.js | molecule/load/state | moleculeState | Neutral (non-protocol UI) | Pass |
| tests/moleculeSwitchHighlight.spec.js | molecule/load/state | molecule switch clears highlight and hides empty masters | Neutral (non-protocol UI) | Pass |
| tests/molecules.spec.js | molecule/load/state | molecule XYZ files | Neutral (non-protocol UI) | Pass |
| tests/no-legacy-imports.spec.js | ui/core | no legacy import references inside mlipviewer2 | Neutral (non-protocol UI) | Pass |
| tests/no-legacy-meshes.spec.js | ui/core | no legacy mesh artifacts | Neutral (non-protocol UI) | Pass |
| tests/noCacheFlag.spec.js | ui/core | noCache debug flag forces legacy endpoints | Legacy REST usage (needs migration) | Fail |
| tests/noStrayOriginPrimitives.spec.js | ui/core | no stray origin primitives | Neutral (non-protocol UI) | Pass |
| tests/optimizer.spec.js | ui/core | optimizer core | Neutral (non-protocol UI) | Pass |
| tests/pageLoad.spec.js | ui/core | index.html smoke load | Neutral (non-protocol UI) | Pass |
| tests/panelRunCompletionProgrammatic.dom.spec.js | ui/core | Panel starts full runs | Neutral (non-protocol UI) | Fail |
| tests/pbc.spec.js | pbc/cell | PBC utilities | Neutral (non-protocol UI) | Pass |
| tests/pbcGhostMasterDisabledOnBondBreak.dom.spec.js | pbc/cell | PBC ghost bond masters disable when bonds disappear | Neutral (non-protocol UI) | Pass |
| tests/pbcLongBondCrossing.dom.spec.js | pbc/cell | PBC long primary bond crossing (10x10x10 cell, C @ 0.1 and 9.9) | Neutral (non-protocol UI) | Pass |
| tests/pbcLongBondCrossingX.dom.spec.js | pbc/cell | PBC along-X primary long bond suppression (10x10x10, C @ 0.1 and 9.9 on X) | Neutral (non-protocol UI) | Pass |
| tests/periodicMonoclinicUI.dom.spec.js | ui/core | Periodic tab monoclinic UI | Neutral (non-protocol UI) | Pass |
| tests/picking.spec.js | ui/core | picking service (logic only) | Neutral (non-protocol UI) | Pass |
| tests/pickingIntegration.spec.js | ui/core | picking integration (selection + view) | Neutral (non-protocol UI) | Pass |
| tests/pointerObservableBug.spec.js | ui/core | pointerObservableBug.spec.js | Neutral (non-protocol UI) | Pass |
| tests/precomputedInjection.spec.js | ui/core | precomputed injection (relax + md) | Legacy REST usage (needs migration) | Fail |
| tests/precomputedResults.spec.js | ui/core | precomputed results seeding | Legacy REST usage (needs migration) | Fail |
| tests/redSphereSwitch.spec.js | ui/core | red sphere artifact regression | Neutral (non-protocol UI) | Pass |
| tests/relaxForceCache.spec.js | relaxation | relax + force cache integration (WS) | WS aligned (protobuf) | Pass |
| tests/relaxForcesUpdate.spec.js | relaxation | relax forces visualization (WS) | WS aligned (protobuf) | Pass |
| tests/relaxSingleStepNetwork.spec.js | relaxation | network: single relax step (WS) | Mixed (migrate to WS-only) | Pass |
| tests/relaxWaterIntegration.spec.js | relaxation | relaxWaterIntegration.spec.js | Legacy REST usage (needs migration) | Fail |
| tests/relax_run_post_init_enable.spec.js | relaxation | post-init feature flag enable | Neutral (non-protocol UI) | Pass |
| tests/relaxationEdgeCases.spec.js | relaxation | relaxationEdgeCases.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/relaxationFairChem.spec.js | relaxation | relaxationFairChem.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/relaxationRunner.spec.js | relaxation | relaxationRunner.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/requestThrottlingLatency.spec.js | ui/core | request pacing and latency breakdown | Legacy REST usage (needs migration) | Fail |
| tests/resetButton.dom.spec.js | ui/core | Reset button | Neutral (non-protocol UI) | Pass |
| tests/resetInvalidation.dom.spec.js | ui/core | Reset invalidation and no reload | Neutral (non-protocol UI) | Fail |
| tests/rotateRecompute.spec.js | ui/core | rotateBond recompute | Neutral (non-protocol UI) | Pass |
| tests/rotationInversion.spec.js | ui/core | rotation inversion regression | Neutral (non-protocol UI) | Pass |
| tests/runCompletionToggles.dom.spec.js | ui/core | Run completion flips toggles Off | Neutral (non-protocol UI) | Fail |
| tests/selectionClear.spec.js | selection | selection clear via empty click | Neutral (non-protocol UI) | Pass |
| tests/selectionNullGuard.spec.js | selection | selectionService null selection robustness | Neutral (non-protocol UI) | Pass |
| tests/selectionPanel.dom.spec.js | selection | Selection panel UI | Neutral (non-protocol UI) | Pass |
| tests/selectionPeriodicClick.dom.spec.js | selection | Periodic table click behavior | Neutral (non-protocol UI) | Pass |
| tests/selectionPersistence.spec.js | selection | selection persistence | Neutral (non-protocol UI) | Pass |
| tests/selectionService.spec.js | selection | selectionService atom/bond mutual exclusivity | Neutral (non-protocol UI) | Pass |
| tests/smilesAndUpload.dom.spec.js | ui/core | System panel SMILES + XYZ upload | Neutral (non-protocol UI) | Pass |
| tests/smilesLoader.unit.spec.js | ui/core | smilesLoader basic | Neutral (non-protocol UI) | Pass |
| tests/sphericalRadialModes.spec.js | ui/core | spherical-drag-core computeRadialDelta | Neutral (non-protocol UI) | Pass |
| tests/stress.spec.js | ui/core | stress.spec (deprecated) | Neutral (non-protocol UI) | Pass |
| tests/temperatureTickPosition.spec.js | ui/core | temperatureTickPosition.spec.js | Neutral (non-protocol UI) | Pass |
| tests/uiPbcBondBreakNoArtifact.dom.spec.js | pbc/cell | UI: PBC enabled, breaking bond does not leave origin artifact | Neutral (non-protocol UI) | Pass |
| tests/uiToggles.dom.spec.js | ui/core | uiToggles.dom.spec.js | Neutral (non-protocol UI) | Pass |
| tests/userInteractionDuringMd.dom.spec.js | md/simulation | USER_INTERACTION during MD | Mixed (migrate to WS-only) | Pass |
| tests/userInteractionVersionIsolation.spec.js | ui/core | userInteractionVersion isolation | Legacy REST usage (needs migration) | Fail |
| tests/userInteractionWhenIdle.dom.spec.js | ui/core | USER_INTERACTION when idle | WS aligned (protobuf) | Pass |
| tests/visualIntegration.spec.js | ui/core | visual integration atoms & bonds | Neutral (non-protocol UI) | Pass |
| tests/vr.setup.spec.js | vr/xr | VR scaffolding | Neutral (non-protocol UI) | Pass |
| tests/vrAtomDragInitiation.spec.js | md/simulation | vr atom drag initiation | Neutral (non-protocol UI) | Pass |
| tests/vrControllerAttach.spec.js | vr/xr | vrControllerAttach placeholder | Neutral (non-protocol UI) | Pass |
| tests/vrDragStability.spec.js | vr/xr | vr drag stability | Neutral (non-protocol UI) | Pass |
| tests/vrDynamicImport.spec.js | vr/xr | VR dynamic import | Neutral (non-protocol UI) | Pass |
| tests/vrHighlight.spec.js | vr/xr | vrHighlight placeholder | Neutral (non-protocol UI) | Pass |
| tests/vrHighlightAlignment.spec.js | vr/xr | vr atom highlight alignment | Neutral (non-protocol UI) | Pass |
| tests/vrHudEnergyPosition.dom.spec.js | energy/plot | VR HUD energy position smoke | Neutral (non-protocol UI) | Pass |
| tests/vrLaserThickness.spec.js | vr/xr | vrLaserThickness.spec.js | Neutral (non-protocol UI) | skipped |
| tests/vrLateBondRotation.spec.js | vr/xr | vr late bond rotation | Neutral (non-protocol UI) | Pass |
| tests/vrModuleAvailability.spec.js | vr/xr | VR module availability | Neutral (non-protocol UI) | Pass |
| tests/vrPostReleaseLatch.spec.js | vr/xr | VR joystick drag release does not snap back | Legacy REST usage (needs migration) | Fail |
| tests/vrReuseParity.spec.js | vr/xr | VR reuse parity | Neutral (non-protocol UI) | Pass |
| tests/vrTrigger.spec.js | vr/xr | vrTrigger.spec.js | Neutral (non-protocol UI) | Pass |
| tests/water_md_lj.spec.js | md/simulation | LJ MD water | Legacy REST usage (needs migration) | Fail |
| tests/water_md_lj_stability.spec.js | md/simulation | LJ MD 250-step stability | Legacy REST usage (needs migration) | Fail |
| tests/water_md_run_stability.spec.js | md/simulation | continuous MD run stability | Neutral (non-protocol UI) | Pass |
| tests/water_md_step_spam.spec.js | md/simulation | water mdStep spam stability | Neutral (non-protocol UI) | Pass |
| tests/water_md_uma.spec.js | md/simulation | UMA MD water | Neutral (non-protocol UI) | Fail |
| tests/water_md_uma_stability.spec.js | md/simulation | UMA MD 250-step stability | Neutral (non-protocol UI) | Fail |
| tests/water_relax_run_parity.spec.js | relaxation | continuous relax run parity | Legacy REST usage (needs migration) | Pass |
| tests/water_relaxation_browser_parity.spec.js | websocket/protocol | water_relaxation_browser_parity.spec.js | Legacy REST usage (needs migration) | Pass |
| tests/wsDoubleConnect.dom.spec.js | websocket/protocol | [mock-browser] WS double connect guard | WS aligned (protobuf) | Pass |
| tests/xrControls.core.spec.js | vr/xr | XR Controls Core Model | Neutral (non-protocol UI) | Pass |
| tests/xrDropdown.dom.spec.js | vr/xr | XR dropdown -> VR API wiring (jsdom) | Neutral (non-protocol UI) | Pass |
| tests/xrEnergyHud.dom.spec.js | energy/plot | XR HUD energy panel | Neutral (non-protocol UI) | Fail |
| tests/xrEnergyPlot.dom.spec.js | energy/plot | sanity placeholder for XR energy plot suite | Neutral (non-protocol UI) | Pass |
| tests/xrHud.dom.spec.js | vr/xr | XR HUD button semantics (simulated) | Neutral (non-protocol UI) | Fail |
| tests/xrHudButtons.spec.js | vr/xr | XR HUD buttons styling | Neutral (non-protocol UI) | Pass |
| tests/xrHudEnergyDepth.spec.js | energy/plot | XR HUD energy depth scaling | Neutral (non-protocol UI) | Pass |
| tests/xrHudEnergyPlotPicking.spec.js | energy/plot | XR HUD energy plot picking & distance | Neutral (non-protocol UI) | Pass |
| tests/xrHudEnergyScale.spec.js | energy/plot | XR HUD top energy scale multiplier | Neutral (non-protocol UI) | Pass |
| tests/xrHudPosition.spec.js | vr/xr | XR HUD vertical positioning | Neutral (non-protocol UI) | Pass |
| tests/xrHudTextureClamp.spec.js | vr/xr | XR HUD texture clamping | Neutral (non-protocol UI) | Pass |
| tests/xyzExtended.spec.js | xyz/parser | Extended XYZ metadata parsing | Neutral (non-protocol UI) | Pass |
| tests/xyzLoader.spec.js | xyz/parser | XYZ loader | Neutral (non-protocol UI) | Pass |
