# Release V1 Fix Proposal

## Findings From The Review (Steps 1 & 2)
- **Bond/Timestep Rendering:** `createBondService` performs periodic bond detection with fractional wrapping, but `molState.bonds` stores only `{ i, j, opacity }`, so downstream renderers lose the periodic metadata needed to rebuild ghost links (`public/domain/bondService.js`, `public/render/moleculeView.js`).
- **Ghost Rendering Pipeline:** Ghost atoms/bonds are generated inside `moleculeView.rebuildGhosts()`. The current augmentation only samples 7 lattice shifts (0, ±a, ±b, ±c) and depends on a fresh `computeBondsNoState` run, which can miss diagonal minimum-image connections once atoms move (`public/render/moleculeView.js:1008‑1160`).
- **State Surfaces:** The frontend routes all mutation through `SessionStateManager` and `stateStore`, while the backend keeps authoritative counters in `SessionState`. Snapshot schema v6 captures mesh assignments, playback config, and WebSocket sequencing, but no field records periodic bond offsets, so ghosts cannot be reconstructed from saved state (`public/core/sessionStateManager.js`, `fairchem_local_server2/session_models.py`).
- **Docs Sync:** The Markdown guides (README, frontend/backend design, state system updates, timeline guide, testing catalog, VR setup) are consistent with the current code; they emphasize the dual-mesh renderer, timeline schema v6, and the requirement to keep logging/searchability via the provided debug flags. No conflicting instructions were found.

## Ghost Bond Regression (Step 3)

### Root-Cause Hypothesis
1. Running a simulation wraps atom positions into the primary cell (`moleculeState.markPositionsChanged`), which is correct for live rendering.
2. Ghost atoms still render because we clone every atom for {0, ±a, ±b, ±c}, but ghost bonds are recomputed by rerunning `computeBondsNoState` over that limited augmentation.
3. Bonds that require a diagonal lattice shift (e.g., `a+b`) vanish: the bond service still flags them as `crossing`, yet the ghost generator never sees a close pair because the necessary shifted copy was not created. After any timestep nudges the geometry, more bonds fall into this category, so only ghost atoms remain visible.

### Option A – Extend Bond Service With Periodic Offsets (Recommended)
- **Implementation**
  - Teach `computePeriodicBonds()` to retain `crossing`, `du/dv/dw`, and the integer image delta used for the minimum-image recompute.
  - Persist the additional fields in `molState.bonds` (`{ i, j, opacity, crossing, imageDelta: [Δu, Δv, Δw] }`).
  - Update `SessionStateManager` schema to include bond `imageDelta` so JSON snapshots round-trip the periodic intent.
  - Replace the ad-hoc augmentation in `rebuildGhosts()` with a helper that builds ghost bonds straight from the stored `imageDelta`, producing a matrix for each `(i, j)` plus its mirrored image without re-running `computeBondsNoState`.
- **Pros**
  - Single source of truth for minimum-image math (bond service) shared by renderer, tests, and snapshots.
  - Ghosts and primary bonds stay in lock-step, even after refactors to the bonding thresholds.
  - Simplifies debugging: logs already emitted by the bond service can be correlated with rendered ghost instances.
- **Cons**
  - Schema bump required (v6), fixtures/tests need updates.
  - Any consumer relying on `molState.bonds` shape must be adjusted (selection tests, snapshot fixtures).
  - Adds extra data to each bond entry (small memory increase).

### Option B – Recompute Ghost Bonds Locally With Fractional Wrapping
- **Implementation**
  - Embed the fractional/`wrap` logic from `bondService` directly inside `rebuildGhosts()`.
  - Expand the lattice augmentation to cover all ± combinations of {a,b,c} (26 offsets) before calling `computeBondsNoState`.
  - Cache the fractional transform to avoid per-frame matrix inversion churn.
- **Pros**
  - No schema change; the renderer remains self-contained.
  - Faster to land if we only need to un-break ghost bonds for release.
- **Cons**
  - Duplicates complex periodic math in two places, increasing maintenance cost.
  - Still heuristic: `computeBondsNoState` thresholds could drift from bond service limits, reintroducing divergence later.
  - Larger augmentation set (26 images) costs more CPU/GPU every time ghosts rebuild.

### Recommendation
Adopt **Option A**. It centralises periodic knowledge, keeps JSON snapshots authoritative, and reduces recomputation. We can still guard the renderer with a fallback (if `imageDelta` is missing, fall back to the existing augmentation) to preserve backward compatibility during the schema transition.

## Selection Status Bar Relocation
- Move the `#timelineEditorStatus` host from bottom-centre to top-centre by updating the injected CSS and ensuring it respects safe area insets.
- Add a responsive breakpoint (e.g., clamp width and adjust `gap`) so it avoids overlapping the VR toolbar.
- Provide a utility on the timeline editor module to re-measure the dock height when the status bar moves, keeping hover hitboxes aligned.
- Tests: update the Playwright authoring-mode suite to assert the bar’s bounding rect is near the top, and ensure existing status text still reflects selection changes.

## State Management Touchpoints
- **Frontend**
  - After Option A, `SessionStateManager.captureSnapshot()` must copy `bond.imageDelta` into `snapshot.viewer.meshAssignments` (or a dedicated `viewer.periodicBonds` array) and reapply it in `loadSnapshot`.
  - `stateStore` should expose a helper to reset `latchedUntil` when rehydrating periodic modes so ghost visibility stays in sync after a load.
- **Backend**
  - `SessionState` can optionally cache the last `imageDelta` per bond to validate client uploads or to inform future server-side ghost exports.
  - Add targeted logging (behind `MLIPVIEW_RESUME_DEBUG` or a new flag) enumerating how many periodic bonds were generated per frame; this will aid regression triage.

## JSON Snapshot Format (Schema v6 Draft)

```jsonc
{
  "schemaVersion": 6,
  "savedAt": "2025-01-11T02:15:43.210Z",
  "source": { "kind": "xyz", "label": "sn2-reactants" },
  "viewer": {
    "elements": ["C", "Cl", "H", "..."],
    "positions": [[-2.14, 0.03, 0.12], "..."],
    "velocities": [[0.0, 0.0, 0.0], "..."],
    "cell": [[6.5, 0, 0], [0, 6.5, 0], [0, 0, 6.5]],
    "showCell": true,
    "showGhostCells": true,
    "meshAssignments": {
      "atoms": ["solid", "solid", "soft", "..."],
      "bonds": ["solid", "soft", "..."]
    },
    "periodicBonds": [
      { "i": 12, "j": 37, "opacity": 0, "imageDelta": [1, 0, -1] },
      { "i": 5, "j": 9, "opacity": 0.18, "imageDelta": [0, 0, 0] }
    ]
  },
  "render": {
    "overrides": {
      "atoms": { "soft": [2, 5, 9] },
      "bonds": { "soft": [1, 4], "solid": [7] }
    }
  },
  "timeline": {
    "capacity": 500,
    "frames": [
      {
        "id": "frame-00495",
        "numericId": 495,
        "kind": "md",
        "seq": 812,
        "simStep": 138,
        "userInteractionCount": 41,
        "timestamp": 1736552663210,
        "energy": -152.43,
        "temperature": 305.4,
        "stress": [[...], [...], [...]],
        "energyIndex": 221,
        "positions": [[...]],
        "velocities": [[...]],
        "forces": [[...]]
      }
    ],
    "lastLiveMode": "md",
    "wasRunning": true,
    "pendingSimParams": { "temperature": 305.4, "friction": 0.02, "timestep_fs": 1 },
    "playback": {
      "autoPlay": false,
      "defaultFps": 20,
      "loop": false,
      "loopRange": null
    },
    "controlMessages": [
      {
        "id": "approach-phase",
        "label": "Approach",
        "actions": [
          { "type": "timeline.playbackSpeed", "fps": 12 },
          {
            "type": "visual.opacityFocus",
            "mode": { "focus": "solid", "background": "soft" },
            "focus": { "atoms": [0, 1, 2], "includeBonds": "connected" }
          }
        ]
      }
    ]
  },
  "energyPlot": {
    "series": [
      { "energy": -150.21, "kind": "idle" },
      { "energy": -152.43, "kind": "md" }
    ],
    "markerIndex": 221
  },
  "websocket": {
    "seq": 812,
    "nextSeq": 813,
    "clientAck": 809,
    "userInteractionCount": 41,
    "totalInteractionCount": 57,
    "simStep": 138
  }
}
```

This structure stays backward-compatible: older snapshots lacking `periodicBonds` will still load, and the renderer can regenerate ghost data via the fallback path.

## Testing & Instrumentation Plan
- **Unit**
  - New tests for `bondService.computePeriodicBonds()` verifying it emits `imageDelta` for crossing bonds and persists zero offsets for non-periodic ones.
  - `moleculeView` spec that feeds in periodic bonds with deltas and asserts ghost links survive a simulated timestep.
  - `SessionStateManager` round-trip covering schema v6 migration and legacy v5/v4 upgrades.
- **Playwright**
  - Extend `ws-ghost-periodic.spec.js` to run a single MD step and confirm `viewerApi.debugGhostSnapshot()` still reports ghost bonds.
  - Update authoring-mode suite to check the relocated status bar and verify selection text updates while at the top.
- **Backend**
  - Pytest ensuring the server accepts snapshots carrying `imageDelta` and ignores them safely when marshaling `ClientAction`.
  - Optional: assert debug logging for periodic bonds toggles correctly under `MLIPVIEW_RESUME_DEBUG`.
- **Logging**
  - When debugging, enable `window.__MLIPVIEW_DEBUG_STRETCH` and capture console output before/after a simulated step to compare periodic bond counts (documented in `testing.md`).

## Next Steps (Step 5 pending approval)
1. Land schema update & bond service refactor (Option A), guarded by feature flag if needed.
2. Update renderer, session manager, fixtures, and tests accordingly.
3. Adjust selection status bar styles and Playwright assertions.
4. Run targeted Jest/Playwright/Pytest suites; share logs for periodic bond counters to confirm parity.
5. Once validations pass, request approval to proceed with implementation and merge.
