# State System Update (Release V1)

## Findings Recap
- Timeline playback still renders atoms and bonds with the exact Babylon masters used in live mode. When a mask asks for translucency, we simply rewrite thin-instance colours (`public/render/moleculeView.js:1317-1476`), forcing every instance through the same mesh regardless of opacity intent.
- Clearing a mask resets all alphas to `1.0`, wiping out bond opacity that was computed by `createBondService()` (periodic crossing suppression, stretch fading). Ghost atoms/bonds avoid the bug because they already rely on separate meshes.
- Prior attempts to solve depth ordering through material flags / depth-prepass tuning created regressions and did not hold up under VR/Quest loads.

## Revised Strategy — Dual Mesh Pipeline
We will render every atom and bond through one of two masters:

| Category | Mesh | Intended State |
| --- | --- | --- |
| Primary opaque | `atom_<el>_solid` / `bond_<key>_solid` | Default “solid” instances (opaque timeline frames, live streaming) |
| Primary translucent | `atom_<el>_soft` / `bond_<key>_soft` | Timeline fades, stretched bonds, periodic crossings |
| Ghost (all) | reuse `*_soft` meshes | Ghost atoms/bonds always live in the translucent pool |

Key rules:
- **Single responsibility** – each mesh has a fixed transparency configuration (material alpha, depth pre-pass, render group). We never toggle transparency flags inside a frame.
- **Instance routing** – when timeline playback (or any overlay) asks for translucency, we migrate the affected instances from the solid master to the soft master. Returning to opaque pushes them back to solid.
- **Ghost parity** – ghost atoms/bonds move to the soft meshes permanently so their behaviour matches the new translucent pipeline.
- **Bond opacity** – baseline opacity from `createBondService()` decides which pool an instance should start in (e.g., periodic-crossing bonds spawn in soft). Timeline masks can still move them, but clearing a mask simply restores the original pool, preserving the periodic logic.

## Implementation Plan

### 1. Mesh Management (`moleculeView`)
- Introduce master registries for `solid` and `soft` variants per element/bond key (`public/render/moleculeView.js`).
- Solid masters retain current material settings; soft masters use preconfigured translucent materials (alpha ≈0.4 for atoms, derived from bond opacity for bonds) with `Material.MATERIAL_ALPHABLEND` and depth-prepass enabled.
- Add helpers:
  - `attachInstance(group, mode, index, matrix)` – writes the matrix and inserts the instance into the correct master array.
  - `migrateInstance(group, fromMode, toMode, index)` – swaps buffers when a timeline mask toggles visibility.

### 2. Baseline Build
- When `rebuildAtoms()` and `rebuildBonds()` run, classify each instance as `solid` or `soft`:
  - Atoms default to `solid`.
  - Bond opacity `< solidThreshold` (≈0.97) routes to `soft`; periodic zero-opacity bonds live exclusively in `soft`.
  - Ghost atoms/bonds are always created in `soft`.
- Maintain `instanceMeta` arrays so we know each instance’s current mesh and can flip it without recomputing geometry.

### 3. Timeline & Overlays
- Refactor `setOpacityMask` so it translates focus/background intents into mesh mode switches rather than mutating alpha in place:
  - Accepts `{ mode: 'solid' | 'soft' }` decisions per atom/bond index (when provided) and derives modes from opacity fallbacks otherwise.
  - Uses the migration helpers to move instances to the soft masters (for timeline fades) or back to solid.
- Update `applyControlActions` (`public/index.js:1948`) so `visual.opacityFocus` resolves to mesh choices:
  - Focus atoms/bonds → `solid`.
  - Background atoms/bonds → `soft`.
  - Optional numeric opacities still map to solid vs soft using thresholds (`>=0.99 → solid`, `<0.99 → soft`).
- When masks clear, consult stored baseline metadata to send instances back to their original pool.

### 4. Ghost Rendering
- Update ghost builders to reuse the soft masters (`public/render/moleculeView.js:694-870`).
- Ensure ghost bond creation marks instances as “soft baseline” so timeline clears won’t push them into the solid pool.

### 5. Session State & Snapshots
- Extend `SessionStateManager` (`public/core/sessionStateManager.js`) to persist which instances default to soft vs solid:
  - Add `viewer.meshAssignments` with arrays `atoms: 'solid' | 'soft'`, `bonds: 'solid' | 'soft'`.
  - Record the most recent timeline override set so playback restores the same distribution on load.
- Bump snapshot `schemaVersion` to **5**. Migrate v4 snapshots by defaulting all atoms/bonds to `solid` unless a bond stored `opacity < 0.99` (then seed `soft`).

### 6. Control Message Schema
- Keep `visual.opacityFocus` actions but allow an optional `mode` field to explicitly request `soft` or `solid`. If absent, derive from provided `focusOpacity/backgroundOpacity` thresholds.
- Example action:
  ```jsonc
  {
    "type": "visual.opacityFocus",
    "mode": {
      "focus": "solid",
      "background": "soft"
    },
    "focus": { "atoms": [0, 1, 2], "includeBonds": "connected" }
  }
  ```

### 7. Testing Plan
- **Unit**
  - `tests/moleculeView.meshModes.spec.js` – ensure build assigns proper masters and migration swaps buffers without leaks.
  - `tests/controlMessageEngine.spec.js` – cover `mode.focus/mode.background` defaults and threshold derivation.
- `tests/sessionStateManager.spec.js` – verify schema v6 round-trips mesh assignments and upgrades v4 snapshots.
- **Playwright**
  - `ws-timeline-mesh-mode.spec.js` – confirm timeline fades move instances between solid/soft pools and return correctly.
  - `ws-ghost-periodic.spec.js` – scrape `viewerApi.debugGhostSnapshot()` ensuring ghost bonds stay soft in both live and timeline modes.
- **Manual**
  - Rotate high-atom-count scenes in timeline playback to confirm draw-order stability.
  - Enable `window.__MLIP_DEBUG_STRETCH = true` (or `?bondStretchDebug=1`) during live-mode drags to confirm stretched bonds migrate to the soft mesh; disable afterwards to keep logs clean.
  - Quest smoke: ensure soft masters do not overwhelm fill-rate (watch frame timing HUD).

## Updated Snapshot Skeleton (schemaVersion 6)
```jsonc
{
  "schemaVersion": 6,
  "viewer": {
    "elements": ["C", "Cl", "N", "..."],
    "positions": [[-2.14, 0.03, 0.12], "..."],
    "bonds": {
      "atomIndexA": [12, 5, "..."],
      "atomIndexB": [37, 9, "..."],
      "cellOffsetA": [[0, 0, 0], [0, 0, 0], "..."],
      "cellOffsetB": [[1, 0, -1], [0, 0, 0], "..."],
      "opacity": [0.0, 0.18, "..."],
      "length": [1.41, 1.09, "..."],
      "flags": { "inRing": [1, 0, "..."], "crossing": [1, 0, "..."] }
    },
    "ghostAtoms": [
      { "atomIndex": 12, "shift": [1, 0, 0], "position": [4.36, 0.03, 0.12] }
    ],
    "ghostBondMeta": [
      {
        "base": { "i": 12, "j": 37 },
        "shiftA": [1, 0, 0],
        "shiftB": [0, 0, 0],
        "imageDelta": [1, 0, -1]
      }
    ],
    "meshAssignments": {
      "atoms": ["solid", "solid", "soft", "..."],
      "bonds": ["solid", "soft", "..."]
    },
    "showCell": true,
    "showGhostCells": true
  },
  "timeline": {
    "frames": ["..."],
    "controlMessages": [
      {
        "id": "approach-phase",
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
    "render": {
      "overrides": {
        "atoms": { "soft": [18, 19, 20] },
        "bonds": {
          "soft": [2, 3, 4],
          "solid": [7]
        }
      }
    },
  "websocket": { "nextSeq": 413, "lastAck": 409, "simStep": 98 }
}
```

Overrides store zero-based atom and bond indices. `solid` arrays are optional and appear only when defaults are overridden to the opaque mesh.

## Risks / Concerns
- **Buffer churn** – Moving thin instances between masters reallocates buffers. We should batch migrations (collect indices first, update buffers once per group) to avoid per-atom thrash.
- **Memory footprint** – Doubling masters slightly increases GPU memory usage. Need to monitor on Quest; if necessary we can lazy-create soft masters only when a translucent request is made.
- **Authoring UX** – Editors may want finer granularity (e.g., multiple translucency levels). We can extend with more than two modes later, but for V1 the binary split keeps rendering deterministic.
- **Snapshot size** – Storing per-instance mesh assignments adds arrays proportional to atom/bond count. Compression is optional; we can delta-encode (e.g., store indices of soft instances only) if file size becomes an issue.

## Next Steps
1. OK the dual-mesh plan (adjust thresholds or schema details if needed).
2. Convert `public/examples/sn2/sn2.json` and other fixtures to schema v6 once implementation lands.
3. Proceed with coding and validation after you sign off.
