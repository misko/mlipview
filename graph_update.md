# Graph & Bond Refactor Proposal

## Current Gaps
- **Cross-cell bonds are filtered out.** The augmentation already considers the primary cell plus the six axis-aligned neighbours, but we drop any candidate when either endpoint lives in a shifted cell because `bondMap` only accepts `isPrimary` pairs ([public/domain/bonding/periodicAugment.js:254](public/domain/bonding/periodicAugment.js:254)). As a result, bonds that should connect atoms across the `±x`, `±y`, or `±z` images never reach rendering.
- **State surfaces don’t expose canonical offsets.** `createBondService()` stores `imageDelta` without the explicit `(atomIndex, cellOffset)` tuples ([public/domain/bondService.js:67](public/domain/bondService.js:67)), so `rebuildGhosts()` reconstructs transforms indirectly ([public/render/moleculeView.js:1090](public/render/moleculeView.js:1090)) and other consumers cannot reason about which cell an atom belongs to.
- **State sync fragmentation.** Snapshot capture/migration lives in `SessionStateManager` ([public/core/sessionStateManager.ts:4](public/core/sessionStateManager.ts:4)), yet backend `SessionState` still mirrors only flat geometry ([fairchem_local_server2/session_models.py:63](fairchem_local_server2/session_models.py:63)), so periodic graph metadata is lost when sessions round-trip through JSON.

## Proposed `compute_bonds` Core
We will replace the current augmentation helper with a pure function that matches the requested signature:

```ts
type Vec3 = readonly [number, number, number];
type Cell = readonly [Vec3, Vec3, Vec3] | null;

export function compute_bonds(positions: Vec3[], cell: Cell | null): {
  atom_idxA: Int32Array;
  atom_idxB: Int32Array;
  cell_offsetsA: Int16Array; // flattened triplets
  cell_offsetsB: Int16Array;
  metrics: {
    length: Float32Array;
    weight: Float32Array;
    opacity: Float32Array;
    inRing: Uint8Array;
    crossing: Uint8Array;
  };
};
```

Implementation sketch (periodic):
1. Recreate the existing augmentation but capture metadata: generate seven images (primary + `±x`, `±y`, `±z`) and record the shift applied to each atom. No diagonal images are produced.
2. Run the base neighbour search (`computeBaseBonds`) over the augmented list. For each candidate, read the base indices and their corresponding shift vectors (`shiftA`, `shiftB`).
3. Canonicalise by comparing `(atomIndex, shift)` tuples; swap ends when needed so every bond is emitted once even if both endpoints live outside the primary cell.
4. Populate flat arrays for indices, offsets, and metrics (length, weight, opacity, flags) using the canonical tuples. Deduplicate using a map keyed by `(i, shiftA, j, shiftB)`.
5. Emit `ghostBondMeta` using the canonical offsets so the renderer can seed ghost thin instances without recomputing augmentation.

Non-periodic callers skip steps 1–4, derive weights directly from Euclidean distance, and set both offsets to `(0,0,0)`.

### Options
1. **Option A – Canonical axis augmentation (recommended)**
   - *Pros:* Respects the existing seven-image policy; fixes the filtering bug by admitting cross-cell pairs; produces canonical tuples reusable by all consumers.
   - *Cons:* Requires careful canonical ordering to avoid duplicate entries and explicit guards that reject any offset outside `{-1,0,1}`.
2. **Option B – Minimal patch**
   - Keep the current augmentation but simply drop the `isPrimary` guard and add tuple storage for the renderer.
   - *Pros:* Lowest implementation risk; quick fallback if timelines demand.
   - *Cons:* Leaves canonical ordering implicit and makes later refactors harder because other systems still need to reconstruct offsets on their own.

## Graph & State Refactors
- **Shared graph module.** Extract weight/ring helpers into `public/domain/bonding/graphCore.js` so both periodic and non-periodic computations reuse the same primitives. Export utilities to build adjacency lists (needed for selection and analytics downstream).
- **Front-end state integration.**
  - Store the flattened offset arrays (`bondOffsetsA/B`) on `molState` alongside the high-level bond list so `moleculeView` can feed them directly into ghost builders without reconstructing shifts.
  - Extend `SessionStateManager.captureSnapshot()` to persist the new arrays under `viewer.bonds`, replacing the existing `periodicBonds` structure. Migrations v6→v7 will synthesise offsets from legacy `imageDelta` when possible.
- **Backend parity.**
  - Mirror `bond_offsets_a/b` fields on `SessionState` so saved sessions survive a round-trip through `/serve` (add optional validation when dense snapshots arrive).
  - If we ever ship backend-rendered graphs, the same data is ready for protobuf extensions.
- **Session manager consolidation.** Keep interaction counters and drag gating under the session manager so JSON snapshots remain self-contained, avoiding hidden singletons that might skip resets during restore.

## Rendering & Ghosts
- With per-end offsets available, `rebuildGhosts()` can compute ghost transforms directly from the canonical tuples instead of walking `ghostBondMeta`. That reduces duplication and ensures ghost meshes stay in sync with whatever `compute_bonds` emits.
- Keep `ghostBondMeta` as a derived artefact for debugging/tests (fill it from the canonical arrays during rebuild). This maintains compatibility with Playwright helpers that call `viewerApi.debugGhostSnapshot()`.

## JSON Snapshot (schema v7 draft)
```jsonc
{
  "schemaVersion": 7,
  "viewer": {
    "elements": ["C", "H", "..."],
    "positions": [[-2.14, 0.03, 0.12], "..."],
    "cell": [[6.5, 0, 0], [0, 6.5, 0], [0, 0, 6.5]],
    "showCell": true,
    "showGhostCells": true,
    "bonds": {
      "atomIndexA": [12, 5, "..."],
      "atomIndexB": [37, 9, "..."],
      "cellOffsetA": [[0, 0, 0], [0, 0, 0], "..."],
      "cellOffsetB": [[1, 0, 0], [0, 0, 0], "..."],
      "opacity": [0.0, 0.18, "..."],
      "length": [1.41, 1.09, "..."],
      "flags": { "inRing": [1, 0, "..."], "crossing": [1, 0, "..."] }
    },
    "meshAssignments": { "atoms": ["solid", "..."], "bonds": ["soft", "..."] }
  },
  "timeline": {
    "frames": [...],
    "controlMessages": [...],
    "playback": { "defaultFps": 20, "autoPlay": false }
  },
  "websocket": {
    "nextSeq": 813,
    "clientAck": 809,
    "userInteractionCount": 41,
    "simStep": 138
  }
}
```

Migration strategy:
- **v5/v6 → v7:** promote `viewer.periodicBonds` into the new arrays by assigning `cellOffsetA=(0,0,0)` and `cellOffsetB=imageDelta`. Any missing opacity defaults to `1.0`.
- **v7 → runtime:** `SessionStateManager` hydrates `molState` bonds/offsets directly so `rebuildGhosts()` can update thin-instance buffers without recomputing periodic images.

## Testing Matrix
- **Unit**
  - `tests/bonding/compute_bonds.spec.js`: cover non-periodic cases plus each axis shift (`±x`, `±y`, `±z`) where one or both atoms sit in non-zero cells. Assert canonical tuples and enforce the `{ -1, 0, 1 }` constraint.
  - `tests/bonding/ringMetrics.spec.js`: verify extracted helpers still classify aromatic rings identically to the previous pipeline.
  - `tests/sessionStateManager.spec.js`: migrate v6 snapshots, ensure offsets survive capture/restore, and confirm counters move with the rest of the state.
- **Integration**
  - `tests/moleculeView.ghostBonds.spec.js`: feed synthetic bond data with canonical axis offsets, assert ghost meshes render whenever a bond references a non-zero cell, and guard the regression where `rebuildBonds()` must refresh ghosts without an extra viewer call.
  - `tests/moleculeView.periodicToggle.spec.js`: place two carbons near one another inside a large box; verify one opaque bond appears when periodic is off, then assert periodic mode yields the same primary bond plus six ghost bonds.
  - `tests/moleculeView.largeCellEdge.spec.js`: use a 10×10×10 cell with carbons at `(-4.5, 0, 0)` and `(4.5, 0, 0)`; confirm no bonds render when periodic is disabled and exactly two ghost bonds appear (no primary bond) once periodic is enabled.
  - `tests/apiCellPayload.spec.js`: extend parity checks so backend snapshots accept the new bond structure.
- **Playwright**
  - Update `ws-ghost-periodic.spec.js` with a scenario where a bond crosses into a `+x` (or other axis) cell and verify ghost bonds persist through timeline playback and live resume.
  - Add `ws-bond-offsets-json.spec.js`: save/load a session containing periodic offsets and confirm the viewer resumes continuous MD without losing bonds.
- **Manual/Logging**
  - Enable `window.__MLIP_DEBUG_STRETCH` and capture the new per-offset summaries during MD, verifying the canonical tuples match expectations.

## Next Steps
1. Finalise the API (`compute_bonds`), factoring shared math utilities and documenting canonical ordering.
2. Update `createBondService`, `moleculeView`, and `SessionStateManager` to consume the new payload and populate the snapshot schema.
3. Mirror the structure on the backend (`SessionState`, optional protobuf extension for future streaming).
4. Migrate fixtures / JSON sessions to schema v7 and refresh affected Playwright baselines.
5. Land the test suite updates, capture diagnostic logs for periodic systems, and review performance impact before requesting implementation approval.
