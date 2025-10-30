# Bond Computation Refactor Plan

## Current Pain Points
- **Duplicate logic**: Periodic handling lives in `bondService.js` while the base detection heuristics (covalent radii, ring boosts) live in `bond_render.js`. Ghost reconstruction re-implements parts of the periodic math in `moleculeView.js`.
- **Stateful side effects**: `bondService.recomputeAndStore()` mutates `molState` directly, making it difficult to reuse outside the viewer or to unit-test in isolation.
- **Recompute pressure**: The frontend calls `recomputeAndStore()` on every position change, even if topology cannot change (e.g., minor drags) and even for components that only need the derived metadata.
- **Inconsistent outputs**: The renderer requires `{ crossing, imageDelta }` to keep ghost bonds alive, but the pure helper (`computeBondsNoState`) does not emit that metadata.
- **Testing gaps**: No dedicated unit/benchmark suites exist for periodic vs non-periodic inputs; coverage relies on viewer integration tests that are expensive to run.

## Target: Single Modular Entry Point
Create a pure utility `computeBonds()` that accepts atom positions, element symbols, and optional cell metadata, and returns a normalized `BondResult` payload. All consumers (viewer, tests, future CLI exporters) use this function.

### Proposed API
```ts
type Vec3 = readonly [number, number, number];

interface Cell {
  enabled: boolean;
  a: Vec3;
  b: Vec3;
  c: Vec3;
  origin?: Vec3;
}

interface BondResult {
  bonds: Array<{
    i: number;
    j: number;
    length: number;
    weight: number;
    opacity: number;
    inRing: boolean;
    crossing: boolean;
    imageDelta: [number, number, number];
  }>;
  ghostImages?: Array<{
    atomIndex: number;
    shift: [number, number, number];
    position: Vec3;
  }>;
  diagnostics?: {
    neighborCounts: number[];
    prunedPairs: number;
    periodicAugmentations: number;
  };
}

export function computeBonds(input: {
  elements: readonly string[];
  positions: readonly Vec3[];
  cell?: Cell | null;
  options?: { maxNeighbors?: number; opacityGamma?: number };
}): BondResult;
```

### Module Layout
```
public/domain/bonding/
  ├─ covalent.ts           # radii, heuristics (ported from bond_render)
  ├─ neighborSearch.ts     # shared neighbor pruning
  ├─ ringMetrics.ts        # aromatic/planarity scoring
  ├─ periodicAugment.ts    # cell wrapping & image delta calculation
  ├─ computeBonds.ts       # orchestrates the full pipeline (exported API)
  └─ index.ts              # re-export convenience
```

## Implementation Steps
1. **Extract base detection**
   - Move the current `computeBondsNoState` logic into `neighborSearch.ts` + `ringMetrics.ts`.
   - Return raw edges (`{ i, j, distance, weight, metadata }`) without opacity/crossing applied.
2. **Implement periodic augmentation**
   - Factor the matrix inversion + fractional wrapping (`bondService`) into `periodicAugment.ts`.
   - Produce both primary bonds and periodic metadata in a single pass (no double traversal).
   - Output ghost atom seeds for the renderer so it no longer rebuilds them ad hoc.
3. **Compose in `computeBonds`**
   - Merge primary and periodic paths; apply opacity/crossing via a dedicated formatter.
   - Accept typed arrays (e.g., `Float32Array`) seamlessly by normalizing input upfront.
   - Support optional tuning (max neighbors, distance cutoff overrides) for future CLI tooling.
4. **Update consumers**
   - `createBondService` becomes a thin wrapper that calls `computeBonds` and stores results on `molState`.
   - `moleculeView` reads `ghostImages` from the bond result instead of recomputing periodic images.
   - `SessionStateManager` snapshot logic stores the `BondResult` metadata (bond list + ghost deltas).
   - Clean up legacy helpers (`computeBondsNoState`), replacing usage with the new module.
5. **Performance instrumentation (stateless recompute)**
   - Recompute bonds every frame to keep rendering independent of prior results.
   - Offer optional accuracy knobs (e.g., skip ring detection) that callers can toggle explicitly when they accept approximate visuals.
   - Add lightweight timing hooks (dev-only) so automated tests can record runtimes without affecting production builds.
6. **Documentation & migration**
   - Document the new API in `frontend_design.md` and note the schema changes (JSON snapshots now store full `BondResult` metadata).
   - Provide migration path for tests referencing old fields (`bonds`, `crossing`, etc.).

## Testing Plan
### Unit
- `tests/bonding/neighborSearch.spec.js` – verify max-neighbor pruning and distance thresholds.
- `tests/bonding/periodicAugment.spec.js` – check `imageDelta` outputs for cubic, triclinic, and skew cells.
- `tests/bonding/computeBonds.spec.js` – snapshot `BondResult` for representative molecules (benzene, SN2, water box) with and without periodic cells.
- `tests/bonding/performance.spec.js` – benchmark systems of different sizes (e.g., N≈256, N≈1024); assert average runtime stays under agreed thresholds and emit timing logs for trend tracking.

### Integration
- Update existing viewer tests (`x-bondOpacityIsolation`, ghost-related DOM specs) to assert the new ghost data comes from the API (no renderer recompute).
- Add a new Playwright check that toggles periodic cells mid-simulation and ensures ghost bonds persist through a single-step run.

### Backend Compatibility
- Ensure the backend `SessionState` accepts the new `BondResult` payload (even if unused) to keep JSON snapshots symmetric between client/server.

## Cleanup & Follow-up
- Remove legacy globals (`window.__MLIP_DEBUG_STRETCH`) in favor of structured diagnostics returned by `computeBonds`.
- Replace scattered numeric literals (opacity thresholds, distance caps) with shared constants exported from `covalent.ts`.
- Purge the duplicated ghost reconstruction logic from `moleculeView`, using the new `ghostImages` array directly.
- Consider exposing `computeBonds` through a small CLI tool (e.g., `npm run bonds:inspect`) for regression triage.

## Review Checklist
- [ ] API signature and module structure approved.
- [ ] Back-compat strategy for JSON snapshots acceptable.
- [ ] Testing scope sufficient (unit + integration + perf).
- [ ] Recompute/timing strategy agreed (no caching, perf instrumentation in place).
- [ ] Rollout notes captured (docs/tests updated).
