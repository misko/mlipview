# JSON Library Format (Schema v6/v7 Overview)

This document captures the full structure of timeline/library exports that the MLIP viewer reads and writes.  The format is JSON with a single top‑level object.  Unless otherwise noted, numeric arrays are in **Ångström** units and use right‑handed world coordinates (Babylon’s X, Y, Z).

## 1. Top-Level Envelope

```jsonc
{
  "schemaVersion": 6,
  "savedAt": "2025-10-30T14:07:17.779Z",
  "source": { "kind": "xyz-inline", "label": "inline.xyz" },
  "viewer": { ... },
  "render": { ... },
  "timeline": { ... },
  "energyPlot": { ... },
  "websocket": { ... }
}
```

| Field | Type | Notes |
| --- | --- | --- |
| `schemaVersion` | integer | Current writer is **6**.  Future versions bump the number; a schema ≥7 uses the same layout but adds canonical bond offsets (`atomIndexA/B`, `cellOffsetA/B`). |
| `savedAt` | ISO-8601 string | Timestamp when the snapshot was captured. |
| `source` | object | Optional provenance.  `kind` is a short identifier (`"xyz"`, `"xyz-inline"`, `"session"`); `label` is human-readable. |

## 2. `viewer` Block

Contains the instantaneous rendering state.

```jsonc
"viewer": {
  "elements": ["C", "H", "Cl"],
  "positions": [[2.28, 0.69, 1.03], "..."],
  "velocities": [[0.0, 0.0, 0.0], "..."],
  "cell": [[6.5, 0, 0], [0, 6.5, 0], [0, 0, 6.5]],
  "showCell": true,
  "showGhostCells": true,
  "bonds": {
    "atomIndexA": [12, 5],
    "atomIndexB": [37, 9],
    "cellOffsetA": [[0, 0, 0], [0, 0, 0]],
    "cellOffsetB": [[1, 0, 0], [0, 0, 0]],
    "opacity": [0.0, 0.18],
    "length": [1.41, 1.09],
    "flags": {
      "inRing": [1, 0],
      "crossing": [1, 0]
    }
  },
  "ghostAtoms": [
    { "atomIndex": 12, "shift": [1, 0, 0], "position": [4.36, 0.03, 0.12] }
  ],
  "ghostBondMeta": [
    {
      "base": { "i": 12, "j": 37 },
      "shiftA": [1, 0, 0],
      "shiftB": [0, 0, 0],
      "imageDelta": [-1, 0, 0]
    }
  ],
  "meshAssignments": {
    "atoms": ["solid", "soft", "..."],
    "bonds": ["solid", "soft", "..."]
  }
}
```

### Required arrays
- `elements`: element symbols for each atom.
- `positions`: list of `[x, y, z]` real-space coordinates in Å.  Indices align with `elements`.
- `cell`: optional 3×3 cell vectors.  If `null` or omitted, the system is non-periodic.  When present, `showCell` toggles whether it is drawn.

### Bond payload (`schemaVersion >= 6`)
`bonds` is the canonical, deduplicated list.  Each index `k` describes bond `(atomIndexA[k], atomIndexB[k])`.  Offsets and flags use parallel arrays:

| Field | Meaning |
| --- | --- |
| `cellOffsetA/B` | Integer cell coordinates (Δu, Δv, Δw). A primary-cell bond uses `[0,0,0]`; a periodic crossing stores the lattice image for each endpoint. |
| `opacity` | Default opacity used to place the bond into “solid” or “soft” render buckets. Crossing bonds normally have `0`. |
| `flags.inRing` | `1` when the bond participates in an aromatic/small ring. |
| `flags.crossing` | `1` when the bond crosses the periodic boundary. |
| `length` | Cached bond length in Å. |
| `weight` | Optional numeric weight used by the bond service (typically `null` unless diagnostics are enabled). |
| `imageDelta` | Integer lattice delta (`shiftB - shiftA`) retained for diagnostics and future tooling. |

Snapshots must already use schema v6; earlier schema versions are rejected during load.

### Ghost metadata
- `ghostAtoms`: instanced clones generated for periodic visualization.  `shift` is the integer lattice translation; `position` is the already-shifted world coordinate used for rendering.
- `ghostBondMeta`: describes thin-instance ghost bonds.  `base` references the underlying atom indices; `shiftA/B` are per-end translations; `imageDelta` equals `shiftB - shiftA`.

### Mesh assignments
`meshAssignments` stores which render bucket each atom or bond currently occupies (`"solid"` or `"soft"`).  When absent, everything defaults to `"solid"`. For narrated sequences, keep the baseline snapshot fully `"solid"` and rely on `visual.opacityFocus` actions to introduce translucency. Persisting `"soft"` overrides in the JSON forces the opening frames to appear semi-transparent before the control message triggers.

## 3. `render` Overrides

Optional user overrides (e.g. opacity masks) are collected under `render`.

```jsonc
"render": {
  "overrides": {
    "atoms": { "soft": [2, 5, 9] },
    "bonds": { "soft": [1, 4], "solid": [7] }
  }
}
```

## 4. Timeline block

Snapshots of the frame buffer, control messages, and playback settings.

```jsonc
"timeline": {
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
      "stress": [[ ... ], [ ... ], [ ... ]],
      "energyIndex": 221,
      "positions": [[ ... ]],
      "velocities": [[ ... ]],
      "forces": [[ ... ]]
    }
  ],
  "playback": {
    "defaultFps": 20,
    "autoPlay": false,
    "loop": false,
    "loopRange": null
  },
  "controlMessages": [
    {
      "id": "approach-phase",
      "priority": 10,
      "range": {
        "start": { "frameId": "frame-00480" },
        "end":   { "frameId": "frame-00510", "inclusive": true }
      },
      "actions": [
        { "type": "timeline.playbackSpeed", "fps": 12 },
        {
          "type": "overlay.callout",
          "text": "Approach the chloride ion",
          "anchor": { "mode": "atom", "atoms": [12] },
          "offset": [0.0, 2.0, 0.0],
          "panelSize": { "width": 1.6, "height": 0.8 },
          "textSize": 0.9,
          "style": { "background": "rgba(12,16,24,0.85)", "color": "#f6f8ff" }
        }
      ]
    }
  ]
}
```

### Frame entries
Each `frames[]` element lodges the full molecular state for a specific history offset.  Fields align with live WebSocket payloads (`seq`, `simStep`, `positions`, `forces`, etc.).

### Playback settings
- `defaultFps`: timeline playback rate when no override is active.
- `autoPlay`: enter timeline mode and start playback immediately on load.
- `loop` / `loopRange`: loop playback between specific frames.
- `lastLiveMode`: remembers which continuous mode (`"idle"`, `"md"`, or `"relax"`) the viewer was in before entering timeline mode so that resuming live view can restart the same mode automatically.
- `wasRunning`: boolean snapshot indicating whether the continuous loop was actively running when timeline mode was entered (`true` means resume the run on exit; `false` means stay paused).

### Control messages
- `id`, `priority`, `notes` (optional) identify the message.
- `range.start` / `range.end` specify when it’s active.  Each reference may use `frameId`, `frameIndex`, or `offset`.  `inclusive: false` on the end reference excludes the final frame.
- `actions`: ordered list; categories below.

#### 4.1 `timeline.playbackSpeed`
```jsonc
{ "type": "timeline.playbackSpeed", "fps": 12 }
```
Alternative: `{"type":"timeline.playbackSpeed","speedMultiplier":0.5,"transitionMs":500,"easing":"linear"}`.

#### 4.2 `visual.opacityFocus`
Applies opacity masks by moving atoms/bonds between render buckets.
```jsonc
{
  "type": "visual.opacityFocus",
  "focus": { "atoms": [5, 6], "includeBonds": "connected" },
  "focusOpacity": { "atoms": 1.0, "bonds": 1.0 },
  "backgroundOpacity": { "atoms": 0.2, "bonds": 0.05 },
  "mode": { "focus": "solid", "background": "soft" },
  "transitionMs": 300
}
```

#### 4.3 `overlay.callout`
Displays a billboarded text panel anchored to an atom, bond, or world coordinate.

```jsonc
{
  "type": "overlay.callout",
  "text": "Labile chloride",
  "anchor": { "mode": "atom", "atoms": [12] },
  "offset": [0.0, 2.0, 0.0],
  "panelSize": { "width": 1.4, "height": 0.7 },
  "textSize": 0.9,
  "style": {
    "background": "rgba(15,22,34,0.85)",
    "color": "#f4f9ff",
    "cornerRadius": 12,
    "opacity": 0.95
  }
}
```

**Anchor modes**

| Mode | Fields | Behaviour |
| --- | --- | --- |
| `"atom"` | `atoms: [index]` | Anchors to the chosen atom in world coordinates. |
| `"bond"` | `atoms: [i, j]` or `i`, `j` | Anchors to the midpoint of the two atoms. |
| `"world"` | `position: [x, y, z]` | Anchors to an absolute world-space position. |
| `"xyz"` (legacy) | `xyz` or `value` array | Treated like `"world"`. |

**Offset semantics**

- `offset` is optional `[dx, dy, dz]`.  When specified, the final position is `anchorPoint + offset`.
- Axes align with the Babylon scene: X → right, Y → up, Z → toward/away from the camera.  Changing **Y** moves the callout vertically; varying **Z** changes depth (often subtle because the plane billboards toward the camera).
- If `offset` is omitted, the raw anchor position is used.

## 5. `energyPlot`

Stored energy history (used for the timeline plot marker).

```jsonc
"energyPlot": {
  "series": [
    { "energy": -150.21, "kind": "idle" },
    { "energy": -152.43, "kind": "md" }
  ],
  "markerIndex": 221
}
```

## 6. `websocket` State

Tracks counters/cursors needed to resume streaming without desynchronising ACKs.

```jsonc
"websocket": {
  "seq": 812,
  "nextSeq": 813,
  "clientAck": 809,
  "userInteractionCount": 41,
  "totalInteractionCount": 57,
  "simStep": 138
}
```

## 7. Minimal Example

```jsonc
{
  "schemaVersion": 6,
  "savedAt": "2025-01-01T00:00:00Z",
  "viewer": {
    "elements": ["H", "H"],
    "positions": [[0, 0, 0], [0.74, 0, 0]],
    "showCell": false,
    "showGhostCells": false,
    "bonds": {
      "atomIndexA": [0],
      "atomIndexB": [1],
      "cellOffsetA": [[0, 0, 0]],
      "cellOffsetB": [[0, 0, 0]],
      "opacity": [1.0],
      "length": [0.74],
      "flags": { "inRing": [0], "crossing": [0] }
    },
    "ghostAtoms": [],
    "ghostBondMeta": []
  },
  "timeline": {
    "frames": [],
    "playback": { "defaultFps": 20, "autoPlay": false, "loop": false, "loopRange": null },
    "controlMessages": []
  },
  "energyPlot": { "series": [], "markerIndex": null },
  "websocket": {}
}
```

> **Tip**: When authoring new library entries, start with the minimal example and progressively add sections.  Always update `schemaVersion` to the highest structure the viewer supports (6+ for canonical bond offsets).  Remember to call `Save Changes` in the timeline editor so UI adjustments (including callout offsets) propagate into the JSON.

## 8. Loading a Library Entry from the URL

The viewer can jump straight into a library session without manual UI interaction by appending `?library=<id>` to the page URL. For example:

```
http://localhost:5174/index.html?library=sn2
```

During startup the viewer looks up the requested `id` (using `library.json`), streams the corresponding JSON snapshot, and hydrates timeline/control state automatically. If the download fails, the viewer falls back to the default molecule and surfaces an error banner so narration sessions never leave the page blank.
