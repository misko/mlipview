# MLIP Viewer

Interactive molecular viewer with desktop and WebXR (VR) support.

## New: Molecule Selection & URL Parameter

The app now defaults to loading `roy.xyz` (ROY molecule). You can switch molecules dynamically using the Molecules button (desktop HUD or VR HUD) or by specifying a URL parameter.

### URL Parameter

```
?molecule=roy
?molecule=benzene
?molecule=someOther
```

The value is case-insensitive and you can omit the `.xyz` extension. The app searches for a matching `public/molecules/<name>.xyz` file.

### Fallback Order
1. Requested molecule (if provided via `?molecule=`)
2. `roy.xyz`
3. `benzene.xyz`
4. Procedural benzene (hardcoded) if all file loads fail

### Listing Available Molecules
The server now exposes an endpoint:

```
GET /api/molecules
```
Returns JSON:
```json
{ "molecules": [ { "file": "roy.xyz", "name": "roy" }, { "file": "benzene.xyz", "name": "benzene" } ] }
```

### UI Selector
Click the `Molecules` button to open an overlay listing all `.xyz` files in `public/molecules/`. Selecting one reloads the page with the `?molecule=` parameter applied.

Works in both:
- Desktop (`index.html`)
- VR (`vr.html`) via a VR HUD button (renders an HTML overlay when pressed)

### In-Headset (XR) Molecule Picker
Inside an active WebXR session, a dedicated 3D panel can be opened using the `Molecules` HUD button that appears in the VR overlay. Selecting a molecule from this panel immediately reloads the page with the new `?molecule=` parameter so the chosen structure loads after re-entry.

### Development Notes
- Selector UI defined in `public/ui/molecule-selector.js`
- Desktop HUD button injected via `hud.js` / wired in `main.js`
- VR button created in `vr-ui.js` and wired in `vr/main-vr.js`
- Molecule load logic refactored in `public/molecules/molecule-loader.js`
- Server endpoint added in `server.js`

### Adding New Molecules
Drop additional `.xyz` files into `public/molecules/` and they will automatically appear in the selector (no restart required unless server-side caching is added later).

---
For VR setup details see `VR_SETUP_README.md`.

### Spherical Atom Drag (VR/AR)
Advanced spherical drag with multiple radial push/pull modes (adaptive default) is documented in `docs/VR_SphericalDrag.md`. Runtime tuning via `window.vrSpherical*` flags enables rapid experimentation (gain, smoothing, working radius cap, mode switching). See that doc for formulas, defaults, and troubleshooting.

## Force Field / MLIP Interface

The app now exposes a unified force field interface supporting both the local Lennard-Jones mock potential and remote FAIR-Chem models on HuggingFace.

### Selecting Backend

Use a URL parameter:

```
?ff=lj          # (default) local LJ + harmonic bonds (scene-cached)
?ff=fairchem    # remote FAIR-Chem (requires network & HF token)
```

If creation fails (e.g. invalid kind or network error), the viewer falls back to `lj`.

### Interface Contract

Factory returns an object implementing:

```
{
	kind: string,             // 'lj' | 'fairchem' | ...
	units: { energy, length },
	meta: { ... backend info ... },
	compute(): { energy:number, forces: BABYLON.Vector3[] },
	computeRaw({ Z:int[], xyz:number[][] }): Promise<{ energy:number, forces:number[][] }>
}
```

Helper utilities & registry live in `public/forcefield/`.

### FAIR-Chem Configuration

File: `public/forcefield/fairchem.js`

By default it targets:

```
https://api-inference.huggingface.co/models/facebook/fairchem_mace_large
```

Provide a HuggingFace token at runtime (before initiating compute):

```
window.HF_TOKEN = 'hf_xxxYourToken';
```

or fork `fairchem.js` to hard-code one (not recommended). You may also pass a different endpoint via query parameters in future (simple extension: parse `?hf_endpoint=` and forward to factory).

#### Temporary Embedded Token (Development Only)

There is currently a temporary token embedded in `public/forcefield/fairchem.js` (`__EMBEDDED_HF_TOKEN`). This is ONLY for local experimentation. You MUST:

1. Replace it with `window.HF_TOKEN = 'hf_new_token';` or environment proxy before sharing.
2. Rotate / revoke the embedded token on HuggingFace once you finish testing.
3. Avoid committing real secrets—use a server proxy for production.

Response mapping defaults to fields `{ total_energy, forces }`. Adjust via `responseMap` if your Space outputs different keys.

### Extending

Add a new backend by editing `public/forcefield/registry.js`:

```
FACTORIES.mybackend = (opts) => createMyBackendForceField(opts);
```

Then load with `?ff=mybackend`.

### Raw vs Scene Compute

`compute()` returns Babylon `Vector3` forces for immediate visualization.
`computeRaw({Z, xyz})` is backend-neutral numeric data (no Babylon dependency) and can be used for batch or headless analysis.

## Local UMA FAIR-Chem Server (Simplified)

When developing locally over `http://` the FAIR-Chem adapter talks to `http://localhost:8000` (default UMA). If you access the viewer itself over HTTPS, use an external reverse proxy or tunnel (e.g. Caddy, Traefik, nginx, ngrok) to terminate TLS and forward to the Node server + UMA backend. The previous dual HTTP+HTTPS helper script has been removed to reduce maintenance.

Mixed-content avoidance is now handled by:
1. Same-origin proxy endpoints in `server.js` (`/simple_calculate`, `/calculate`).
2. Automatic upgrade attempt inside `public/forcefield/fairchem.js` (logs a single warning on fallback).

If you need HTTPS for VR quickly:
```
caddy reverse-proxy --from https://localhost:8443 --to :3000
```
(Accept the local trust prompt / add a dev certificate.)

Production-grade deployment should still use a hardened reverse proxy with proper certificates.


## Test Pyramid Strategy

The test suite is layered for speed and fidelity:

1. Core (Jest / Node) – logic & state (`tests`, `mlipviewer2/tests`). Uses lightweight Babylon mocks.
2. Browser (Jest / jsdom) – DOM integration / highlight & canvas behavior (`tests-browser`). Adds minimal canvas + pointer event polyfills.
3. Engine (NullEngine harness) – optional Babylon scene graph checks without real GL (`tests/engine`).
4. Future (Playwright) – real headless browser for WebGL rendering regression (bonds present, no stray primitives). Not yet added, but structure is prepared.

Run all:
```
npm test
```
Only browser project:
```
npx jest --selectProjects browser
```
Only core project:
```
npx jest --selectProjects core
```
Engine-only sample:
```
npx jest tests/engine/nullEngine.spec.js --selectProjects core
```

Adding a browser test: drop a spec in `tests-browser/*.spec.js` and it will run under `jsdom` automatically.

## Static Build (Vite)

You can produce a fully static production bundle (HTML + ES modules + assets) using Vite.

Scripts added:

```
npm run build    # Generates dist/ using Vite (root = public/)
npm run preview  # Serves the built dist/ locally on port 5174
```

Build output goes to `dist/` and can be deployed behind any static file host or CDN (Netlify, GitHub Pages, nginx, S3 + CloudFront, etc.). Requires Node 20.19+ (see `.nvmrc`).

Notes:
1. External Babylon CDN scripts referenced in `public/index.html` are preserved as-is (they are not bundled). This keeps the bundle lean and leverages browser caching. For full offline packaging you can vendor these scripts by removing the `<script src=...>` tags and importing modules instead.
2. Environment-specific settings (e.g. remote FAIR-Chem endpoints) should be passed via query params or a small runtime config file; no `.env` integration has been added yet.
3. If you only need a static viewer (no Node server API routes), you can deploy just the contents of `dist/`. The express server remains useful for local dev and any dynamic endpoints (molecule listing, force proxy, etc.).
4. Source maps are enabled for easier debugging of production issues (`dist/**/*.js.map`).

Minimal end-to-end check:

```
npm install
npm run build
npm run preview
# open http://localhost:5174
```

If you see the viewer and can load the default molecule, the static build is functioning.

