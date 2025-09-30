/**
 * Minimal standalone VR entry module (no dependency on deprecated legacy folders).
 * Exports initVRApp(scene?) which sets up a lightweight XR experience,
 * energy placeholder, and exposes a few globals expected by existing UI code.
 * This is a trimmed version – advanced HUD, selection, physics loops, and
 * bond rotation logic from the legacy implementation can be incrementally
 * reintroduced here (or split into smaller modules under public/vr/).
 */

// NOTE: We purposely avoid importing anything from deleted legacy folder.
// We rely only on Babylon (already loaded globally on the page) and the
// existing non‑VR app setup utilities if present.

// Ensure we reuse any existing viewer engine/scene; never import deprecated setup code.
async function ensureBaseScene() {
	// Prefer an existing global viewer API (desktop mode already running)
	try {
		if (window._viewer && window._viewer.engine && window._viewer.scene) {
			return { engine: window._viewer.engine, scene: window._viewer.scene, reuse: true };
		}
	} catch {}
	// Secondary: previously stored vrEngine/vrScene
	if (window && window.vrEngine && window.vrScene) {
		return { engine: window.vrEngine, scene: window.vrScene, reuse: true };
	}
	// Fallback: create minimal fresh engine/scene
	const canvas = document.getElementById('viewer') || document.querySelector('canvas') || (function(){
		const c = document.createElement('canvas');
		c.id = 'vrCanvasFallback';
		document.body.appendChild(c);
		return c;
	})();
	const engine = new BABYLON.Engine(canvas, false, { antialias: false });
	const scene = new BABYLON.Scene(engine);
	scene.clearColor = new BABYLON.Color4(0.043,0.059,0.078,1.0); // match desktop dark background (#0b0f14)
	new BABYLON.HemisphericLight('vrHemi', new BABYLON.Vector3(0,1,0), scene);
	const cam = new BABYLON.ArcRotateCamera('vrCam', Math.PI/4, Math.PI/3, 8, BABYLON.Vector3.Zero(), scene);
	cam.attachControl(canvas, true);
	return { engine, scene, canvas, reuse: false };
}

async function setupXR(engine, scene, { sessionMode='immersive-vr', referenceSpaceType='local-floor', optionalFeatures=[] } = {}) {
	// Reuse helper only if same sessionMode already active
	if (window.vrHelper && window.vrHelper.baseExperience && window.vrHelper._sessionMode === sessionMode) {
		console.log('[XR] Reusing existing XR helper ('+sessionMode+')');
		return window.vrHelper;
	}
	try {
		const helper = await scene.createDefaultXRExperienceAsync({
			uiOptions: { sessionMode, referenceSpaceType },
			disableTeleportation: true,
			optionalFeatures
		});
		helper._sessionMode = sessionMode; // tag for reuse discrimination
		return helper;
	} catch (e) {
		console.warn('[XR] WebXR initialization failed ('+sessionMode+'):', e.message || e);
		return null;
	}
}

function createOverlay(scene) {
	if (!BABYLON.GUI || !BABYLON.GUI.AdvancedDynamicTexture) return null;
	const adt = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI('vrLiteHUD', false, scene);
	try { adt.useInvalidateRectOptimization = false; } catch {}
	const panel = new BABYLON.GUI.Rectangle('vrLiteEnergyPanel');
	panel.width = '260px';
	panel.height = '100px';
	panel.cornerRadius = 12;
	panel.thickness = 0;
	panel.background = 'rgba(15,18,24,0.85)';
	panel.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_RIGHT;
	panel.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_TOP;
	panel.top = '20px';
	panel.left = '-20px';
	adt.addControl(panel);
	const title = new BABYLON.GUI.TextBlock('vrLiteEnergyTitle', 'Energy');
	title.color = '#a4b0c0';
	title.fontSize = 18;
	title.fontWeight = 'bold';
	title.top = '-28px';
	panel.addControl(title);
	const value = new BABYLON.GUI.TextBlock('vrLiteEnergyValue', '—');
	value.color = '#f5ffd1';
	value.fontSize = 32;
	panel.addControl(value);
	return { adt, value, panel };
}

async function initXRLike({ sessionMode='immersive-vr', transparent=false } = {}) {
	console.log('[XR] initXRLike', { sessionMode, transparent });
	const { engine, scene, reuse } = await ensureBaseScene();
	if (transparent) {
		try {
			// Make scene / canvas transparent (AR passthrough or pseudo-AR fallback)
			scene.autoClear = false;
			scene.clearColor = new BABYLON.Color4(0,0,0,0);
			const canvas = engine.getRenderingCanvas();
			if (canvas) canvas.style.background = 'transparent';
			if (document && document.body) document.body.style.background = 'transparent';
		} catch {}
	}
	const xr = await setupXR(engine, scene, { sessionMode });
	if (!xr) {
		console.warn('[XR] '+sessionMode+' not available');
	}
	const hud = createOverlay(scene);

	// Try to compute one energy value if mlip is available
	try {
		const mlip = window.vrMlip || window.mlip;
		if (mlip && typeof mlip.compute === 'function') {
			const r = mlip.compute();
			if (r && typeof r.then === 'function') {
				r.then(val => {
					if (val && Number.isFinite(val.energy) && hud?.value) hud.value.text = val.energy.toFixed(3);
				}).catch(()=>{});
			} else if (r && Number.isFinite(r.energy) && hud?.value) {
				hud.value.text = r.energy.toFixed(3);
			}
		}
	} catch {}

	window.vrEngine = engine;
	window.vrScene = scene;
	window.vrHelper = xr;
	window.VR_DEBUG = true; // enable verbose debug globally by default

	// Basic render loop (no physics / advanced controls yet)
	if (!reuse) {
		engine.runRenderLoop(() => {
			try { scene.render(); } catch (e) { /* keep loop alive */ }
		});
		window.addEventListener('resize', () => engine.resize());
	}

 	return { engine, scene, xr, hud, sessionMode };
}

// Backwards compatible VR-specific init
export async function initVRApp() { return initXRLike({ sessionMode:'immersive-vr' }); }
// New AR (transparent) init; will attempt immersive-ar, fallback to transparent pseudo-AR if unsupported.
export async function initARApp() {
	if (navigator?.xr && await navigator.xr.isSessionSupported?.('immersive-ar')) {
		return initXRLike({ sessionMode:'immersive-ar', transparent:true });
	} else {
		console.warn('[AR] immersive-ar not supported; enabling transparent background fallback');
		return initXRLike({ sessionMode:'immersive-vr', transparent:true });
	}
}

export async function initXRApp(opts){ return initXRLike(opts||{}); }

export default { initVRApp, initARApp, initXRApp };
