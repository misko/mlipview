// New VR implementation replacing legacy logic. Provides setupVR, setupVRFeatures, and backward-compatible createVRSupport.
import {
  getMoleculeMasters,
  getMoleculeDiag,
  transformLocalToWorld,
  getAnyMaster,
  refreshMoleculeMasters,
} from './vr-utils.js';
import { computeRadialDelta } from './spherical-drag-core.js';
import {
  applyIncrement as vrApplyRotIncrement,
  identityQ as vrIdentityQ,
} from './rotation-core.js';
import { createVRPicker } from './vr-picker.js';
import {
  createEmptySelection,
  applyBondClick,
  bondOrientationColor,
  selectAtom as modelSelectAtom,
} from '../selection-model.js';
import { createVRLaserManager } from './vr-laser.js';
import {
  ensureWorldHUD as ensureXRHudBars,
  forceWorldHUD as forceXRHudBars,
} from './xr-hud-bars.js';

// Minimal guarded logger (kept lightweight for potential future enablement via window.vrLog)
function vrLog() {
  try {
    if (typeof window !== 'undefined' && window.vrLog) console.log.apply(console, arguments);
  } catch (_) {}
}
function vrDebug() {
  try {
    if (typeof window !== 'undefined' && window.vrLog)
      (console.debug || console.log).apply(console, arguments);
  } catch (_) {}
}

// Legacy accumulators replaced by quaternion-based approach for stable orientation.
let _accRotQ = null; // {x,y,z,w}
let _rotMigratedToRoot = false; // one-time hoist of any pre-root master rotation
// Regression guard: rotations must apply only to 'molecule_root' so late-created masters inherit orientation.
// Spherical drag defaults (override via window.*):
//  vrSphericalDrag=true
//  vrSphericalRadialMode='adaptive'
//  vrSphericalWorkingRadiusCap=20
//  vrSphericalAdaptiveGain=7.0
//  vrSphericalAdaptiveMaxDeltaFraction=1.2
//  vrSphericalMinRadius=0.05
//  vrSphericalMaxScaleFactor=4
//  vrSphericalRadialSmoothing=0.25
//  vrSphericalForwardGain=6.0 (forward mode)
//  vrSphericalExpK=1.4 / vrSphericalExpGain=0.9 (exp mode)
//  vrSphericalDebugFrames=12

export async function setupVR(engine, scene, picking) {
  vrLog('[VR Setup] Initializing...');
  if (!scene?.createDefaultXRExperienceAsync) {
    console.warn('[VR Setup] Missing createDefaultXRExperienceAsync');
    return null;
  }
  try {
    engine?.setHardwareScalingLevel?.(1.0);
    scene.skipPointerMovePicking = true;
    scene.autoClear = true;
    scene.autoClearDepthAndStencil = true;
    const xrHelper = await scene.createDefaultXRExperienceAsync({
      floorMeshes: [],
      optionalFeatures: ['hand-tracking', 'local-floor', 'hit-test', 'dom-overlay', 'unbounded'],
      uiOptions: { sessionMode: 'immersive-vr', referenceSpaceType: 'local-floor' },
      // Explicitly enable pointer selection (Babylon normally adds it, but we state intent for clarity)
      inputOptions: {
        doNotLoadControllerMeshes: false,
        disableControllerAnimation: false,
        forceInputProfile: true,
      },
      disableTeleportation: true,
      useStablePlugins: false,
    });
    // Remove Babylon's default floating XR button; we expose XR entry via in-app controls instead.
    try {
      const suppressDefaultButton = () => {
        try {
          const overlay = document.querySelector('.xr-button-overlay');
          if (overlay?.parentElement) overlay.parentElement.removeChild(overlay);
        } catch {}
        try {
          const icon = document.querySelector('.babylonVRicon');
          if (icon?.parentElement) icon.parentElement.removeChild(icon);
        } catch {}
      };
      suppressDefaultButton();
      setTimeout(suppressDefaultButton, 0);
      if (typeof requestAnimationFrame === 'function') {
        requestAnimationFrame(() => suppressDefaultButton());
      }
      if (xrHelper?.baseExperience?.enterExitUI) {
        try {
          xrHelper.baseExperience.enterExitUI.setEnabled?.(false);
          xrHelper.baseExperience.enterExitUI.useCustomButtons = true;
        } catch {}
      }
      if (typeof MutationObserver !== 'undefined' && typeof document !== 'undefined') {
        const observer = new MutationObserver(() => suppressDefaultButton());
        try {
          observer.observe(document.body, { childList: true, subtree: true });
        } catch {}
        if (typeof window !== 'undefined' && Array.isArray(window.__MLIPVIEW_CLEANUP)) {
          window.__MLIPVIEW_CLEANUP.push(() => observer.disconnect());
        }
      }
    } catch {}
    try {
      if (xrHelper?.pointerSelection) {
        console.log('[XR][Input] pointerSelection already present');
      } else if (xrHelper?.featuresManager) {
        const ps = xrHelper.featuresManager.enableFeature(
          BABYLON.WebXRFeatureName.POINTER_SELECTION,
          'latest',
          { xrInput: xrHelper.input, enablePointerSelectionOnAllControllers: true }
        );
        console.log('[XR][Input] pointerSelection enabled via featuresManager', !!ps);
      }
    } catch (e) {
      console.warn('[XR][Input] pointerSelection enable failed', e);
    }
    return setupVRFeatures(xrHelper, scene, picking);
  } catch (e) {
    console.error('[VR Setup] failed', e?.message || e);
    return null;
  }
}

export function setupVRFeatures(xrHelper, scene, picking) {
  const selection = createEmptySelection();
  if (typeof window !== 'undefined') window.vrSelection = selection;
  const molRef = () => {
    // Prefer provided picking.molState (authoritative simulation state) over global fallbacks
    if (picking?.molState) return picking.molState;
    if (typeof window !== 'undefined') return window.vrMol || window.mol || null;
    return null;
  };
  let picker = null;
  try {
    if (window?._viewer?.view) picker = createVRPicker({ scene, view: window._viewer.view });
  } catch {}
  // VR no longer manages its own highlight meshes; relies on canonical highlight meshes
  // created and managed by moleculeView.js (highlight_atom_mesh, highlight_bond_mesh).
  // These are parented to molecule masters so they follow VR rotations/scales.
  function clearSelection() {
    selection.kind = null;
    selection.data = null;
  }
  xrHelper?.baseExperience?.sessionManager?.onXRSessionInit.add(() => {
    _accRotQ = vrIdentityQ();
    vrLog('VR session start');
    clearSelection();
  });
  xrHelper?.baseExperience?.sessionManager?.onXRSessionEnded.add(() => {
    vrLog('VR session end');
    clearSelection();
  });
  // Always expose a lightweight manual snapshot helper (works even if debug sampling disabled)
  try {
    if (typeof window !== 'undefined' && !window.xrDbgSnap) {
      window.xrDbgSnap = function (tag) {
        try {
          const masters = getMoleculeMasters(scene) || [];
          // Aggregate bounding box
          let bbMin = null,
            bbMax = null;
          const temp = {};
          for (const m of masters) {
            try {
              m.refreshBoundingInfo?.();
              const bi = m.getBoundingInfo();
              if (bi) {
                const mn = bi.boundingBox.minimumWorld,
                  mx = bi.boundingBox.maximumWorld;
                if (mn && mx) {
                  if (!bbMin) {
                    bbMin = { x: mn.x, y: mn.y, z: mn.z };
                    bbMax = { x: mx.x, y: mx.y, z: mx.z };
                  } else {
                    bbMin.x = Math.min(bbMin.x, mn.x);
                    bbMin.y = Math.min(bbMin.y, mn.y);
                    bbMin.z = Math.min(bbMin.z, mn.z);
                    bbMax.x = Math.max(bbMax.x, mx.x);
                    bbMax.y = Math.max(bbMax.y, mx.y);
                    bbMax.z = Math.max(bbMax.z, mx.z);
                  }
                }
              }
            } catch {}
          }
          let center = null,
            diag = null;
          if (bbMin && bbMax) {
            center = {
              x: (bbMin.x + bbMax.x) / 2,
              y: (bbMin.y + bbMax.y) / 2,
              z: (bbMin.z + bbMax.z) / 2,
            };
            const dx = bbMax.x - bbMin.x,
              dy = bbMax.y - bbMin.y,
              dz = bbMax.z - bbMin.z;
            diag = Math.sqrt(dx * dx + dy * dy + dz * dz);
          }
          const snap = masters.map((m) => ({
            name: m.name,
            pos: m.position && { x: m.position.x, y: m.position.y, z: m.position.z },
            scaling: m.scaling && { x: m.scaling.x, y: m.scaling.y, z: m.scaling.z },
            rotQ: m.rotationQuaternion && {
              x: m.rotationQuaternion.x,
              y: m.rotationQuaternion.y,
              z: m.rotationQuaternion.z,
              w: m.rotationQuaternion.w,
            },
          }));
          const cam = scene.activeCamera;
          const camPos = cam?.position && {
            x: cam.position.x,
            y: cam.position.y,
            z: cam.position.z,
          };
          const camInfo = cam && {
            name: cam.name,
            type: cam.getClassName && cam.getClassName(),
            pos: camPos,
            target: cam?.target && { x: cam.target.x, y: cam.target.y, z: cam.target.z },
            fov: cam.fov,
          };
          const xrCam = xrHelper?.baseExperience?.camera;
          const xrPos = xrCam && { x: xrCam.position.x, y: xrCam.position.y, z: xrCam.position.z };
          const xrPose = xrCam && { position: xrPos, realWorldHeight: xrCam.realWorldHeight };
          let dist = null;
          if (center && xrPos) {
            const dx = center.x - xrPos.x,
              dy = center.y - xrPos.y,
              dz = center.z - xrPos.z;
            dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
          }
          // Derived recommendation: targetDiag of ~0.35m; propose scale factor if diag known
          let recommendedScale = null;
          if (diag && diag > 1e-4) {
            const target = 0.35;
            recommendedScale = target / diag;
          }
          const rec = {
            t: performance.now ? performance.now() : Date.now(),
            tag: tag || 'manual',
            sessionMode: (xrHelper && xrHelper._sessionMode) || 'unknown',
            masters: snap,
            bb: { min: bbMin, max: bbMax, center, diag },
            cam: camInfo,
            xrPose,
            distCameraToCenter: dist,
            recommendedScale,
          };
          if (!window._xrDbgSamples) window._xrDbgSamples = [];
          window._xrDbgSamples.push(rec);
          console.log('[XRDBG][manual]', rec);
          return rec;
        } catch (e) {
          console.warn('[XRDBG] manual snapshot failed', e);
        }
      };
      // Calibration helper: re-center & re-scale masters in front of XR camera
      window.xrCalibrate = function ({ targetDiag = 0.35, forward = 0.8 } = {}) {
        try {
          const xrCam = xrHelper?.baseExperience?.camera || scene.activeCamera;
          if (!xrCam) {
            console.warn('[XRDBG] calibrate: no camera');
            return false;
          }
          const masters = getMoleculeMasters(scene) || [];
          if (!masters.length) {
            console.warn('[XRDBG] calibrate: no masters');
            return false;
          }
          window.xrDbgSnap('preCalibrate');
          // Compute current combined BB
          let bbMin = null,
            bbMax = null;
          for (const m of masters) {
            try {
              m.refreshBoundingInfo?.();
              const bi = m.getBoundingInfo();
              if (!bi) continue;
              const mn = bi.boundingBox.minimumWorld,
                mx = bi.boundingBox.maximumWorld;
              if (!bbMin) {
                bbMin = { x: mn.x, y: mn.y, z: mn.z };
                bbMax = { x: mx.x, y: mx.y, z: mx.z };
              } else {
                bbMin.x = Math.min(bbMin.x, mn.x);
                bbMin.y = Math.min(bbMin.y, mn.y);
                bbMin.z = Math.min(bbMin.z, mn.z);
                bbMax.x = Math.max(bbMax.x, mx.x);
                bbMax.y = Math.max(bbMax.y, mx.y);
                bbMax.z = Math.max(bbMax.z, mx.z);
              }
            } catch {}
          }
          if (!bbMin || !bbMax) {
            console.warn('[XRDBG] calibrate: no BB');
            return false;
          }
          const center = {
            x: (bbMin.x + bbMax.x) / 2,
            y: (bbMin.y + bbMax.y) / 2,
            z: (bbMin.z + bbMax.z) / 2,
          };
          const dx = bbMax.x - bbMin.x,
            dy = bbMax.y - bbMin.y,
            dz = bbMax.z - bbMin.z;
          const diag = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1;
          const scaleFactor = targetDiag / diag;
          // Place new center forward of camera
          const fwd = xrCam.getDirection
            ? xrCam.getDirection(BABYLON.Axis.Z).clone()
            : new BABYLON.Vector3(0, 0, 1);
          const camPos = xrCam.position.clone();
          const newCenter = camPos.add(fwd.scale(forward));
          for (const m of masters) {
            // Scale relative to center: translate to origin -> scale -> translate to new center
            if (m.position) {
              m.position.x = newCenter.x + (m.position.x - center.x) * scaleFactor;
              m.position.y = newCenter.y + (m.position.y - center.y) * scaleFactor;
              m.position.z = newCenter.z + (m.position.z - center.z) * scaleFactor;
            }
            if (m.scaling) {
              m.scaling.x *= scaleFactor;
              m.scaling.y *= scaleFactor;
              m.scaling.z *= scaleFactor;
            }
          }
          window.xrDbgSnap('postCalibrate');
          console.log('[XRDBG] calibrate applied', { scaleFactor, targetDiag, forward });
          return true;
        } catch (e) {
          console.warn('[XRDBG] calibrate failed', e);
          return false;
        }
      };
    }
  } catch {}
  // XR debug instrumentation (enabled when window.XR_DEBUG or url contains xrdebug=1)
  try {
    const xrDbgEnabled =
      typeof window !== 'undefined' &&
      (window.XR_DEBUG || (window.location && /[?&]xrdebug=1/.test(window.location.search)));
    if (xrDbgEnabled) {
      if (typeof window !== 'undefined') window._xrDbgSamples = window._xrDbgSamples || [];
      function dbgMoleculeState(tag) {
        try {
          const masters = getMoleculeMasters(scene) || [];
          const diag = getMoleculeDiag(scene);
          const snap = masters.map((m) => {
            let bb = null;
            try {
              m.refreshBoundingInfo?.();
              const bi = m.getBoundingInfo();
              if (bi) {
                bb = {
                  min: bi.boundingBox.minimumWorld && {
                    x: bi.boundingBox.minimumWorld.x,
                    y: bi.boundingBox.minimumWorld.y,
                    z: bi.boundingBox.minimumWorld.z,
                  },
                  max: bi.boundingBox.maximumWorld && {
                    x: bi.boundingBox.maximumWorld.x,
                    y: bi.boundingBox.maximumWorld.y,
                    z: bi.boundingBox.maximumWorld.z,
                  },
                };
              }
            } catch {}
            return {
              name: m.name,
              pos: m.position && { x: m.position.x, y: m.position.y, z: m.position.z },
              scaling: m.scaling && { x: m.scaling.x, y: m.scaling.y, z: m.scaling.z },
              rotQ: m.rotationQuaternion && {
                x: m.rotationQuaternion.x,
                y: m.rotationQuaternion.y,
                z: m.rotationQuaternion.z,
                w: m.rotationQuaternion.w,
              },
              bb,
            };
          });
          const cam = scene.activeCamera;
          const camInfo = cam && {
            name: cam.name,
            type: cam.getClassName && cam.getClassName(),
            pos: cam.position && { x: cam.position.x, y: cam.position.y, z: cam.position.z },
            fov: cam.fov,
            minZ: cam.minZ,
            maxZ: cam.maxZ,
            target: cam.target && { x: cam.target.x, y: cam.target.y, z: cam.target.z },
          };
          const xrCam = xrHelper?.baseExperience?.camera;
          const xrPose = xrCam && {
            xrCamera: true,
            realWorldHeight: xrCam.realWorldHeight,
            position: { x: xrCam.position.x, y: xrCam.position.y, z: xrCam.position.z },
          };
          const rec = {
            t: performance.now ? performance.now() : Date.now(),
            tag,
            diag,
            masterCount: masters.length,
            masters: snap,
            cam: camInfo,
            xrPose,
            sessionMode: (xrHelper && xrHelper._sessionMode) || 'unknown',
          };
          window._xrDbgSamples.push(rec);
          console.log('[XRDBG][' + tag + ']', rec);
        } catch (e) {
          console.warn('[XRDBG] snapshot failed', e);
        }
      }
      // Expose manual trigger
      if (typeof window !== 'undefined')
        window.xrDbgSnap = (tag) => dbgMoleculeState(tag || 'manual');
      let frameSamplesRemaining = 90; // capture first 90 frames after session start
      xrHelper?.baseExperience?.sessionManager?.onXRSessionInit.add(() => {
        dbgMoleculeState('sessionInit');
        frameSamplesRemaining = 90;
      });
      scene?.onBeforeRenderObservable?.add(() => {
        if (frameSamplesRemaining > 0) {
          const tag = 'frame' + (90 - frameSamplesRemaining);
          // Only sample while an XR session is active; avoids spamming desktop frames.
          try {
            const active = xrHelper?.baseExperience?.sessionManager?._xrSession;
            if (active) {
              dbgMoleculeState(tag);
              frameSamplesRemaining--;
            } else frameSamplesRemaining = 0;
          } catch {
            frameSamplesRemaining = 0;
          }
        }
      });
      // Capture after XR session end to compare any transform resets
      xrHelper?.baseExperience?.sessionManager?.onXRSessionEnded.add(() => {
        dbgMoleculeState('sessionEnded');
      });
    }
  } catch {}
  // Expose inspection helpers when debug flags enabled
  try {
    if (typeof window !== 'undefined') {
      window.vrInspectRotation = () => {
        const ms = getMoleculeMasters(scene);
        return ms.map((m) => ({
          name: m.name,
          q: m.rotationQuaternion && {
            x: m.rotationQuaternion.x,
            y: m.rotationQuaternion.y,
            z: m.rotationQuaternion.z,
            w: m.rotationQuaternion.w,
          },
          s: m.scaling && { x: m.scaling.x, y: m.scaling.y, z: m.scaling.z },
          pos: m.position && { x: m.position.x, y: m.position.y, z: m.position.z },
        }));
      };
      window.vrDumpMasters = () => {
        console.log('[VR][masters]', window.vrInspectRotation());
      };
      window.vrRefreshMasters = () => {
        const ms = refreshMoleculeMasters(scene, { force: true });
        console.log(
          '[VR][masters][refresh]',
          ms.map((m) => m.name)
        );
        return ms;
      };
    }
  } catch {}
  const controllerState = new WeakMap();
  const HOLD_DELAY = typeof window !== 'undefined' && window.vrHoldDelay ? window.vrHoldDelay : 220;
  const getControllerNode = (c) => c.pointer || c.grip || c.pointer?.parent || null;
  const getYawPitch = (c) => {
    if (typeof BABYLON === 'undefined') return { yaw: 0, pitch: 0 };
    const n = getControllerNode(c);
    if (n?.getDirection) {
      const f = n.getDirection(BABYLON.Vector3.Forward()).normalize();
      return { yaw: Math.atan2(f.x, f.z), pitch: Math.asin(-Math.max(-1, Math.min(1, f.y))) };
    }
    return { yaw: 0, pitch: 0 };
  };
  const isTrigger = (c) => {
    try {
      const mc = c.motionController;
      if (mc?.getComponent) {
        const t = mc.getComponent('xr-standard-trigger');
        if (typeof t?.pressed === 'boolean') return t.pressed;
      }
      const gp = c.inputSource?.gamepad;
      if (gp?.buttons?.length) return !!gp.buttons[0]?.pressed;
    } catch {}
    return false;
  };
  xrHelper?.input?.onControllerAddedObservable.add((c) => {
    controllerState.set(c, {
      pressed: false,
      lastYaw: 0,
      lastPitch: 0,
      downTime: 0,
      tapHandled: false,
    });
    laserMgr?.ensureForController(c);
  });
  try {
    xrHelper?.input?.onControllerRemovedObservable?.add((c) => {
      const st = controllerState.get(c);
      if (st?.laserMesh?.dispose)
        try {
          st.laserMesh.dispose();
        } catch {}
      controllerState.delete(c);
    });
  } catch {}
  function pickAtomOrBondWithRay(ray) {
    if (!ray) return null;
    // Try existing picker first
    if (picker) {
      try {
        const bond = picker.pickBondWithRay(ray);
        if (bond) {
          // Attempt to get distance via a raw scene pick (best effort)
          let dist = null;
          try {
            if (scene.pickWithRay) {
              const pr = scene.pickWithRay(ray);
              if (pr?.hit) dist = pr.distance;
            }
          } catch {}
          return { kind: 'bond', bond, dist };
        }
      } catch {}
      try {
        const atom = picker.pickAtomWithRay(ray);
        if (atom) {
          let dist = null;
          try {
            if (scene.pickWithRay) {
              const pr = scene.pickWithRay(ray);
              if (pr?.hit) dist = pr.distance;
            }
          } catch {}
          return { kind: 'atom', atom, dist };
        }
      } catch {}
    }
    // Fallback: raw scene pick + external resolver
    try {
      if (scene.pickWithRay) {
        const pr = scene.pickWithRay(ray);
        if (pr && pr.hit) {
          const resolver = picking?.view || picking;
          if (resolver) {
            const bond = resolver.resolveBondPick?.(pr);
            if (bond) {
              return { kind: 'bond', bond, dist: pr.distance };
            }
            const atom = resolver.resolveAtomPick?.(pr);
            if (atom) {
              return { kind: 'atom', atom, dist: pr.distance };
            }
          }
        }
      }
    } catch {}
    return null;
  }
  function externalSelectAtom(idx) {
    try {
      picking?.selectionService?.clickAtom?.(idx);
    } catch {}
    try {
      picking?.selection?.clickAtom?.(idx);
    } catch {}
  }
  function externalSelectBond(b) {
    try {
      picking?.selectionService?.clickBond?.(b);
    } catch {}
    try {
      picking?.selection?.clickBond?.(b);
    } catch {}
  }
  function startAtomDrag(st, controller, atomIdx) {
    if (st.dragging) return false;
    const manipulation = picking?.manipulation;
    const mol =
      picking?.molState || (typeof window !== 'undefined' && (window.vrMol || window.mol));
    if (!manipulation || !mol || !mol.positions || !mol.positions[atomIdx]) return false;
    const localStart = { ...mol.positions[atomIdx] };
    const baseMaster = getAnyMaster(scene);
    let scale = { x: 1, y: 1, z: 1 },
      rotQ = null,
      masterPos = { x: 0, y: 0, z: 0 };
    if (baseMaster) {
      if (baseMaster.scaling)
        scale = { x: baseMaster.scaling.x, y: baseMaster.scaling.y, z: baseMaster.scaling.z };
      if (baseMaster.rotationQuaternion) rotQ = baseMaster.rotationQuaternion;
      if (baseMaster.position)
        masterPos = {
          x: baseMaster.position.x,
          y: baseMaster.position.y,
          z: baseMaster.position.z,
        };
    }
    const controllerNode = controller.pointer || controller.grip || controller.pointer?.parent;
    const controllerOriginStart = controllerNode?.getAbsolutePosition
      ? controllerNode.getAbsolutePosition().clone()
      : new BABYLON.Vector3();
    const SPHERICAL_ENABLED =
      typeof window !== 'undefined' ? window.vrSphericalDrag !== false : true; // default on
    const DRAG_GAIN = typeof window !== 'undefined' && window.vrDragGain ? window.vrDragGain : 50.0; // legacy fallback
    function applyInvRotation(vec) {
      if (!rotQ) return;
      const iq = new BABYLON.Quaternion(-rotQ.x, -rotQ.y, -rotQ.z, rotQ.w);
      vec.applyRotationQuaternionInPlace(iq);
    }

    let sphericalData = null;
    if (SPHERICAL_ENABLED && typeof BABYLON !== 'undefined') {
      try {
        const cam = xrHelper?.baseExperience?.camera || scene.activeCamera;
        if (cam?.position) {
          // Compute atom's world position from localStart using master transform snapshot
          let atomWorld = new BABYLON.Vector3(localStart.x, localStart.y, localStart.z);
          // Apply scale then rotation then translation (mirrors transformLocalToWorld)
          atomWorld.x *= scale.x;
          atomWorld.y *= scale.y;
          atomWorld.z *= scale.z;
          if (rotQ) atomWorld.applyRotationQuaternionInPlace(rotQ);
          atomWorld.addInPlace(new BABYLON.Vector3(masterPos.x, masterPos.y, masterPos.z));
          const center = cam.position.clone();
          const atomVec = atomWorld.subtract(center);
          const initialRadius = atomVec.length();
          const ctrlPos = controllerNode?.getAbsolutePosition
            ? controllerNode.getAbsolutePosition()
            : controllerOriginStart;
          const ctrlDir = controllerNode?.getDirection
            ? controllerNode.getDirection(BABYLON.Vector3.Forward()).normalize()
            : atomVec.normalizeToNew();
          const controllerDist0 = ctrlPos.subtract(center).length();
          const invRotQ = rotQ ? new BABYLON.Quaternion(-rotQ.x, -rotQ.y, -rotQ.z, rotQ.w) : null;
          const ctrlVec0 = ctrlPos.subtract(center); // baseline controller vector
          sphericalData = {
            center,
            initialRadius,
            controllerDist0,
            ctrlDir0: ctrlDir.clone(),
            rotQ,
            invRotQ,
            scale: { ...scale },
            masterPos: { ...masterPos },
            framesLogged: 0,
            // dynamic current radius for smoothing
            currentRadius: initialRadius,
            // precompute to speed update
            atomLocalStart: { ...localStart },
            ctrlVec0,
          };
          // (spherical drag init debug suppressed)
        }
      } catch {}
    }

    function sphericalIntersector() {
      if (!sphericalData) {
        return legacyIntersector();
      }
      try {
        const cam = xrHelper?.baseExperience?.camera || scene.activeCamera;
        if (!cam) return legacyIntersector();
        const node = controller.pointer || controller.grip || controller.pointer?.parent;
        const ctrlPos = node?.getAbsolutePosition
          ? node.getAbsolutePosition()
          : controllerOriginStart;
        const forward = node?.getDirection
          ? node.getDirection(BABYLON.Vector3.Forward()).normalize()
          : sphericalData.ctrlDir0.clone();
        // radial adjustment based on controller distance change
        const ctrlVec = ctrlPos.subtract(sphericalData.center);
        const distCtrl = ctrlVec.length();
        const radialGain =
          typeof window !== 'undefined' && window.vrSphericalRadialGain != null
            ? window.vrSphericalRadialGain
            : 1.0;
        // Radial mode selection:
        //   'distance'  : change in controller distance from head.
        //   'projection': change projected onto baseline controller vector.
        //   'hybrid'    : average(distance, projection).
        //   'adaptive'  : ratio-based scaling relative to initial radius.
        //   'forward'   : projection of (controllerPos - center) onto camera forward delta vs baseline (orientation-invariant push/pull).
        //   'exp'       : exponential mapping of controller distance ratio for large radii compression.
        const mode =
          typeof window !== 'undefined' && window.vrSphericalRadialMode
            ? window.vrSphericalRadialMode
            : 'adaptive';
        const core = computeRadialDelta(
          {
            initialRadius: sphericalData.initialRadius,
            controllerDist0: sphericalData.controllerDist0,
            ctrlVec0: [
              sphericalData.ctrlVec0.x,
              sphericalData.ctrlVec0.y,
              sphericalData.ctrlVec0.z,
            ],
            initialCtrlForward: [
              sphericalData.ctrlDir0.x,
              sphericalData.ctrlDir0.y,
              sphericalData.ctrlDir0.z,
            ],
          },
          {
            controllerPos: [ctrlPos.x, ctrlPos.y, ctrlPos.z],
            controllerForward: [forward.x, forward.y, forward.z],
            cameraPos: [sphericalData.center.x, sphericalData.center.y, sphericalData.center.z],
          },
          {
            mode,
            forwardGain:
              typeof window !== 'undefined' && window.vrSphericalForwardGain != null
                ? window.vrSphericalForwardGain
                : 6.0,
            adaptiveGain:
              typeof window !== 'undefined' && window.vrSphericalAdaptiveGain != null
                ? window.vrSphericalAdaptiveGain
                : 7.0,
            adaptiveMaxFrac:
              typeof window !== 'undefined' && window.vrSphericalAdaptiveMaxDeltaFraction != null
                ? window.vrSphericalAdaptiveMaxDeltaFraction
                : 1.2,
            expK:
              typeof window !== 'undefined' && window.vrSphericalExpK != null
                ? window.vrSphericalExpK
                : 1.4,
            expGain:
              typeof window !== 'undefined' && window.vrSphericalExpGain != null
                ? window.vrSphericalExpGain
                : 0.9,
          }
        );
        let deltaR = core.deltaR;
        // Working radius cap: treat enormous initial radii as a smaller effective workspace for responsiveness.
        const workingCap =
          typeof window !== 'undefined' && window.vrSphericalWorkingRadiusCap != null
            ? window.vrSphericalWorkingRadiusCap
            : 20.0; // meters (new default)
        const effectiveInitialRadius = Math.min(sphericalData.initialRadius, workingCap);
        // Amplification (disabled for adaptive/forward/exp which already reshape delta)
        const amplifyCfg =
          typeof window !== 'undefined' && window.vrSphericalRadialAmplify != null
            ? window.vrSphericalRadialAmplify
            : 0.0;
        const amplificationMultiplier =
          mode === 'adaptive' || mode === 'forward' || mode === 'exp'
            ? 1
            : 1 + amplifyCfg * (effectiveInitialRadius / 1.0);
        let newRadiusTarget =
          sphericalData.initialRadius +
          deltaR *
            radialGain *
            amplificationMultiplier *
            (effectiveInitialRadius / sphericalData.initialRadius);
        // Joystick forward/backward: treat as continuous per-frame speed so holding the stick keeps moving
        let __joyActive = false;
        try {
          const gp = controller?.inputSource?.gamepad;
          if (gp && Array.isArray(gp.axes) && gp.axes.length) {
            // Prefer right stick Y if present, else left stick Y
            const y = gp.axes[3] != null ? gp.axes[3] : gp.axes[1] != null ? gp.axes[1] : 0;
            const dead =
              typeof window !== 'undefined' && window.vrJoystickDeadzone != null
                ? window.vrJoystickDeadzone
                : 0.08;
            if (Math.abs(y) > dead) {
              __joyActive = true;
              const speed =
                -y *
                (typeof window !== 'undefined' && window.vrJoystickSpeed != null
                  ? window.vrJoystickSpeed
                  : 0.12); // meters per frame approx (doubled)
              // Base on current radius so the target advances each frame (continuous motion)
              newRadiusTarget = sphericalData.currentRadius + speed;
            }
          }
        } catch {}
        // Smoothing (exponential): currentRadius = lerp(currentRadius, newRadiusTarget, alpha)
        const smoothing =
          typeof window !== 'undefined' && window.vrSphericalRadialSmoothing != null
            ? window.vrSphericalRadialSmoothing
            : 0.25; // 0..1
        let newRadius =
          sphericalData.currentRadius +
          (newRadiusTarget - sphericalData.currentRadius) * Math.max(0, Math.min(1, smoothing));
        const minR =
          typeof window !== 'undefined' && window.vrSphericalMinRadius != null
            ? window.vrSphericalMinRadius
            : 0.05;
        let maxFactor =
          typeof window !== 'undefined' && window.vrSphericalMaxScaleFactor != null
            ? window.vrSphericalMaxScaleFactor
            : 4.0;
        // While joystick is actively pushed, allow much larger travel (or configurable)
        if (__joyActive) {
          const joyMax =
            typeof window !== 'undefined' && window.vrJoystickMaxScaleFactor != null
              ? window.vrJoystickMaxScaleFactor
              : Infinity;
          if (Number.isFinite(joyMax)) maxFactor = joyMax;
          else maxFactor = 1e9;
        }
        const maxR = sphericalData.initialRadius * maxFactor;
        if (newRadius < minR) newRadius = minR;
        else if (newRadius > maxR) newRadius = maxR;
        sphericalData.currentRadius = newRadius; // persist
        // New world position on sphere (or adjusted radius) along controller forward direction.
        const newWorld = sphericalData.center.add(forward.scale(newRadius));
        // Convert world -> local (inverse of master transform at drag start)
        const local = newWorld.clone();
        // subtract master translation
        local.subtractInPlace(
          new BABYLON.Vector3(
            sphericalData.masterPos.x,
            sphericalData.masterPos.y,
            sphericalData.masterPos.z
          )
        );
        // inverse rotation
        if (sphericalData.invRotQ) local.applyRotationQuaternionInPlace(sphericalData.invRotQ);
        // inverse scale
        local.x /= sphericalData.scale.x;
        local.y /= sphericalData.scale.y;
        local.z /= sphericalData.scale.z;
        // (per-frame spherical drag debug suppressed)
        return { x: local.x, y: local.y, z: local.z };
      } catch (e) {
        return legacyIntersector();
      }
    }

    function legacyIntersector() {
      const node = controller.pointer || controller.grip || controller.pointer?.parent;
      const cur = node?.getAbsolutePosition ? node.getAbsolutePosition() : controllerOriginStart;
      const worldDelta = cur.subtract(controllerOriginStart).scale(DRAG_GAIN);
      const d = worldDelta.clone();
      applyInvRotation(d);
      d.x /= scale.x;
      d.y /= scale.y;
      d.z /= scale.z;
      return { x: localStart.x + d.x, y: localStart.y + d.y, z: localStart.z + d.z };
    }

    const intersector = SPHERICAL_ENABLED ? sphericalIntersector : legacyIntersector;
    const started = manipulation.beginDrag(intersector, { source: 'vr' });
    if (started) {
      st.dragging = true;
      st.dragAtomIndex = atomIdx;
      st.dragIntersector = intersector;
      return true;
    }
    return false;
  }
  function updateDrag(st) {
    if (!st.dragging) return;
    try {
      picking?.manipulation?.updateDrag?.(st.dragIntersector);
    } catch {}
  }
  function endDrag(st) {
    if (!st.dragging) return;
    try {
      picking?.manipulation?.endDrag?.();
    } catch {}
    st.dragging = false;
    st.dragAtomIndex = null;
    st.dragIntersector = null;
  }

  // Helper: build a picking ray from a controller node
  function rayFromController(c) {
    try {
      const n = c && (c.pointer || c.grip || c.pointer?.parent);
      if (n?.getDirection && n.getAbsolutePosition) {
        const o = n.getAbsolutePosition();
        const d = n.getDirection(BABYLON.Vector3.Forward()).normalize();
        return new BABYLON.Ray(o, d, 120);
      }
    } catch {}
    return null;
  }

  // Helper: apply selection result (atom or bond) reused for press + release flows
  function applySelectionResult(res, st) {
    if (!res) {
      return;
    }
    if (res.kind === 'atom') {
      const idx = res.atom.idx ?? res.atom.index;
      // Set initialPick BEFORE any optional visualization (removed) so drag hold logic always sees it.
      st.initialPick = st.initialPick || { kind: 'atom', index: idx };
      modelSelectAtom(selection, { idx });
      externalSelectAtom(idx);
      // Highlight visualization now handled entirely by moleculeView.js; removed legacy VR marker call.
    } else if (res.kind === 'bond') {
      const b = res.bond;
      const prev = selection.kind;
      const state = applyBondClick(selection, { i: b.i, j: b.j, key: b.key, index: b.index });
      if (state === 'cleared') {
        clearSelection();
      } else if (selection.kind === 'bond') {
        externalSelectBond(b);
        const col = bondOrientationColor(selection.data.orientation);
        const m = molRef();
        let ai = m?.atoms?.[selection.data.i]?.pos;
        let aj = m?.atoms?.[selection.data.j]?.pos;
        if ((!ai || !aj) && m?.positions) {
          const pi = m.positions[selection.data.i];
          const pj = m.positions[selection.data.j];
          if (pi && pj) {
            ai = pi;
            aj = pj;
          }
        }
        // Bond highlight cylinder handled by moleculeView.js selection system.
        if (prev !== 'bond')
          vrDebug('[VR Select] Bond ' + selection.data.i + '-' + selection.data.j);
        st.initialPick = st.initialPick || { kind: 'bond', ...b };
      }
    }
  }

  const laserMgr = createVRLaserManager({
    scene,
    xrHelper,
    controllerState,
    getControllerNode,
    pickAtomOrBondWithRay,
  });
  scene?.onBeforeRenderObservable?.add(() => {
    // Bond highlight frame updates delegated to moleculeView.js
    // Refresh masters if we only have highlight meshes cached
    let masters = getMoleculeMasters(scene);
    if (masters.length && masters.every((m) => /^highlight_/.test(m.name))) {
      masters = refreshMoleculeMasters(scene, { force: true });
    }
    const noMasters = !masters?.length;
    // Even if no masters (e.g., early init or test harness), still perform laser updates below.
    // Frame counter no longer needed after debug removal
    const inputs = xrHelper?.input?.controllers || [];
    const pressedWithPos = [];
    for (const c of inputs)
      if (isTrigger(c)) {
        const node = getControllerNode(c);
        if (node?.getAbsolutePosition)
          pressedWithPos.push({ c, p: node.getAbsolutePosition().clone() });
      }
    // Update lasers every frame regardless of masters
    try {
      laserMgr.updateFrame(inputs);
    } catch {}
    if (noMasters) {
      // Minimal laser diameter relaxation when nothing loaded
      for (const c of inputs) {
        let st = controllerState.get(c);
        if (!st) {
          st = {
            pressed: false,
            lastYaw: 0,
            lastPitch: 0,
            downTime: 0,
            tapHandled: false,
            dragging: false,
          };
          controllerState.set(c, st);
        }
        if (typeof st.laserDiameterCurrent !== 'number') st.laserDiameterCurrent = LASER_BASE_DIAM;
        st.laserDiameterCurrent =
          st.laserDiameterCurrent + (LASER_BASE_DIAM - st.laserDiameterCurrent) * LASER_LERP;
      }
    }
    // If no masters we still proceed for laser hover & controller state but skip molecule transform logic below via guards.
    const twoHandMode = pressedWithPos.length >= 2;
    // Two-hand stretch (pinch zoom) based on distance change between controllers
    if (twoHandMode) {
      try {
        const p0 = pressedWithPos[0].p,
          p1 = pressedWithPos[1].p;
        const dist = p0.subtract(p1).length();
        if (!scene._vrTwoHandLastDist) {
          scene._vrTwoHandLastDist = dist;
          scene._vrTwoHandBaseScales = masters.map(
            (m) => m.scaling?.clone?.() || new BABYLON.Vector3(1, 1, 1)
          );
          // scale start
        } else {
          const baseDist = scene._vrTwoHandLastDist;
          if (baseDist > 1e-5) {
            const ratio = dist / baseDist; // uniform scaling ratio
            const maxRatio = 5,
              minRatio = 0.05;
            const clamped = Math.max(minRatio, Math.min(maxRatio, ratio));
            let before = null;
            masters.forEach((m, i) => {
              const bs = scene._vrTwoHandBaseScales[i] || m.scaling || BABYLON.Vector3.One();
              if (m.scaling) {
                m.scaling.x = bs.x * clamped;
                m.scaling.y = bs.y * clamped;
                m.scaling.z = bs.z * clamped;
              }
            });
            // scale update
          }
        }
      } catch (e) {
        /* silent */
      }
    } else {
      scene._vrTwoHandLastDist = null;
      scene._vrTwoHandBaseScales = null;
    }
    // Determine authoritative rotation target: prefer single molecule_root if present in masters
    let moleculeRoot = null;
    try {
      moleculeRoot = masters.find((m) => m && m.name === 'molecule_root');
    } catch {}
    // One-time migration: if we have already rotated individual masters before root existed, consolidate
    if (moleculeRoot && !_rotMigratedToRoot) {
      try {
        const donor = masters.find((m) => m.rotationQuaternion);
        if (donor) {
          moleculeRoot.rotationQuaternion = donor.rotationQuaternion.clone();
          _accRotQ = {
            x: moleculeRoot.rotationQuaternion.x,
            y: moleculeRoot.rotationQuaternion.y,
            z: moleculeRoot.rotationQuaternion.z,
            w: moleculeRoot.rotationQuaternion.w,
          };
        } else if (!moleculeRoot.rotationQuaternion) {
          moleculeRoot.rotationQuaternion = BABYLON.Quaternion.Identity();
        }
        for (const m of masters) {
          if (m !== moleculeRoot && m.rotationQuaternion) {
            m.rotationQuaternion = null;
          }
        }
        _rotMigratedToRoot = true;
      } catch {}
    }
    for (const c of inputs) {
      let st = controllerState.get(c);
      if (!st) {
        st = {
          pressed: false,
          lastYaw: 0,
          lastPitch: 0,
          downTime: 0,
          tapHandled: false,
          dragging: false,
        };
        controllerState.set(c, st);
        // If controller added outside of onControllerAddedObservable (test shim), initialize laser state
        if (typeof st.laserDiameterCurrent === 'undefined') {
          st.laserDiameterCurrent = LASER_BASE_DIAM;
          st.hoverKind = null;
        }
      }
      // Lazy-create laser mesh if missing (e.g., controller injected after observable)
      if (typeof BABYLON !== 'undefined' && !st.laserMesh) {
        try {
          const node = getControllerNode(c);
          if (node) {
            const laser = BABYLON.MeshBuilder.CreateCylinder(
              'xrLaser_' + (c.uniqueId || 'u'),
              { height: LASER_LENGTH, diameter: LASER_BASE_DIAM, tessellation: 6, subdivisions: 1 },
              scene
            );
            laser.isPickable = false;
            laser.parent = node;
            laser.position = new BABYLON.Vector3(0, 0, LASER_LENGTH / 2);
            const mat = new BABYLON.StandardMaterial('xrLaserMat_' + (c.uniqueId || 'u'), scene);
            mat.disableLighting = true;
            mat.emissiveColor = new BABYLON.Color3(
              LASER_BASE_COLOR.r,
              LASER_BASE_COLOR.g,
              LASER_BASE_COLOR.b
            );
            mat.alpha = 0.95;
            laser.material = mat;
            st.laserMesh = laser;
            if (typeof st.laserDiameterCurrent !== 'number')
              st.laserDiameterCurrent = LASER_BASE_DIAM;
          }
        } catch {}
      }
      const pressed = isTrigger(c);
      if (noMasters) {
        // Skip selection/drag logic entirely without masters; still perform hover ray & diameter update
      }
      if (pressed && !st.pressed) {
        const { yaw, pitch } = getYawPitch(c);
        st.lastYaw = yaw;
        st.lastPitch = pitch;
        st.pressed = true;
        st.downTime = performance.now ? performance.now() : Date.now();
        st.initialPick = null;
        try {
          const ray = rayFromController(c);
          if (ray) {
            const res = pickAtomOrBondWithRay(ray);
            applySelectionResult(res, st);
          }
        } catch {}
      } else if (!pressed && st.pressed) {
        const upTime = performance.now ? performance.now() : Date.now();
        const duration = upTime - (st.downTime || upTime);
        if (st.dragging) {
          endDrag(st);
        } else if (duration <= 300) {
          // tap release selects (already selected on press for atoms; ensure bond tap works if not applied)
          try {
            if (!st.initialPick) {
              const ray = rayFromController(c);
              if (ray) {
                const res = pickAtomOrBondWithRay(ray);
                applySelectionResult(res, st);
              }
            }
          } catch {}
        }
        st.pressed = false;
        if (st.dragging) endDrag(st);
      }
      // Hold-to-drag detection
      if (st.pressed && !st.dragging && st.initialPick?.kind === 'atom') {
        const now = performance.now ? performance.now() : Date.now();
        if (!twoHandMode && now - st.downTime >= HOLD_DELAY) {
          startAtomDrag(st, c, st.initialPick.index);
        }
      }
      if (st.pressed && !twoHandMode && typeof BABYLON !== 'undefined') {
        // If dragging, skip only the accumulation/application, but still keep yaw/pitch baseline current so rotation won't jump afterwards.
        if (st.dragging) {
          const { yaw, pitch } = getYawPitch(c);
          st.lastYaw = yaw;
          st.lastPitch = pitch;
          // rotation suppressed during drag
        } else {
          const sens = 0.9;
          const { yaw, pitch } = getYawPitch(c);
          let dYaw = yaw - st.lastYaw;
          let dPitch = pitch - st.lastPitch;
          if (dYaw > Math.PI) dYaw -= 2 * Math.PI;
          else if (dYaw < -Math.PI) dYaw += 2 * Math.PI;
          if (dPitch > Math.PI) dPitch -= 2 * Math.PI;
          else if (dPitch < -Math.PI) dPitch += 2 * Math.PI;
          const max = 0.05;
          dYaw = Math.max(-max, Math.min(max, dYaw));
          dPitch = Math.max(-max, Math.min(max, dPitch));
          // Initialize accumulator if missing
          if (!_accRotQ) _accRotQ = vrIdentityQ();
          // Camera-based right axis for pitch so vertical drag feels consistent after yaw rotations
          let cam = scene.activeCamera;
          let rightAxis = [1, 0, 0];
          if (cam?.getDirection) {
            try {
              const r = cam.getDirection(BABYLON.Axis.X);
              rightAxis = [r.x, r.y, r.z];
            } catch {}
          }
          // Apply incremental update in a quaternion accumulator (math in rotation-core.js)
          _accRotQ = vrApplyRotIncrement(_accRotQ, -dYaw, -dPitch, rightAxis, {
            sens,
            maxStep: max,
          });
          const q = new BABYLON.Quaternion(_accRotQ.x, _accRotQ.y, _accRotQ.z, _accRotQ.w);
          if (moleculeRoot) {
            if (!moleculeRoot.rotationQuaternion)
              moleculeRoot.rotationQuaternion = BABYLON.Quaternion.Identity();
            moleculeRoot.rotationQuaternion.copyFrom(q);
          } else {
            for (const m of masters) m.rotationQuaternion = q.clone();
          }
          st.lastYaw = yaw;
          st.lastPitch = pitch;
        }
      }
      // Joystick bond rotation (continuous): if a bond is selected, map stick X to delta angle per frame.
      // Uses the same rotation logic as desktop by delegating to manipulation.rotateBond(deltaAngle).
      // Active whenever a controller is present (no trigger hold required), skipped during two-hand scaling or atom drag.
      try {
        const sel = picking?.molState?.selection;
        if (sel && sel.kind === 'bond' && !twoHandMode && !st.dragging) {
          const gp = c?.inputSource?.gamepad;
          if (gp && Array.isArray(gp.axes)) {
            // Prefer right stick X if present, otherwise left stick X
            const x = gp.axes[2] != null ? gp.axes[2] : gp.axes[0] != null ? gp.axes[0] : 0;
            const dead =
              typeof window !== 'undefined' && window.vrBondRotateDeadzone != null
                ? window.vrBondRotateDeadzone
                : 0.08;
            if (Math.abs(x) > dead) {
              const gain =
                typeof window !== 'undefined' && window.vrBondRotateGain != null
                  ? window.vrBondRotateGain
                  : 0.08; // rad/frame @ full deflection
              const delta = x * gain;
              try {
                picking?.manipulation?.rotateBond?.(delta);
              } catch {}
            }
          }
        }
      } catch {}
      // Suppress drag updates entirely in two-hand mode (scaling takes precedence)
      if (!noMasters && st.dragging && !twoHandMode) updateDrag(st);
    }
    // If we entered two-hand mode this frame and any controller was dragging, end drags
    if (twoHandMode) {
      for (const c of inputs) {
        const s = controllerState.get(c);
        if (s?.dragging) endDrag(s);
      }
    }
    if (twoHandMode && typeof BABYLON !== 'undefined') {
      const center = pressedWithPos[0].p.add(pressedWithPos[1].p).scale(0.5);
      if (!scene._vrLastCenter) scene._vrLastCenter = center.clone();
      else {
        let delta = center.subtract(scene._vrLastCenter);
        const cam = scene.activeCamera;
        if (cam?.getDirection) {
          const right = cam.getDirection(BABYLON.Axis.X);
          const up = cam.getDirection(BABYLON.Axis.Y);
          const fwd = cam.getDirection(BABYLON.Axis.Z);
          const dx = BABYLON.Vector3.Dot(delta, right);
          const dy = BABYLON.Vector3.Dot(delta, up);
          const dz = BABYLON.Vector3.Dot(delta, fwd);
          delta = right
            .scale(dx * 16)
            .addInPlace(up.scale(dy * 16))
            .addInPlace(fwd.scale(dz * 60));
        } else delta.scaleInPlace(5);
        const maxStep = getMoleculeDiag(scene) * 2.0;
        if (delta.length() > maxStep) delta.normalize().scaleInPlace(maxStep);
        for (const m of masters) m.position.addInPlace(delta);
        scene._vrLastCenter.copyFrom(center);
      }
    } else scene._vrLastCenter = null;
    try {
      const mol = molRef();
      // Frame-level updates for highlights are handled centrally by moleculeView.js now.
      // (VR still responsible for rotation/scale/translation of molecule masters.)
      if (mol && selection.kind === 'bond') {
        /* debug hook retained */ if (bondDbg)
          console.log('[VR][sentinel] frame bond branch (delegated)');
      }
    } catch {}
  });
  return xrHelper;
}

function enableStandardXRInteractions(scene, xrHelper) {
  if (typeof BABYLON === 'undefined') return false;
  try {
    const fm = xrHelper?.baseExperience?.featuresManager;
    if (!fm) return false;
    try {
      fm.enableFeature(BABYLON.WebXRFeatureName.POINTER_SELECTION, 'latest', {
        xrInput: xrHelper.input,
        forceGazeMode: false,
        disableScenePointerVectorUpdate: false,
      });
    } catch {}
    try {
      fm.enableFeature(BABYLON.WebXRFeatureName.NEAR_INTERACTION, 'latest');
    } catch {}
    const masters = getMoleculeMasters(scene);
    if (!masters?.length) return false;
    let min = new BABYLON.Vector3(
      Number.POSITIVE_INFINITY,
      Number.POSITIVE_INFINITY,
      Number.POSITIVE_INFINITY
    );
    let max = new BABYLON.Vector3(
      Number.NEGATIVE_INFINITY,
      Number.NEGATIVE_INFINITY,
      Number.NEGATIVE_INFINITY
    );
    for (const m of masters) {
      try {
        m.refreshBoundingInfo?.();
        const bi = m.getBoundingInfo();
        const bmin = bi.boundingBox.minimumWorld;
        const bmax = bi.boundingBox.maximumWorld;
        min = BABYLON.Vector3.Minimize(min, bmin);
        max = BABYLON.Vector3.Maximize(max, bmax);
      } catch {}
    }
    const center = min.add(max).scale(0.5);
    const size = max.subtract(min);
    const wrapper = BABYLON.MeshBuilder.CreateBox(
      'molWrapper',
      { width: (size.x || 0.1) * 1.1, height: (size.y || 0.1) * 1.1, depth: (size.z || 0.1) * 1.1 },
      scene
    );
    wrapper.position.copyFrom(center);
    wrapper.rotationQuaternion = BABYLON.Quaternion.Identity();
    wrapper.isPickable = true;
    wrapper.visibility = 0.02;
    wrapper.alwaysSelectAsActiveMesh = true;
    const drag = new BABYLON.SixDofDragBehavior();
    wrapper.addBehavior(drag);
    const scaler = new BABYLON.MultiPointerScaleBehavior();
    wrapper.addBehavior(scaler);
    return true;
  } catch (e) {
    console.warn('[VR Behaviors] failed', e?.message || e);
    return false;
  }
}

// Backwards compatible API expected by existing code/tests
export function createVRSupport(scene, { picking } = {}) {
  let xrHelper = null,
    supported = false,
    engineRef = null;
  const controllerStates = new Map(); // compatibility for legacy tests
  let _shimControllers = [];
  function isSupported() {
    return supported;
  }
  function controllers() {
    const real = xrHelper?.input?.controllers;
    if (real && real.length) return real;
    if (real && !real.length && _shimControllers.length) return _shimControllers;
    return real || _shimControllers;
  }
  async function supportsSession(mode) {
    try {
      if (!navigator?.xr?.isSessionSupported) return false;
      return await navigator.xr.isSessionSupported(mode);
    } catch {
      return false;
    }
  }
  function _rawSession() {
    try {
      return (
        xrHelper?.baseExperience?.sessionManager?.session ||
        xrHelper?.baseExperience?.sessionManager?._session ||
        xrHelper?.baseExperience?.sessionManager?._xrSession ||
        null
      );
    } catch {
      return null;
    }
  }
  function _sessionModeInternal() {
    try {
      if (xrHelper?._sessionMode) return xrHelper._sessionMode;
      const rs = _rawSession();
      if (rs && rs.environmentBlendMode) return 'immersive';
    } catch {}
    return 'unknown';
  }
  function setTransparent(on) {
    try {
      if (!scene) return;
      const ccBefore = scene.clearColor
        ? {
            r: scene.clearColor.r,
            g: scene.clearColor.g,
            b: scene.clearColor.b,
            a: scene.clearColor.a,
          }
        : null;
      console.log('[VR][AR toggle] applying', {
        on,
        autoClearBefore: scene.autoClear,
        clearColorBefore: ccBefore,
      });
      if (on) {
        scene.autoClear = false;
        if (scene.clearColor && typeof BABYLON !== 'undefined' && BABYLON.Color4)
          scene.clearColor = new BABYLON.Color4(0, 0, 0, 0);
        const canvas =
          scene.getEngine?.()?.getRenderingCanvas?.() || document.getElementById('viewer');
        if (canvas) canvas.style.background = 'transparent';
        if (typeof document !== 'undefined' && document.body)
          document.body.style.background = 'transparent';
      } else {
        scene.autoClear = true;
        if (scene.clearColor && typeof BABYLON !== 'undefined' && BABYLON.Color4)
          scene.clearColor = new BABYLON.Color4(0.043, 0.059, 0.078, 1.0);
        const canvas =
          scene.getEngine?.()?.getRenderingCanvas?.() || document.getElementById('viewer');
        if (canvas) canvas.style.background = '#0b0f14';
        if (typeof document !== 'undefined' && document.body)
          document.body.style.background = '#0b0f14';
      }
      const ccAfter = scene.clearColor
        ? {
            r: scene.clearColor.r,
            g: scene.clearColor.g,
            b: scene.clearColor.b,
            a: scene.clearColor.a,
          }
        : null;
      console.log('[VR][AR toggle] applied', {
        on,
        autoClearAfter: scene.autoClear,
        clearColorAfter: ccAfter,
      });
    } catch (e) {
      /* silent */
    }
  }
  function activeSession() {
    try {
      return _rawSession();
    } catch {
      return null;
    }
  }
  function sessionMode() {
    try {
      const s = activeSession();
      if (!s) return _sessionModeInternal();
      if (xrHelper?._sessionMode) return xrHelper._sessionMode;
      return _sessionModeInternal();
    } catch {
      return 'unknown';
    }
  }
  function environmentBlend() {
    try {
      return activeSession()?.environmentBlendMode || 'unknown';
    } catch {
      return 'unknown';
    }
  }
  // --- BEGIN cleanup logging helpers (non-breaking) ---
  const __vrLogVerbose = () => {
    try {
      return (
        typeof window !== 'undefined' &&
        (window.VR_VERBOSE || /[?&]vrverbose=1/.test(window.location?.search || ''))
      );
    } catch {
      return false;
    }
  };
  function vrLogInfo(msg, meta) {
    try {
      if (__vrLogVerbose())
        meta !== undefined ? console.log('[VR]', msg, meta) : console.log('[VR]', msg);
    } catch {}
  }
  function vrLogWarn(msg, meta) {
    try {
      meta !== undefined ? console.warn('[VR]', msg, meta) : console.warn('[VR]', msg);
    } catch {}
  }
  // --- END cleanup logging helpers ---
  async function enterAR() {
    console.log('[VR][AR] enterAR requested');
    const ok = await supportsSession('immersive-ar');
    console.log('[VR][AR] immersive-ar support', ok);
    if (!ok) {
      console.warn('[VR][AR] immersive-ar not supported; using transparent fallback only');
      return false;
    }
    // If we already have an immersive session (vr), end it first to avoid InvalidStateError
    if (activeSession()) {
      console.log('[VR][AR] existing session detected; exiting before AR');
      try {
        await exitXR();
      } catch {}
    }
    try {
      await xrHelper?.baseExperience?.enterXRAsync?.('immersive-ar', 'local-floor');
      if (xrHelper) xrHelper._sessionMode = 'immersive-ar';
      console.log('[VR][AR] entered immersive-ar session', { blend: environmentBlend() });
      // Auto-normalization for AR: if molecule too close / large, re-scale & offset forward.
      try {
        const { computeMoleculeWorldBounds, getMoleculeMasters } = await import('./vr-utils.js');
        const AUTO_TAG = '[XRAC]';
        const targetDiag = 0.35; // meters desired diag
        const idealForward = 0.8; // meters in front of camera
        const minForward = 0.55; // if closer than this, push out
        const maxForward = 1.25; // clamp forward placement
        let attempts = 0,
          applied = false;
        const applyIfNeeded = () => {
          attempts++;
          const masters = getMoleculeMasters(scene) || [];
          if (!masters.length) {
            if (attempts < 120) return;
            else {
              console.log(AUTO_TAG, 'no masters after retries');
              return;
            }
          }
          // Compute combined bb (robust to root-only)
          let bounds = computeMoleculeWorldBounds(scene, masters);
          if (!bounds) {
            if (attempts < 120) return;
            else {
              console.log(AUTO_TAG, 'empty bb after retries');
              return;
            }
          }
          const bbMin = bounds.min,
            bbMax = bounds.max;
          const center = {
            x: (bbMin.x + bbMax.x) / 2,
            y: (bbMin.y + bbMax.y) / 2,
            z: (bbMin.z + bbMax.z) / 2,
          };
          const dx = bbMax.x - bbMin.x,
            dy = bbMax.y - bbMin.y,
            dz = bbMax.z - bbMin.z;
          const diag = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1;
          const xrCam = xrHelper?.baseExperience?.camera || scene.activeCamera;
          if (!xrCam) {
            console.log(AUTO_TAG, 'no camera yet');
            if (attempts < 120) return;
            else return;
          }
          const camPos = xrCam.position.clone();
          const camToCenter = new BABYLON.Vector3(
            center.x - camPos.x,
            center.y - camPos.y,
            center.z - camPos.z
          );
          const dist = camToCenter.length();
          const scaleFactor = targetDiag / diag;
          const needScale = diag > targetDiag * 1.25 || diag < targetDiag * 0.55; // outside tolerance window
          const needPush = dist < minForward;
          if (!(needScale || needPush)) {
            if (attempts < 90) return;
            else {
              console.log(AUTO_TAG, 'within tolerances; no change');
              applied = true;
              return;
            }
          }
          const fwd = xrCam.getDirection
            ? xrCam.getDirection(BABYLON.Axis.Z).clone()
            : new BABYLON.Vector3(0, 0, 1);
          fwd.normalize();
          const targetDist = Math.min(maxForward, Math.max(idealForward, dist, minForward));
          const desiredCenter = camPos.add(fwd.scale(targetDist));
          const appliedScale = needScale ? scaleFactor : 1.0;
          // If root is present, apply to root only; else fallback to masters
          const root =
            masters.length === 1 && masters[0]?.name === 'molecule_root' ? masters[0] : null;
          if (root) {
            if (root.position) {
              root.position.x = desiredCenter.x + (root.position.x - center.x) * appliedScale;
              root.position.y = desiredCenter.y + (root.position.y - center.y) * appliedScale;
              root.position.z = desiredCenter.z + (root.position.z - center.z) * appliedScale;
            }
            if (needScale && root.scaling) {
              root.scaling.x *= appliedScale;
              root.scaling.y *= appliedScale;
              root.scaling.z *= appliedScale;
            }
          } else {
            for (const m of masters) {
              if (m.position) {
                m.position.x = desiredCenter.x + (m.position.x - center.x) * appliedScale;
                m.position.y = desiredCenter.y + (m.position.y - center.y) * appliedScale;
                m.position.z = desiredCenter.z + (m.position.z - center.z) * appliedScale;
              }
              if (needScale && m.scaling) {
                m.scaling.x *= appliedScale;
                m.scaling.y *= appliedScale;
                m.scaling.z *= appliedScale;
              }
            }
          }
          console.log(AUTO_TAG, 'applied', {
            attempts,
            originalDiag: diag,
            scaleApplied: needScale ? appliedScale : 1.0,
            originalDist: dist,
            targetDist,
            desiredCenter,
          });
          try {
            window.__XR_AUTO_NORM = {
              enabled: true,
              scale: needScale ? appliedScale : 1.0,
              desiredCenter,
              bbCenter: center,
            };
          } catch {}
          applied = true;
        };
        // Retry for up to ~2 seconds (120 frames at 60fps) until masters populated
        let frames = 0;
        const obs = scene?.onBeforeRenderObservable?.add(() => {
          if (applied) {
            scene.onBeforeRenderObservable.remove?.(obs);
            return;
          }
          frames++;
          applyIfNeeded();
          if (applied || frames > 130) {
            scene.onBeforeRenderObservable.remove?.(obs);
          }
        });
      } catch (e) {
        console.warn('[XRAC] auto-normalization failed', e);
      }
      return true;
    } catch (e) {
      console.warn('[VR][AR] enter immersive-ar failed', e?.message || e);
      return false;
    }
  }
  async function enterVR() {
    console.log('[VR] enterVR requested');
    try {
      const support = await supportsSession('immersive-vr');
      console.log('[VR][diag] isSessionSupported(immersive-vr)=', support);
      if (!support) {
        console.warn('[VR] immersive-vr not supported');
        return false;
      }
      if (!xrHelper?.baseExperience?.enterXRAsync) {
        console.warn('[VR] baseExperience.enterXRAsync missing');
        return false;
      }
      const pre = { modeBefore: sessionMode(), activeBefore: !!activeSession() };
      console.log('[VR][diag] pre-enter state', pre);
      await xrHelper.baseExperience.enterXRAsync('immersive-vr', 'local-floor');
      if (xrHelper) xrHelper._sessionMode = 'immersive-vr';
      const post = {
        modeAfter: sessionMode(),
        activeAfter: !!activeSession(),
        blend: environmentBlend(),
      };
      console.log('[VR] entered immersive-vr', post);
      setTimeout(() => {
        try {
          console.log('[VR][diag] post-enter delayed check', {
            active: !!activeSession(),
            mode: sessionMode(),
            blend: environmentBlend(),
          });
        } catch {}
      }, 300);
      return true;
    } catch (e) {
      console.warn('[VR] enter immersive-vr failed', e?.message || e);
      return false;
    }
  }
  async function exitXR() {
    try {
      await xrHelper?.baseExperience?.exitXR?.();
      console.log('[VR] XR session exited');
    } catch (e) {
      console.warn('[VR] exit XR failed', e?.message || e);
    }
  }
  async function switchXR(target) {
    // target: 'vr' | 'ar' | 'none'
    console.log('[XR] switch request', target);
    try {
      console.log('[XR][diag] pre-switch', {
        target,
        active: !!activeSession(),
        mode: sessionMode(),
        blend: environmentBlend(),
        supportedFlag: supported,
      });
    } catch {}
    if (target === 'none') {
      await exitXR();
      return true;
    }
    const wantsAR = target === 'ar';
    const active = !!activeSession();
    if (wantsAR) {
      if (active) {
        await exitXR();
      }
      const ok = await enterAR();
      if (ok) setTransparent(true); // auto-enable transparency in AR
      if (ok)
        try {
          console.log('[XR][HUD] attempting create (AR branch)');
          ensureXRHUD();
        } catch (e) {
          console.warn('[XR][HUD] create failed', e);
        }
      if (ok) {
        try {
          setTimeout(() => {
            try {
              forceWorldHUD();
            } catch (e) {
              console.warn('[XR][HUD] auto worldHUD failed', e);
            }
          }, 600);
        } catch {}
      }
      return ok;
    } else {
      if (active) {
        console.log('[XR] existing session before VR; exiting first');
        await exitXR();
      }
      const okVR = await enterVR();
      console.log('[XR][diag] after enterVR', {
        okVR,
        active: !!activeSession(),
        mode: sessionMode(),
        blend: environmentBlend(),
      });
      // Restore opaque background for VR unless user manually toggled otherwise
      setTransparent(false);
      try {
        console.log('[XR][HUD] attempting create (VR branch)');
        ensureXRHUD();
      } catch (e) {
        console.warn('[XR][HUD] create failed', e);
      }
      try {
        setTimeout(() => {
          try {
            forceWorldHUD();
          } catch (e) {
            console.warn('[XR][HUD] auto worldHUD failed', e);
          }
        }, 600);
      } catch {}
      return true;
    }
  }
  // XR 2D HUD (bottom overlay) with primary simulation controls for Quest3
  let _xrHudCreated = false;
  function ensureXRHUD() {
    if (_xrHudCreated) return;
    try {
      const dbgInfo = {
        hasScene: !!scene,
        hasBABYLON: typeof BABYLON !== 'undefined',
        hasGUI: !!(typeof BABYLON !== 'undefined' && BABYLON.GUI),
        hasADT: !!(
          typeof BABYLON !== 'undefined' &&
          BABYLON.GUI &&
          BABYLON.GUI.AdvancedDynamicTexture
        ),
        sessionMode: sessionMode(),
        existingHUD: !!window.__XR_HUD_ADT,
      };
      console.log('[XR][HUD] ensureXRHUD pre-check', dbgInfo);
      // Delegate to extracted builder (top/bottom world-space bars)
      if (
        !scene ||
        typeof BABYLON === 'undefined' ||
        !BABYLON.GUI ||
        !BABYLON.GUI.AdvancedDynamicTexture
      ) {
        console.warn('[XR][HUD] prerequisites missing');
        return;
      }
      ensureXRHudBars({
        scene,
        getViewer: () => {
          try {
            return window._viewer || window.viewerApi;
          } catch {}
          return null;
        },
      });
      _xrHudCreated = !!window.__XR_HUD_FALLBACK;
    } catch (e) {
      console.warn('[XR][HUD] failed', e);
    }
  }
  // spawnWorldPanelFallback now provided by xr-hud-bars.js via ensure/force helpers

  // Global pointer ray debug (controller forward ray pick) to verify intersection with HUD
  try {
    if (scene && !scene.__xrGlobalRayDebug) {
      scene.__xrGlobalRayDebug = true;
      scene.onBeforeRenderObservable.add(() => {
        try {
          const xr =
            window._viewer && window._viewer.vr && window._viewer.vr.getXR
              ? window._viewer.vr.getXR()
              : null;
          const helper = xrHelper || xr;
          const input = helper?.input;
          if (!input) return;
          input.controllers?.forEach((c) => {
            const n = c.pointer || c.grip || c.pointer?.parent;
            if (!n?.getDirection || !n.getAbsolutePosition) return;
            const origin = n.getAbsolutePosition();
            const dir = n.getDirection(BABYLON.Vector3.Forward()).normalize();
            const ray = new BABYLON.Ray(origin, dir, 10);
            const pr = scene.pickWithRay(ray, (m) => m && (!m.isPickable ? false : true));
            if (
              window.XR_HUD_DEBUG &&
              pr?.hit &&
              pr.pickedMesh &&
              pr.pickedMesh.name &&
              /xrHudPanel|btn/i.test(pr.pickedMesh.name)
            ) {
              console.log('[XR][HUD][RayHit]', pr.pickedMesh.name, { dist: pr.distance });
            }
          });
        } catch {}
      });
    }
  } catch {}
  async function init() {
    const nodeLike = typeof window === 'undefined' || typeof navigator === 'undefined';
    const lacksRenderLoop = !scene?.onBeforeRenderObservable;
    if (nodeLike || lacksRenderLoop) {
      // Fabricate shim environment for unit tests (no real XR / render loop)
      supported = true;
      xrHelper = xrHelper || { input: { controllers: [] }, baseExperience: { sessionManager: {} } };

      return { supported: true, shim: true };
    }
    try {
      if (!engineRef && scene?.getEngine) engineRef = scene.getEngine();
      console.log('[VR][init] calling setupVR (createDefaultXRExperienceAsync)...');
      xrHelper = await setupVR(engineRef, scene, picking);
      supported = !!xrHelper;
      console.log('[VR][init] setupVR result', {
        supported,
        hasHelper: !!xrHelper,
        hasBase: !!xrHelper?.baseExperience,
        hasSM: !!xrHelper?.baseExperience?.sessionManager,
      });
      // If helper missing or unsupported features, fall back to shim controllers for tests to still drive laser logic
      // If unsupported we do not fabricate controllers here (legacy behavior retained)
    } catch (e) {
      // Fall back to shim on failure
      console.warn('[VR][init] setupVR failed; falling back to shim', e?.message || e);
      supported = true;
      xrHelper = { input: { controllers: [] }, baseExperience: { sessionManager: {} } };
    }
    return { supported };
  }
  // Legacy API names kept: enterVR/exitVR
  function legacyEnterVR() {
    enterVR();
  }
  function legacyExitVR() {
    exitXR();
  }
  // Expose diagnostic helper
  function debugInfo() {
    return {
      supported,
      hasXRHelper: !!xrHelper,
      sessionMode: sessionMode(),
      environmentBlendMode: environmentBlend(),
      sceneAutoClear: scene?.autoClear,
      clearColor: scene?.clearColor && {
        r: scene.clearColor.r,
        g: scene.clearColor.g,
        b: scene.clearColor.b,
        a: scene.clearColor.a,
      },
      canvasBg: (scene?.getEngine?.()?.getRenderingCanvas?.() || {}).style?.background,
      bodyBg:
        (typeof document !== 'undefined' && document.body && document.body.style.background) ||
        null,
    };
  }
  try {
    if (typeof window !== 'undefined') window.vrDebugInfo = debugInfo;
  } catch {}
  function forceWorldHUD() {
    return forceXRHudBars({
      scene,
      getViewer: () => {
        try {
          return window._viewer || window.viewerApi;
        } catch {}
        return null;
      },
    });
  }
  function getLaserInfo() {
    try {
      return typeof laserMgr !== 'undefined' && laserMgr?.getInfo ? laserMgr.getInfo() : [];
    } catch {
      return [];
    }
  }
  return {
    init,
    isSupported,
    enterVR: legacyEnterVR,
    exitVR: legacyExitVR,
    enterAR,
    exitXR,
    switchXR,
    controllers,
    controllerStates,
    setTransparent,
    debugInfo,
    ensureXRHUD,
    forceWorldHUD,
    getLaserInfo,
    _debugHUD: () => ({
      created: _xrHudCreated,
      pollAttempts: _hudPollAttempts,
      session: sessionMode(),
      hasADT: !!window.__XR_HUD_ADT,
      hasFallback: !!window.__XR_HUD_FALLBACK,
    }),
  };
}

export default { setupVR, setupVRFeatures, createVRSupport };
