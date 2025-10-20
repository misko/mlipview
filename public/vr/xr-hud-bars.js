// XR HUD Bars (world-space panels): encapsulates creation of top and bottom bars
// Provides ensureWorldHUD(scene, options) and forceWorldHUD(scene) helpers.
// This was extracted from setup.js to simplify VR setup concerns.

// Contract
// - Inputs: { scene, getViewer? }
// - Side effects: creates two mesh planes with AdvancedDynamicTexture GUIs that follow the camera
// - Globals for diagnostics/back-compat: window.__XR_HUD_FALLBACK, window.__XR_HUD_FALLBACK_TOP
// - Returns: { plane, tex, followObs, planeTop, texTop, followObsTop, buttons:[...], energy:true }

// Configuration defaults (fractions of half frustum height). Positive fraction means below center
// when used directly (we subtract up*downward). Top row uses negative verticalFrac to invert.
const DEFAULT_TOP_FRAC_MAG = 0.2; // 20% of halfHeight above center
const DEFAULT_BOTTOM_FRAC_MAG = 0.45; // 45% of halfHeight below center (raised from legacy 60%)

function clamp01(v) {
  return isFinite(v) ? Math.max(0, Math.min(1, v)) : 0;
}
function resolveFracMag(kind) {
  try {
    if (typeof window !== 'undefined') {
      if (kind === 'top' && window.XR_HUD_TOP_FRAC != null)
        return clamp01(Math.abs(Number(window.XR_HUD_TOP_FRAC)));
      if (kind === 'bottom' && window.XR_HUD_BOTTOM_FRAC != null)
        return clamp01(Math.abs(Number(window.XR_HUD_BOTTOM_FRAC)));
      const qs = (window.location && window.location.search) || '';
      if (kind === 'top') {
        const m = /[?&]xrtopfrac=([0-9.]+)/.exec(qs);
        if (m) return clamp01(Math.abs(parseFloat(m[1])));
      }
      if (kind === 'bottom') {
        const m = /[?&]xrbotfrac=([0-9.]+)/.exec(qs);
        if (m) return clamp01(Math.abs(parseFloat(m[1])));
      }
    }
  } catch {}
  return kind === 'top' ? DEFAULT_TOP_FRAC_MAG : DEFAULT_BOTTOM_FRAC_MAG;
}

export function ensureWorldHUD({
  scene,
  getViewer,
  topEnergyDistanceMult = 2,
  topEnergyScaleMult = 2,
  topEnergyIsPickable = false,
  maxHUDTextureDim,
} = {}) {
  try {
    if (typeof window !== 'undefined' && window.__XR_HUD_FALLBACK) return window.__XR_HUD_FALLBACK;
    if (
      !scene ||
      typeof BABYLON === 'undefined' ||
      !BABYLON.GUI ||
      !BABYLON.GUI.AdvancedDynamicTexture
    ) {
      console.warn('[XR][HUD] world HUD prerequisites missing');
      return null;
    }
    const cam = scene.activeCamera;
    if (!cam) return null;

    const dist = 1.1; // base distance for bottom row
    // Sanitize user-provided multipliers
    topEnergyDistanceMult =
      Number.isFinite(topEnergyDistanceMult) && topEnergyDistanceMult > 0
        ? topEnergyDistanceMult
        : 2;
    topEnergyScaleMult =
      Number.isFinite(topEnergyScaleMult) && topEnergyScaleMult > 0 ? topEnergyScaleMult : 2;
    const topDist = dist * topEnergyDistanceMult;
    const fov = cam.fov || Math.PI / 2; // vertical fov
    const halfHeight = Math.tan(fov / 2) * dist;

    // Debug flags helpers
    const wpEnergyDebug = () =>
      typeof window !== 'undefined' &&
      (window.XR_ENERGY_DEBUG ||
        window.XR_HUD_DEBUG ||
        /[?&](xrenergydebug|xrhuddebug)=1/.test(window.location.search));
    const hudDebug = () =>
      typeof window !== 'undefined' &&
      (window.XR_HUD_DEBUG || /[?&]xrhuddebug=1/.test(window.location.search));
    function logBtn(phase, name, extra) {
      if (!hudDebug()) return;
      try {
        console.log('[XR][HUD][PANEL][BTN]', phase, name, extra || '');
      } catch {}
    }

    function btn(stack, label, cb) {
      const id = 'w' + label + Math.random().toString(36).slice(2, 7);
      const b = BABYLON.GUI.Button.CreateSimpleButton(id, label);
      b.width = '240px';
      b.height = '260px';
      b.color = '#d8e6f3';
      b.thickness = 0;
      b.background = 'rgba(60,70,82,0.85)';
      b.fontSize = 72;
      b.onPointerDownObservable.add(() => logBtn('down', label));
      b.onPointerUpObservable.add(() => {
        logBtn('up', label);
        try {
          cb();
          logBtn('invoke', label, 'OK');
        } catch (e) {
          logBtn('invokeErr', label, e?.message || e);
        }
      });
      b.onPointerEnterObservable.add(() => {
        b.background = 'rgba(90,140,190,0.9)';
        logBtn('enter', label);
      });
      b.onPointerOutObservable.add(() => {
        b.background = 'rgba(60,70,82,0.85)';
        logBtn('leave', label);
      });
      stack.addControl(b);
      return b;
    }

    // Resolve maximum HUD texture dimension (engine cap or provided). We fetch once per call.
    let resolvedMaxHUDTex = 0;
    try {
      if (Number.isFinite(maxHUDTextureDim) && maxHUDTextureDim > 0) {
        resolvedMaxHUDTex = maxHUDTextureDim | 0;
      } else {
        const eng = scene.getEngine?.();
        const cap = eng?.getCaps?.().maxTextureSize;
        resolvedMaxHUDTex = Number.isFinite(cap) && cap > 0 ? cap : 4096; // conservative default for standalone headsets
      }
    } catch {
      resolvedMaxHUDTex = 4096;
    }

    function buildWorldRow({
      name = 'xrHudPanel',
      verticalFrac = 0.6,
      top = false,
      includeEnergy = true,
      energyScaleX = 1,
      energyScaleY = 1,
      planeWidth,
      planeHeight,
    }) {
      // Determine scale from explicit plane size when provided; otherwise fall back to energy scale
      const basePW = 1.35,
        basePH = 0.38;
      const hasPW = typeof planeWidth === 'number' && planeWidth > 0;
      const hasPH = typeof planeHeight === 'number' && planeHeight > 0;
      const sX = hasPW ? planeWidth / basePW : includeEnergy ? Math.max(1, energyScaleX || 1) : 1;
      const sY = hasPH ? planeHeight / basePH : includeEnergy ? Math.max(1, energyScaleY || 1) : 1;
      const pw = hasPW ? planeWidth : basePW * sX;
      const ph = hasPH ? planeHeight : basePH * sY;
      const plane = BABYLON.MeshBuilder.CreatePlane(name, { width: pw, height: ph }, scene);
      const forward = cam.getDirection(BABYLON.Axis.Z).normalize();
      const up = cam.getDirection
        ? cam.getDirection(BABYLON.Axis.Y).normalize()
        : BABYLON.Vector3.Up();
      const downward = halfHeight * verticalFrac;
      const basePos = cam.position.add(forward.scale(dist)).subtract(up.scale(downward));
      plane.position.copyFrom(basePos);
      plane.billboardMode = BABYLON.Mesh.BILLBOARDMODE_NONE;
      plane.isPickable = true;
      try {
        plane.metadata = { hudPanel: true, row: top ? 'top' : 'bottom' };
      } catch {}

      let texW = Math.round(1300 * sX);
      let texH = Math.round(340 * sY);
      // Clamp oversized GUI textures to avoid black planes on devices with lower max texture size (e.g. mobile XR headsets)
      let rawW = texW,
        rawH = texH;
      let clamped = false;
      try {
        const maxDim = resolvedMaxHUDTex;
        if (maxDim && (texW > maxDim || texH > maxDim)) {
          const scale = Math.min(maxDim / texW, maxDim / texH);
          texW = Math.max(16, Math.floor(texW * scale));
          texH = Math.max(16, Math.floor(texH * scale));
          clamped = true;
        }
      } catch {}
      const tex = BABYLON.GUI.AdvancedDynamicTexture.CreateForMesh(plane, texW, texH, false);
      try {
        tex._rootContainer?.children?.forEach((c) => {
          c.metadata = c.metadata || {};
          c.metadata.hudRoot = true;
        });
      } catch {}
      try {
        plane.metadata = plane.metadata || {};
        plane.metadata.hudTexSize = {
          w: texW,
          h: texH,
          rawW,
          rawH,
          clamped,
          max: resolvedMaxHUDTex,
        };
      } catch {}

      const bg = new BABYLON.GUI.Rectangle();
      bg.thickness = 0;
      bg.cornerRadius = 20;
      bg.background = top ? 'transparent' : 'rgba(25,32,42,0.78)';
      tex.addControl(bg);

      const stack = new BABYLON.GUI.StackPanel();
      stack.isVertical = false;
      stack.height = 300 * sY + 'px';
      stack.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
      bg.addControl(stack);

      // Energy panel (line plot via offscreen canvas -> Image source)
      if (includeEnergy) {
        const energyBox = new BABYLON.GUI.StackPanel();
        const __eScaleX = Math.max(1, energyScaleX || 1);
        const __eScaleY = Math.max(1, energyScaleY || 1);
        energyBox.isVertical = true;
        energyBox.width = 360 * __eScaleX + 'px';
        energyBox.height = 260 * __eScaleY + 'px';
        energyBox.paddingRight = '10px';
        let energyText = null;
        try {
          const w = 300 * __eScaleX,
            h = 150 * __eScaleY;
          const canvas2 = document.createElement('canvas');
          canvas2.width = w;
          canvas2.height = h;
          if (wpEnergyDebug())
            console.log('[XR][HUD][ENERGY][WP] canvas created', {
              w,
              h,
              row: top ? 'top' : 'bottom',
            });
          let url02 = '';
          try {
            url02 = canvas2.toDataURL('image/png');
            if (wpEnergyDebug())
              console.log('[XR][HUD][ENERGY][WP] primed dataURL', {
                len: (url02 && url02.length) || 0,
              });
          } catch (e) {
            if (wpEnergyDebug())
              console.warn('[XR][HUD][ENERGY][WP] toDataURL failed (init)', e?.message || e);
          }
          const img = new BABYLON.GUI.Image(
            'xrWorldEnergyImg_' + (top ? 'top' : 'bottom'),
            url02 || undefined
          );
          img.width = w + 'px';
          img.height = h + 'px';
          img.stretch = BABYLON.GUI.Image.STRETCH_UNIFORM;
          energyBox.addControl(img);
          energyText = new BABYLON.GUI.TextBlock(
            'xrWorldEnergyValue_' + (top ? 'top' : 'bottom'),
            'E: —'
          );
          energyText.color = '#d8e6f3';
          energyText.fontSize = 28 * __eScaleY;
          energyText.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_LEFT;
          energyBox.addControl(energyText);
          import('../plot/line-plot-core.js')
            .then((p) => {
              try {
                const creator = p.createLinePlot || (p.default && p.default.createLinePlot);
                if (!creator) {
                  if (wpEnergyDebug()) console.warn('[XR][HUD][ENERGY][WP] missing creator');
                  return;
                }
                let updTick = 0;
                const plot = creator({
                  width: w,
                  height: h,
                  getContext: () => canvas2.getContext('2d'),
                  updateTexture: () => {
                    try {
                      const url = canvas2.toDataURL('image/png');
                      img.source = url;
                      img.markAsDirty?.();
                      updTick++;
                      if (wpEnergyDebug() && (updTick <= 3 || updTick % 60 === 0))
                        console.log('[XR][HUD][ENERGY][WP] updateTexture tick', updTick, {
                          urlLen: (url && url.length) || 0,
                          row: top ? 'top' : 'bottom',
                        });
                    } catch (e) {
                      if (wpEnergyDebug())
                        console.warn('[XR][HUD][ENERGY][WP] updateTexture failed', e?.message || e);
                    }
                  },
                  labels: { x: 'step', y: 'E' },
                  maxPoints: 200,
                });
                const v =
                  typeof getViewer === 'function'
                    ? getViewer()
                    : typeof window !== 'undefined'
                      ? window._viewer || window.viewerApi
                      : null;
                function push() {
                  try {
                    const E = v?.state?.dynamics?.energy;
                    if (Number.isFinite(E)) {
                      plot.addPoint(E);
                      if (energyText) energyText.text = 'E: ' + E.toFixed(3) + ' eV';
                    }
                  } catch {}
                }
                try {
                  v?.state?.bus?.on?.('forcesChanged', push);
                } catch (e) {}
                try {
                  push();
                } catch {}
              } catch (e) {
                if (wpEnergyDebug())
                  console.warn('[XR][HUD][ENERGY][WP] build failed', e?.message || e);
              }
            })
            .catch(() => {});
        } catch {}
        stack.addControl(energyBox);
      }

      return { plane, tex, stack };
    }

    // Controls model shared by both rows
    const gv = () => {
      try {
        return getViewer?.() || window._viewer || window.viewerApi;
      } catch {}
      return null;
    };
    let controlsModel;
    let btnSets = { relax: [], md: [], off: [], pbc: [], forces: [] };
    function syncButtons() {
      try {
        const sel = controlsModel?.simSelection?.() || { relax: false, md: false, off: true };
        const activeBg = 'rgba(40,140,100,0.9)';
        const normalBg = 'rgba(60,70,82,0.85)';
        btnSets.relax.forEach((b) => (b.background = sel.relax ? activeBg : normalBg));
        btnSets.md.forEach((b) => (b.background = sel.md ? activeBg : normalBg));
        btnSets.off.forEach((b) => (b.background = sel.off ? activeBg : normalBg));
        const forcesOn = controlsModel?.isForcesOn?.();
        btnSets.forces.forEach((b) => (b.background = forcesOn ? activeBg : normalBg));
        const pbcOn = controlsModel?.isPBCOn?.();
        btnSets.pbc.forEach((b) => (b.background = pbcOn ? activeBg : normalBg));
        // Static labels now (no On/Off suffix)
        btnSets.forces.forEach((b) => {
          try {
            if (b.textBlock) b.textBlock.text = 'Forces';
            else if (b.children && b.children[0]) b.children[0].text = 'Forces';
          } catch {}
        });
        btnSets.pbc.forEach((b) => {
          try {
            if (b.textBlock) b.textBlock.text = 'PBC';
            else if (b.children && b.children[0]) b.children[0].text = 'PBC';
          } catch {}
        });
      } catch {}
    }
    try {
      import('./xr-controls-core.js')
        .then((p) => {
          try {
            const build = p.buildXRControlsModel || (p.default && p.default.buildXRControlsModel);
            if (build) {
              controlsModel = build({
                getViewer: gv,
                onStateChange: syncButtons,
                reloadPage: () => {
                  try {
                    location.reload();
                  } catch {}
                },
              });
              try {
                controlsModel.refresh && controlsModel.refresh();
              } catch {}
              syncButtons();
              try {
                if (window.__XR_HUD_FALLBACK)
                  window.__XR_HUD_FALLBACK.controlsModel = controlsModel;
              } catch {}
            }
          } catch {}
        })
        .catch(() => {});
    } catch {}

    // Bottom row
    // Widen bottom bar to fit: Relax | MD | Off | [sp] | PBC | [sp] | Forces | [sp] | Reset
    // Reduce planeHeight by ~30% (0.38 -> ~0.266) for sleeker bottom bar
    const bottom = buildWorldRow({
      name: 'xrHudPanel',
      verticalFrac: resolveFracMag('bottom'),
      top: false,
      includeEnergy: false,
      planeWidth: 2.0,
      planeHeight: 0.38 * 0.7,
    });
    const b_relax = btn(bottom.stack, 'Relax', () => {
      controlsModel?.setSimulation('relax');
      syncButtons();
    });
    btnSets.relax.push(b_relax);
    const b_md = btn(bottom.stack, 'MD', () => {
      controlsModel?.setSimulation('md');
      syncButtons();
    });
    btnSets.md.push(b_md);
    const b_off = btn(bottom.stack, 'Off', () => {
      controlsModel?.setSimulation('off');
      syncButtons();
    });
    btnSets.off.push(b_off);
    // Spacer, then PBC toggle, then spacer
    try {
      const sp1 = new BABYLON.GUI.Rectangle();
      sp1.width = '30px';
      sp1.height = '260px';
      sp1.thickness = 0;
      sp1.background = 'transparent';
      bottom.stack.addControl(sp1);
    } catch {}
    const b_pbc = btn(bottom.stack, 'PBC', () => {
      controlsModel?.togglePBC?.();
      syncButtons();
    });
    btnSets.pbc.push(b_pbc);
    try {
      const sp1b = new BABYLON.GUI.Rectangle();
      sp1b.width = '30px';
      sp1b.height = '260px';
      sp1b.thickness = 0;
      sp1b.background = 'transparent';
      bottom.stack.addControl(sp1b);
    } catch {}
    const b_forces = btn(bottom.stack, 'Forces', () => {
      controlsModel?.toggleForces?.();
      syncButtons();
    });
    btnSets.forces.push(b_forces);
    // Spacer between Forces and Reset
    try {
      const sp2 = new BABYLON.GUI.Rectangle();
      sp2.width = '30px';
      sp2.height = '260px';
      sp2.thickness = 0;
      sp2.background = 'transparent';
      bottom.stack.addControl(sp2);
    } catch {}
    const b_reset = btn(bottom.stack, 'Reset', () => {
      controlsModel?.reset?.();
      syncButtons();
    });

    // Top row mirror
    const topRow = buildWorldRow({
      name: 'xrHudPanelTop',
      verticalFrac: -resolveFracMag('top'),
      top: true,
      includeEnergy: true,
      // Scale the energy plot internals proportionally so the graph remains readable at farther distance
      energyScaleX: 3 * topEnergyScaleMult,
      energyScaleY: 2.5 * topEnergyScaleMult,
      planeWidth: 1.35 * 2 * 2 * topEnergyScaleMult,
      planeHeight: 0.38 * 1.5 * 1.5 * topEnergyScaleMult,
    });
    // Make energy panel non-pickable by default so it doesn't block atom/bond selection.
    try {
      topRow.plane.isPickable = !!topEnergyIsPickable;
    } catch {}

    try {
      controlsModel?.ensureDefaultActive?.();
    } catch {}
    syncButtons();

    // Camera-follow behaviors for both rows
    const followBottom = scene.onBeforeRenderObservable.add(() => {
      try {
        const c = scene.activeCamera;
        if (!c) return;
        const mag = resolveFracMag('bottom');
        const hh = Math.tan((c.fov || Math.PI / 2) / 2) * dist;
        const fw = c.getDirection(BABYLON.Axis.Z).normalize();
        const up2 = c.getDirection
          ? c.getDirection(BABYLON.Axis.Y).normalize()
          : BABYLON.Vector3.Up();
        const target = c.position.add(fw.scale(dist)).subtract(up2.scale(hh * mag));
        bottom.plane.position.copyFrom(target);
        const toCam = c.position.subtract(bottom.plane.position);
        toCam.y = 0;
        toCam.normalize();
        const yaw = Math.atan2(toCam.x, toCam.z);
        bottom.plane.rotation.y = yaw + Math.PI;
      } catch {}
    });
    const followTop = scene.onBeforeRenderObservable.add(() => {
      try {
        const c = scene.activeCamera;
        if (!c) return;
        const mag = resolveFracMag('top'); // Keep vertical offset based on base dist so Y band unchanged
        const baseHH = Math.tan((c.fov || Math.PI / 2) / 2) * dist;
        const fw = c.getDirection(BABYLON.Axis.Z).normalize();
        const up2 = c.getDirection
          ? c.getDirection(BABYLON.Axis.Y).normalize()
          : BABYLON.Vector3.Up();
        const target = c.position.add(fw.scale(topDist)).subtract(up2.scale(baseHH * -mag));
        topRow.plane.position.copyFrom(target);
        const toCam = c.position.subtract(topRow.plane.position);
        toCam.y = 0;
        toCam.normalize();
        const yaw = Math.atan2(toCam.x, toCam.z);
        topRow.plane.rotation.y = yaw + Math.PI;
      } catch {}
    });

    console.log('[XR][HUD] world bars created (camera-follow top & bottom)');
    const result = {
      plane: bottom.plane,
      tex: bottom.tex,
      followObs: followBottom,
      planeTop: topRow.plane,
      texTop: topRow.tex,
      followObsTop: followTop,
      buttons: ['Relax', 'MD', 'Off', 'Forces', 'Reset'],
      energy: true,
      getControlsModel: () => controlsModel,
      topEnergyDistanceMult,
      topEnergyScaleMult,
    };
    try {
      window.__XR_HUD_FALLBACK = result;
      window.__XR_HUD_FALLBACK_TOP = { plane: topRow.plane, tex: topRow.tex };
    } catch {}
    return result;
  } catch (e) {
    console.warn('[XR][HUD] world HUD create failed', e);
    return null;
  }
}

export function forceWorldHUD({ scene, getViewer } = {}) {
  if (typeof window !== 'undefined' && window.__XR_HUD_FALLBACK) return window.__XR_HUD_FALLBACK;
  return ensureWorldHUD({ scene, getViewer });
}

export default { ensureWorldHUD, forceWorldHUD };
