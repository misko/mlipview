// Core logic for XR HUD controls (VR/AR) with no Babylon dependencies.
// Provides event handlers and label/selection state for Simulation (Relax/MD/Off),
// Forces toggle, and Reset. Designed to be wired to either DOM buttons or Babylon GUI buttons.

/**
 * @typedef {{
 *   getViewer: () => any,
 *   reloadPage?: () => void,
 *   onStateChange?: (state: ReturnType<typeof getState>) => void,
 * }} XRControlsOptions
 */

function getState(viewer) {
  const running = (() => {
    try {
      return viewer?.getMetrics?.().running || null;
    } catch {
      return null;
    }
  })();
  const showForces = !!(viewer && viewer.state && viewer.state.showForces);
  return { running, showForces };
}

/**
 * Build XR controls model: event handlers and label helpers.
 * Consumers should call model.refresh() to sync selection from viewer state when needed.
 *
 * @param {XRControlsOptions} opts
 */
export function buildXRControlsModel(opts) {
  const getViewer = opts?.getViewer || (() => null);
  const reload =
    opts?.reloadPage ||
    (() => {
      try {
        if (typeof location !== 'undefined' && typeof location.reload === 'function')
          location.reload();
      } catch {}
    });
  const notify = (s) => {
    try {
      opts?.onStateChange && opts.onStateChange(s);
    } catch {}
  };

  function setSimulation(mode) {
    const v = getViewer();
    if (!v) return false;
    try {
      const current = (() => {
        try {
          return v?.getMetrics?.().running || null;
        } catch {
          return null;
        }
      })();
      // Toggle semantics: pressing active relax/md switches to off; pressing off while off stays off
      if (mode === 'relax') {
        if (current === 'relax') {
          v.stopSimulation?.();
        } else {
          v.stopSimulation?.();
          v.startRelaxContinuous?.({});
        }
      } else if (mode === 'md') {
        if (current === 'md') {
          v.stopSimulation?.();
        } else {
          v.stopSimulation?.();
          v.startMDContinuous?.({});
        }
      } else {
        // off
        v.stopSimulation?.();
      }
      notify(getState(v));
      return true;
    } catch {
      return false;
    }
  }
  function simSelection() {
    const v = getViewer();
    const r = getState(v).running;
    return { relax: r === 'relax', md: r === 'md', off: !r };
  }
  function toggleForces() {
    const v = getViewer();
    if (!v) return false;
    try {
      if (typeof v.setForceVectorsEnabled === 'function') {
        const cur = !!(v.state && v.state.showForces);
        v.setForceVectorsEnabled(!cur);
      } else {
        v.state?.toggleForceVectorsVisibility?.();
      }
      notify(getState(v));
      return true;
    } catch {
      return false;
    }
  }
  function isForcesOn() {
    const v = getViewer();
    return !!getState(v).showForces;
  }
  function isPBCOn() {
    const v = getViewer();
    try {
      return !!(v && v.state && v.state.showCell && v.state.cell && v.state.cell.enabled);
    } catch {
      return false;
    }
  }
  function togglePBC() {
    const v = getViewer();
    if (!v) return false;
    try {
      const st = v.state;
      if (!st) return false;
      const fn = st.toggleCellVisibilityEnhanced || st.toggleCellVisibility;
      if (typeof fn === 'function') fn.call(st);
      const wantGhosts = !!st.showCell;
      if (!!st.showGhostCells !== wantGhosts && typeof st.toggleGhostCells === 'function')
        st.toggleGhostCells();
      try {
        if (typeof v.recomputeBonds === 'function') v.recomputeBonds('xrTogglePBC');
      } catch {}
      try {
        Promise.resolve().then(() => {
          try {
            v.view?.rebuildGhosts?.();
          } catch {}
        });
      } catch {}
      notify(getState(v));
      return true;
    } catch {
      return false;
    }
  }
  // Retain isPBCOn helper (used for button highlight); label now static 'PBC'
  function reset() {
    const v = getViewer();
    const when =
      typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now();
    try {
      console.log('[XR][HUD] Reset pressed', { hasViewer: !!v, t: new Date().toISOString() });
    } catch {}
    if (!v) return false;
    try {
      if (v && typeof v.resetToInitialPositions === 'function') {
        const t0 =
          typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now();
        const res = v.resetToInitialPositions();
        Promise.resolve(res)
          .then((ok) => {
            const t1 =
              typeof performance !== 'undefined' && performance.now
                ? performance.now()
                : Date.now();
            try {
              console.log('[XR][Reset] completed', {
                ok: !!ok,
                ms: Math.round((t1 - t0) * 100) / 100,
              });
            } catch {}
          })
          .catch((e) => {
            try {
              console.warn('[XR][Reset] failed', e?.message || e);
            } catch {}
          });
        notify(getState(v));
        return true;
      }
      try {
        console.warn('[XR][Reset] viewer missing resetToInitialPositions(), reloading page');
      } catch {}
      if (typeof location !== 'undefined' && location.reload) location.reload();
      return true;
    } catch (e) {
      try {
        console.warn('[XR][Reset] exception', e?.message || e);
      } catch {}
      return false;
    }
  }
  function refresh() {
    const v = getViewer();
    notify(getState(v));
    return simSelection();
  }
  function ensureDefaultActive() {
    const v = getViewer();
    if (!v) return false;
    try {
      const running = (() => {
        try {
          return v?.getMetrics?.().running || null;
        } catch {
          return null;
        }
      })();
      if (!running) {
        v.stopSimulation?.();
        v.startMDContinuous?.({});
        notify(getState(v));
      }
      return true;
    } catch {
      return false;
    }
  }

  return {
    setSimulation,
    simSelection,
    toggleForces,
    isForcesOn,
    togglePBC,
    isPBCOn,
    reset,
    refresh,
    ensureDefaultActive,
  };
}

export default { buildXRControlsModel };
