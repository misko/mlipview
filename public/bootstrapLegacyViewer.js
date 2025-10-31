import { initNewViewer } from './index.js';
import { loadDefault } from './util/moleculeLoader.js';
import { buildDesktopPanel } from './ui/desktopPanel.js';
import { installTouchControls } from './ui/touchControls.js';
import { initForcesToggle } from './ui/forcesToggle.js';
import { initCellToggle } from './ui/cellToggle.js';
import { showErrorBanner } from './ui/errorBanner.js';

function setStatusText(statusEl, text) {
  if (!statusEl) return;
  statusEl.textContent = text;
}

function bindHandler(element, handlerName, fn, restorers) {
  if (!element) return;
  const prev = element[handlerName];
  element[handlerName] = fn;
  if (restorers) restorers.push(() => {
    element[handlerName] = prev ?? null;
  });
}

export async function bootstrapLegacyViewer(options = {}) {
  const {
    canvas = typeof document !== 'undefined' ? document.getElementById('viewer') : null,
    appRoot = typeof document !== 'undefined' ? document.getElementById('app') : null,
    statusElement = typeof document !== 'undefined' ? document.getElementById('status') : null,
    urlParams = new URLSearchParams(typeof location !== 'undefined' ? location.search || '' : ''),
    useReactUi = false,
  } = options;

  if (!canvas) throw new Error('bootstrapLegacyViewer: viewer canvas not found');
  if (!appRoot) throw new Error('bootstrapLegacyViewer: app root not found');

  const DEBUG = urlParams.get('debug') === '1';
  const inlineCleanup = [];
  const teardownFns = [];

  const viewerApi = await initNewViewer(canvas, { elements: [], positions: [], bonds: [] });

  if (typeof window !== 'undefined') {
    window._viewer = viewerApi;
    try {
      window.viewerApi = viewerApi;
      window.dispatchEvent(new CustomEvent('mlipviewer-ready', { detail: { when: Date.now() } }));
    } catch {}
    window.vrScene = viewerApi.scene;
    window.vrEngine = viewerApi.engine;
  }

  try {
    installTouchControls({ canvas, scene: viewerApi.scene, camera: viewerApi.camera, picking: viewerApi.picking });
  } catch {}

  if (DEBUG) console.log('[debug] viewer initialized');

  if (!useReactUi) {
    try {
      buildDesktopPanel({ attachTo: appRoot });
    } catch (err) {
      console.warn('[bootstrapLegacyViewer] desktop panel setup failed', err);
    }
  }

  const libraryId = urlParams.get('library');
  let loadedLabel = null;
  try {
    if (libraryId) {
      await viewerApi.session.loadFromLibrary(libraryId);
      loadedLabel = `library:${libraryId}`;
      if (DEBUG) {
        console.log('[debug] library entry loaded', libraryId, 'atoms=', viewerApi.state.positions.length);
      }
    } else {
      const { file } = await loadDefault(viewerApi);
      try { viewerApi.baselineEnergy?.(); } catch {}
      loadedLabel = file;
      try {
        const sel = typeof document !== 'undefined' ? document.getElementById('moleculeSelect') : null;
        if (sel) sel.value = file;
      } catch {}
      if (DEBUG) {
        console.log('[debug] default molecule loaded', file, 'atoms=', viewerApi.state.positions.length, 'energySeriesLen=1');
      }
    }

    if (!loadedLabel) throw new Error('Unknown load result');
    setStatusText(statusElement, `Loaded ${loadedLabel}`);
  } catch (primaryError) {
    console.error('[loader] load failed', primaryError);
    if (libraryId) {
      try {
        console.warn('[loader] falling back to default molecule after library failure');
        const { file } = await loadDefault(viewerApi);
        try { viewerApi.baselineEnergy?.(); } catch {}
        loadedLabel = file;
        setStatusText(statusElement, `Loaded ${file}`);
        try {
          showErrorBanner(`Library "${libraryId}" failed to load: ${primaryError?.message || primaryError} — reverted to ${file}`);
        } catch {}
      } catch (fallbackError) {
        setStatusText(statusElement, 'Load failed');
        try { showErrorBanner('Load failed: ' + (fallbackError?.message || fallbackError)); } catch {}
      }
    } else {
      setStatusText(statusElement, 'Load failed');
      try { showErrorBanner('Load failed: ' + (primaryError?.message || primaryError)); } catch {}
    }
  } finally {
    if (typeof window !== 'undefined') {
      window.__MLIP_DEFAULT_LOADED = true;
    }
  }

  const timers = [];

  if (!useReactUi) {
    const relaxStepBtn = typeof document !== 'undefined' ? document.getElementById('btnRelax') : null;
    const relaxRunBtn = typeof document !== 'undefined' ? document.getElementById('btnRelaxRun') : null;
    const mdStepBtn = typeof document !== 'undefined' ? document.getElementById('btnMD') : null;
    const mdRunBtn = typeof document !== 'undefined' ? document.getElementById('btnMDRun') : null;
    const forceProviderSel = typeof document !== 'undefined' ? document.getElementById('forceProviderSel') : null;
    const xrSel = typeof document !== 'undefined' ? document.getElementById('xrModeSelect') : null;

    const isRelaxRunning = () => {
      try { return !!(viewerApi && viewerApi.getMetrics().running === 'relax'); } catch { return false; }
    };
    const isMdRunning = () => {
      try { return !!(viewerApi && viewerApi.getMetrics().running === 'md'); } catch { return false; }
    };

    const applyRunningDisable = () => {
      const runningKind = viewerApi && viewerApi.getMetrics().running;
      const disabled = !!runningKind;
      if (relaxStepBtn) relaxStepBtn.disabled = disabled;
      if (mdStepBtn) mdStepBtn.disabled = disabled;
      if (relaxRunBtn) relaxRunBtn.disabled = false;
      if (mdRunBtn) mdRunBtn.disabled = false;
      if (runningKind === 'relax') {
        if (mdRunBtn) mdRunBtn.disabled = true;
        if (mdStepBtn) mdStepBtn.disabled = true;
      } else if (runningKind === 'md') {
        if (relaxRunBtn) relaxRunBtn.disabled = true;
        if (relaxStepBtn) relaxStepBtn.disabled = true;
      }
    };

    bindHandler(relaxStepBtn, 'onclick', async () => {
      if (!viewerApi || isRelaxRunning() || isMdRunning()) return;
      await viewerApi.relaxStep();
      setStatusText(statusElement, 'Relax step');
    }, inlineCleanup);

    bindHandler(relaxRunBtn, 'onclick', async () => {
      if (!viewerApi) return;
      if (!isRelaxRunning()) {
        if (relaxRunBtn) relaxRunBtn.textContent = 'stop';
        setStatusText(statusElement, 'Relax running');
        applyRunningDisable();
        viewerApi.startRelaxContinuous({ maxSteps: 500 }).then((r) => {
          if (relaxRunBtn) relaxRunBtn.textContent = 'run';
          applyRunningDisable();
          setStatusText(statusElement, r?.converged ? 'Relax converged' : 'Relax stopped');
        });
      } else {
        viewerApi.stopSimulation();
        if (relaxRunBtn) relaxRunBtn.textContent = 'run';
        applyRunningDisable();
        setStatusText(statusElement, 'Relax stopped');
      }
    }, inlineCleanup);

    bindHandler(mdStepBtn, 'onclick', () => {
      if (!viewerApi || isMdRunning() || isRelaxRunning()) return;
      const T = typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null
        ? window.__MLIP_TARGET_TEMPERATURE
        : 1500;
      viewerApi.mdStep({ temperature: T });
      const inst = viewerApi.state?.dynamics?.temperature;
      const instLabel = typeof inst === 'number' ? inst.toFixed(1) : 'n/a';
      setStatusText(statusElement, `MD step tgt=${T}K inst=${instLabel}`);
    }, inlineCleanup);

    bindHandler(mdRunBtn, 'onclick', () => {
      if (!viewerApi) return;
      if (!isMdRunning()) {
        const T = typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null
          ? window.__MLIP_TARGET_TEMPERATURE
          : 1500;
        if (mdRunBtn) mdRunBtn.textContent = 'stop';
        setStatusText(statusElement, `MD running @${T}K`);
        applyRunningDisable();
        viewerApi.startMDContinuous({ steps: 500, temperature: T }).then(() => {
          if (!isMdRunning()) {
            if (mdRunBtn) mdRunBtn.textContent = 'run';
            applyRunningDisable();
            setStatusText(statusElement, 'MD stopped');
          }
        });
      } else {
        viewerApi.stopSimulation();
        if (mdRunBtn) mdRunBtn.textContent = 'run';
        applyRunningDisable();
        setStatusText(statusElement, 'MD stopped');
      }
    }, inlineCleanup);

    bindHandler(forceProviderSel, 'onchange', (e) => {
      if (!viewerApi) return;
      const target = e?.target;
      const kind = target && target.value ? target.value : undefined;
      if (kind) {
        viewerApi.setForceProvider(kind);
        setStatusText(statusElement, `Provider ${kind}`);
      }
    }, inlineCleanup);

    const metricsInterval = setInterval(() => {
      if (!viewerApi) return;
      const metrics = viewerApi.getMetrics();
      applyRunningDisable();
      if (!metrics.running) return;
      const parts = [];
      if (metrics.energy != null) parts.push(`E=${metrics.energy.toFixed(3)}`);
      if (metrics.maxForce != null) parts.push(`Fmax=${metrics.maxForce.toExponential(2)}`);
      if (metrics.maxStress != null) parts.push(`Smax=${metrics.maxStress.toExponential(2)}`);
      setStatusText(statusElement, `${metrics.running || 'idle'} ${parts.join(' ')}`);
    }, 500);
    teardownFns.push(() => clearInterval(metricsInterval));
    timers.push(metricsInterval);

    initCellToggle({ getViewer: () => viewerApi });
    initForcesToggle({ getViewer: () => viewerApi });

    bindHandler(xrSel, 'onchange', async () => {
      if (!viewerApi || !xrSel) return;
      const v = xrSel.value;
      if (v === 'none') await viewerApi.vr.switchXR('none');
      else if (v === 'vr') await viewerApi.vr.switchXR('vr');
      else if (v === 'ar') await viewerApi.vr.switchXR('ar');
    }, inlineCleanup);

    const xrPollInterval = setInterval(() => {
      if (!viewerApi || !xrSel) return;
      try {
        const info = viewerApi.vr.debugInfo();
        let mode = 'none';
        if (info.sessionMode.includes('ar')) mode = 'ar';
        else if (info.sessionMode.includes('vr') || (info.hasXRHelper && info.sessionMode !== 'unknown')) mode = 'vr';
        if (xrSel.value !== mode) xrSel.value = mode;
      } catch {}
    }, 2000);
    teardownFns.push(() => clearInterval(xrPollInterval));
    timers.push(xrPollInterval);
  }

  return {
    viewerApi,
    async shutdown() {
      for (const restore of inlineCleanup.splice(0)) {
        try { restore(); } catch {}
      }
      for (const timer of timers.splice(0)) {
        try { clearInterval(timer); } catch {}
      }
      for (const fn of teardownFns.splice(0)) {
        try { fn(); } catch {}
      }
      try { await viewerApi?.shutdown?.(); } catch {}
    },
  };
}
