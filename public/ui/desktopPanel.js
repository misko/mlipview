// Builds the left-side desktop UI panel with collapsible sections.
// Sections: Live Metrics, Simulation, System, Rendering, XR
// Returns { panelEl, containers: { live, simulation, system, rendering, xr } }

import { initTemperatureSlider } from './temperatureSlider.js';
import { initFrictionSlider } from './frictionSlider.js';
import { installMoleculeSelector } from './moleculeSelect.js';

function createSection(id, title, { defaultOpen = false } = {}) {
  const section = document.createElement('div');
  section.className = 'panel-section';
  section.id = id;

  const header = document.createElement('button');
  header.className = 'panel-header';
  header.type = 'button';
  header.setAttribute('aria-expanded', defaultOpen ? 'true' : 'false');
  header.textContent = title;

  const content = document.createElement('div');
  content.className = 'panel-content';
  if (!defaultOpen) content.setAttribute('data-collapsed', 'true');

  header.addEventListener('click', () => {
    const collapsed = content.getAttribute('data-collapsed') === 'true';
    if (collapsed) {
      content.removeAttribute('data-collapsed');
      header.setAttribute('aria-expanded', 'true');
    } else {
      content.setAttribute('data-collapsed', 'true');
      header.setAttribute('aria-expanded', 'false');
    }
  });

  section.appendChild(header);
  section.appendChild(content);
  return { section, content };
}

function makeToggle({ id, labelOn, labelOff, title }) {
  const btn = document.createElement('button');
  btn.id = id;
  btn.type = 'button';
  btn.className = 'toggle';
  btn.setAttribute('role', 'switch');
  btn.setAttribute('aria-checked', 'false');
  btn.setAttribute('data-on', 'false');
  btn.title = title || '';
  btn.textContent = labelOff;

  btn.addEventListener('click', () => {
    const on = btn.getAttribute('data-on') === 'true';
    const next = !on;
    btn.setAttribute('data-on', String(next));
    btn.setAttribute('aria-checked', String(next));
    btn.textContent = next ? labelOn : labelOff;
  });

  return btn;
}

export function buildDesktopPanel({ attachTo } = {}) {
  const host = attachTo || document.getElementById('app') || document.body;

  // Root panel
  const panel = document.createElement('div');
  panel.className = 'hud';
  panel.id = 'controlPanel';
  Object.assign(panel.style, {
    position: 'absolute',
    left: '12px',
    top: '12px',
    maxHeight: 'calc(100% - 24px)',
    width: '300px',
    overflow: 'auto',
    background: '#121212cc', // dark base w/ a touch of translucency
    border: '1px solid rgba(255,255,255,0.07)',
    borderRadius: '10px',
    padding: '8px 10px',
    backdropFilter: 'blur(4px)',
  });

  // Inject CSS once
  if (!document.getElementById('desktopPanelStyles')) {
    const style = document.createElement('style');
    style.id = 'desktopPanelStyles';
    style.textContent = `
      /* Palette & base */
      :root {
        --bg-base: #121212;
        --bg-panel: #1E1E1E;
        --bg-active: #202225;
        --border: #2A2A2A;
        --text: #E5E7EB;
        --muted: #9CA3AF;
        --disabled: #6B7280;
        --accent: #3B82F6;
        --accent-hover: #60A5FA;
      }
      #controlPanel { color: var(--text); font: 13px/1.4 system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, "Helvetica Neue", Arial, sans-serif; }

      .panel-section { margin-bottom: 12px; background: var(--bg-panel); border: 1px solid var(--border); border-radius: 10px; overflow: clip; }
      .panel-header {
        width: 100%;
        text-align: left;
        padding: 8px 10px;
        background: var(--bg-active);
        border: 0;
        color: var(--text);
        font-weight: 600;
        letter-spacing: .02em;
        cursor: pointer;
        transition: background .15s ease;
      }
      .panel-header:hover { background: #25272b; }
      .panel-content { padding: 10px; }
      .panel-content[data-collapsed="true"] { display: none; }

      .panel-subtitle { font-size: 12px; color: var(--muted); text-transform: uppercase; margin: 8px 0 4px; letter-spacing: .06em; }
      .row { display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }
      .row + .row { margin-top: 8px; }
      .block { margin-top: 8px; }
      .form-label { font-size: 12px; color: var(--muted); margin-right: 6px; }
      .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; }

      /* Toggles */
      .toggle {
        appearance: none;
        background: #2A2A2A;
        color: var(--text);
        border: 1px solid var(--border);
        border-radius: 999px;
        padding: 6px 10px;
        font-size: 12px;
        cursor: pointer;
        transition: all .15s ease;
      }
      .toggle:hover { filter: brightness(1.05); }
      .toggle[data-on="true"] { background: var(--accent); border-color: var(--accent); }
      .toggle[disabled] { opacity: .5; cursor: not-allowed; }

      /* Small buttons */
      .btn {
        background: #2A2A2A;
        color: var(--text);
        border: 1px solid var(--border);
        border-radius: 8px;
        padding: 6px 10px;
        font-size: 12px;
        cursor: pointer;
        transition: all .15s ease;
      }
      .btn:hover { border-color: var(--accent); }

      /* Inputs */
      input[type="text"], select {
        background: #171717;
        border: 1px solid var(--border);
        color: var(--text);
        border-radius: 8px;
        padding: 6px 8px;
        font-size: 12px;
        width: 100%;
      }
      input[type="text"]::placeholder { color: #8b8e94; }
      .input-row { display: grid; grid-template-columns: 1fr auto; gap: 8px; }

      /* Energy plot box */
      #energyPlot {
        margin-top: 6px;
        width: 100%;
        height: 90px;
        background: rgba(255,255,255,0.04);
        position: relative;
        overflow: hidden;
        border: 1px solid rgba(255,255,255,0.08);
        border-radius: 8px;
      }
      #energyLabel { position: absolute; left: 6px; top: 4px; font-size: 10px; color: var(--muted); }

      /* Small helpers */
      .stat { font-size: 12px; color: var(--muted); }
      .value { font-size: 12px; }
      .divider { border-top: 1px solid var(--border); margin: 8px 0; }
    `;
    document.head.appendChild(style);
  }

  // Live Metrics (open)
  const live = createSection('section-live-stats', 'Live Metrics', { defaultOpen: true });
  {
    const statsRow = document.createElement('div'); statsRow.className = 'row';
    const status = document.createElement('span'); status.id = 'status'; status.className = 'stat'; status.textContent = 'Ready';
    const instTemp = document.createElement('span'); instTemp.id = 'instTemp'; instTemp.className = 'value mono'; instTemp.textContent = 'T: — K';
    const rps = document.createElement('span'); rps.id = 'rpsLabel'; rps.className = 'value mono'; rps.textContent = 'RPS: —';
    statsRow.append(status, instTemp, rps);
    live.content.appendChild(statsRow);

    const plot = document.createElement('div');
    plot.id = 'energyPlot';
    const canvas = document.createElement('canvas');
    canvas.id = 'energyCanvas';
    canvas.style.position = 'absolute';
    canvas.style.left = '0'; canvas.style.top = '0';
    const label = document.createElement('div'); label.id = 'energyLabel'; label.textContent = 'Energy (eV)';
    plot.append(canvas, label);
    live.content.appendChild(plot);
  }

  // Simulation (open) — Relaxation & MD as toggles
  const sim = createSection('section-simulation', 'Simulation', { defaultOpen: true });
  {
    // Relaxation toggle (mutually exclusive with MD)
    const relaxTitle = document.createElement('div'); relaxTitle.className = 'panel-subtitle'; relaxTitle.textContent = 'Relaxation';
    const relaxRow = document.createElement('div'); relaxRow.className = 'row';
    const relaxToggle = makeToggle({
      id: 'toggleRelax',
      labelOn: 'Relaxation: On',
      labelOff: 'Relaxation: Off',
      title: 'Toggle geometry optimization',
    });
    // MD toggle
    const mdTitle = document.createElement('div'); mdTitle.className = 'panel-subtitle'; mdTitle.textContent = 'Molecular Dynamics';
    const mdRow = document.createElement('div'); mdRow.className = 'row';
    const mdToggle = makeToggle({
      id: 'toggleMD',
      labelOn: 'Molecular Dynamics: On',
      labelOff: 'Molecular Dynamics: Off',
      title: 'Toggle MD',
    });

    // Mutual exclusion: turning one on turns the other off
    const syncExclusion = (src, dst) => {
      src.addEventListener('click', () => {
        const on = src.getAttribute('data-on') === 'true';
        if (on) {
          dst.setAttribute('data-on', 'false');
          dst.setAttribute('aria-checked', 'false');
          dst.textContent = dst.id === 'toggleMD' ? 'Molecular Dynamics: Off' : 'Relaxation: Off';
        }
      });
    };
    syncExclusion(relaxToggle, mdToggle);
    syncExclusion(mdToggle, relaxToggle);

  // Attach behavior to viewer API if available
    function getApi(){ try { return window.viewerApi || window._viewer; } catch { return null; } }
    function setBtnState(btn, on){ btn.setAttribute('data-on', String(!!on)); btn.setAttribute('aria-checked', String(!!on)); btn.textContent = on ? (btn.id==='toggleMD'?'Molecular Dynamics: On':'Relaxation: On') : (btn.id==='toggleMD'?'Molecular Dynamics: Off':'Relaxation: Off'); }
    // Relaxation toggle wiring
    relaxToggle.addEventListener('click', async (e)=>{
      const api = getApi(); if(!api) return;
      const on = relaxToggle.getAttribute('data-on')==='true';
      // If turning on, stop MD and start relax run (single step semantics via run loop)
      if(on){ try { api.stopSimulation(); } catch{} try { await api.startRelaxContinuous({}); } catch{}
      } else { try { api.stopSimulation(); } catch{} }
      // Mirror legacy button text convention for compatibility
      try { const legacy = document.getElementById('btnRelaxRun'); if(legacy) legacy.textContent = on? 'stop':'run'; } catch{}
    });
    // MD toggle wiring
    mdToggle.addEventListener('click', async (e)=>{
      const api = getApi(); if(!api) return;
      const on = mdToggle.getAttribute('data-on')==='true';
      if(on){ try { api.stopSimulation(); } catch{} try { await api.startMDContinuous({}); } catch{}
      } else { try { api.stopSimulation(); } catch{} }
      try {
        const legacy = document.getElementById('btnMDRun'); if(legacy) legacy.textContent = on? 'stop':'run';
        const status = document.getElementById('status'); if(status) status.textContent = on? 'MD running':'MD stopped';
      } catch{}
    });

    relaxRow.appendChild(relaxToggle);
    mdRow.appendChild(mdToggle);

    // Reflect simulation state from API into toggles (handles auto-start MD after load)
    function syncFromApiOnce(){
      try {
        const api = getApi(); if(!api || !api.getMetrics) return;
        const m = api.getMetrics() || {}; const kind = m.running || null;
        if (kind === 'md') {
          // MD running: set MD=On, Relax=Off
          if (mdToggle.getAttribute('data-on') !== 'true') setBtnState(mdToggle, true);
          if (relaxToggle.getAttribute('data-on') !== 'false') setBtnState(relaxToggle, false);
          try { const legacy = document.getElementById('btnMDRun'); if(legacy) legacy.textContent = 'stop'; } catch {}
          try { const status = document.getElementById('status'); if(status) status.textContent = 'MD running'; } catch {}
        } else if (kind === 'relax') {
          // Relax running: set Relax=On, MD=Off
          if (relaxToggle.getAttribute('data-on') !== 'true') setBtnState(relaxToggle, true);
          if (mdToggle.getAttribute('data-on') !== 'false') setBtnState(mdToggle, false);
          try { const legacy = document.getElementById('btnRelaxRun'); if(legacy) legacy.textContent = 'stop'; } catch {}
          try { const status = document.getElementById('status'); if(status) status.textContent = 'Relax running'; } catch {}
        } else {
          // Neither running
          if (mdToggle.getAttribute('data-on') !== 'false') setBtnState(mdToggle, false);
          if (relaxToggle.getAttribute('data-on') !== 'false') setBtnState(relaxToggle, false);
          try { const legacyMD = document.getElementById('btnMDRun'); if(legacyMD) legacyMD.textContent = 'run'; } catch {}
          try { const legacyRX = document.getElementById('btnRelaxRun'); if(legacyRX) legacyRX.textContent = 'run'; } catch {}
          try { const status = document.getElementById('status'); if(status && !/Ready|stopped/i.test(status.textContent||'')) status.textContent = 'Ready'; } catch {}
        }
      } catch {}
    }
    // Initial sync now, then poll briefly to catch auto-start transitions without requiring events
    try { syncFromApiOnce(); } catch {}
    try { 
      let syncPolls = 0; const maxPolls = 60; // ~12s at 200ms
      const t = setInterval(()=>{ syncFromApiOnce(); if(++syncPolls>=maxPolls) clearInterval(t); }, 200);
    } catch {}

    sim.content.append(relaxTitle, relaxRow, mdTitle, mdRow);

    // Sliders (MD-related sliders can always render; your app can disable when MD is Off)
    initTemperatureSlider({ hudEl: sim.content, getViewer: () => null });
    initFrictionSlider({ hudEl: sim.content });
  }

  // System (collapsed)
  const sys = createSection('section-system', 'System', { defaultOpen: false });
  {
    // Preset dropdown (existing helper)
    const presetsRow = document.createElement('div'); presetsRow.className = 'row';
    const label = document.createElement('span'); label.className = 'form-label'; label.textContent = 'Preset:';
    presetsRow.appendChild(label);
    installMoleculeSelector({ hudEl: presetsRow, windowRef: window, documentRef: document });
    sys.content.appendChild(presetsRow);

    // SMILES input (single-line) + Generate button; Enter submits
    const smilesBlock = document.createElement('div'); smilesBlock.className = 'block';
    const smilesLabel = document.createElement('div'); smilesLabel.className = 'form-label'; smilesLabel.textContent = 'SMILES';
    const inputRow = document.createElement('div'); inputRow.className = 'input-row';
    const smilesInput = document.createElement('input'); smilesInput.type = 'text'; smilesInput.id = 'smilesInput'; smilesInput.placeholder = 'Paste or type SMILES here';
    const smilesBtn = document.createElement('button'); smilesBtn.id = 'smilesGoBtn'; smilesBtn.className = 'btn'; smilesBtn.textContent = 'Generate';
    smilesBtn.disabled = true;

    smilesInput.addEventListener('input', () => {
      smilesBtn.disabled = !smilesInput.value.trim();
    });
    smilesInput.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' && !smilesBtn.disabled) {
        smilesBtn.click();
      }
    });

    inputRow.append(smilesInput, smilesBtn);
    smilesBlock.append(smilesLabel, inputRow);
    sys.content.appendChild(smilesBlock);

    // PBC toggle + backend select
    const pbcRow = document.createElement('div'); pbcRow.className = 'row';
    const pbcToggle = makeToggle({
      id: 'togglePBC',
      labelOn: 'PBC: On',
      labelOff: 'PBC: Off',
      title: 'Periodic Boundary Conditions',
    });
    // Wire PBC toggle to viewer state
    pbcToggle.addEventListener('click', ()=>{
      try {
        const api = (window.viewerApi||window._viewer); if(!api) return;
        const st = api.state; if(!st) return;
        const fn = (st.toggleCellVisibilityEnhanced||st.toggleCellVisibility);
        if (typeof fn === 'function') fn.call(st);
        // Keep ghost cells in sync with cell visibility: enable both on, disable both off
        const wantGhosts = !!st.showCell;
        if (!!st.showGhostCells !== wantGhosts && typeof st.toggleGhostCells === 'function') {
          st.toggleGhostCells();
        }
        // Update legacy hidden control to reflect unified state (on if either is shown)
        const legacy = document.getElementById('btnCell');
        if(legacy) legacy.textContent = (st.showCell || st.showGhostCells) ? 'on' : 'off';
        // Optional: status text for clarity
        const status = document.getElementById('status');
        if(status) status.textContent = (st.showCell || st.showGhostCells) ? 'PBC ON' : 'PBC OFF';
      } catch {}
    });
    const backendSel = document.createElement('select'); backendSel.id = 'forceProviderSel';
    backendSel.innerHTML = '<option value="fairchem">UMA (Remote)</option>';
    pbcRow.append(pbcToggle, backendSel);
    sys.content.appendChild(pbcRow);
  }

  // Rendering (collapsed)
  const ren = createSection('section-rendering', 'Rendering', { defaultOpen: false });
  {
    const row1 = document.createElement('div'); row1.className = 'row';
    const forcesT = makeToggle({
      id: 'toggleForces',
      labelOn: 'Show Forces: On',
      labelOff: 'Show Forces: Off',
      title: 'Render per-atom force vectors',
    });
    const energyPlotT = makeToggle({
      id: 'toggleEnergyPlot',
      labelOn: 'Show Energy Plot: On',
      labelOff: 'Show Energy Plot: Off',
      title: 'Toggle energy chart',
    });
    // Default: Energy plot on (matches common default)
    energyPlotT.click(); // set to On visually/semantically
    // Forces wiring: reflect into state toggle and legacy button text
    forcesT.addEventListener('click', ()=>{
      try {
        const api = (window.viewerApi||window._viewer); if(!api) return;
        const on = forcesT.getAttribute('data-on')==='true';
        api.setForceVectorsEnabled(on);
        const legacy = document.getElementById('btnToggleForces'); if(legacy) legacy.textContent = on? 'on':'off';
      } catch{}
    });
    // Energy plot wiring: show/hide plot container and set legacy label
    energyPlotT.addEventListener('click', ()=>{
      const on = energyPlotT.getAttribute('data-on')==='true';
      const box = document.getElementById('energyPlot'); if(box) box.style.display = on? 'block':'none';
      const legacy = document.getElementById('toggleEnergyPlot'); if(legacy) legacy.textContent = on? 'on':'off';
    });
    row1.append(forcesT, energyPlotT);
    ren.content.appendChild(row1);
  }

  // XR (collapsed)
  const xr = createSection('section-xr', 'XR', { defaultOpen: false });
  {
    const row = document.createElement('div'); row.className = 'row';
    const label = document.createElement('span'); label.className = 'form-label'; label.textContent = 'XR mode:';
    const sel = document.createElement('select'); sel.id = 'xrModeSelect';
    sel.innerHTML = '<option value="none">Off</option><option value="vr">VR</option><option value="ar">AR</option>';
    row.append(label, sel);
    xr.content.appendChild(row);
  }

  // Append all sections
  panel.append(live.section, sim.section, sys.section, ren.section, xr.section);
  host.appendChild(panel);

  // Legacy hidden controls for test compatibility (text reflects on/off)
  const legacyIds = [ 'btnRelax','btnRelaxRun','btnMD','btnMDRun','btnToggleForces','toggleEnergyPlot','btnCell' ];
  const legacyBox = document.createElement('div'); legacyBox.style.display='none'; legacyBox.id='legacyControlsHidden';
  for(const id of legacyIds){ const b=document.createElement('button'); b.id=id; b.textContent='run'; legacyBox.appendChild(b); }
  panel.appendChild(legacyBox);

  return {
    panelEl: panel,
    containers: {
      live: live.content,
      simulation: sim.content,
      system: sys.content,
      rendering: ren.content,
      xr: xr.content,
    }
  };
}
