// Builds the left-side desktop UI panel with collapsible sections.
// Sections: Live Metrics, Simulation, System, Rendering, XR
// Returns { panelEl, containers: { live, simulation, system, rendering, xr } }

import { initTemperatureSlider } from './temperatureSlider.js';
import { initFrictionSlider } from './frictionSlider.js';
import {
  installMoleculeSelector,
  buildReloadUrlWithParam,
  base64EncodeUtf8,
} from './moleculeSelect.js';
import { parseXYZ } from '../util/xyzLoader.js';
import { validateParsedXYZ } from '../util/constraints.js';
import { showErrorBanner } from './errorBanner.js';
import { isLikelySmiles } from '../util/smilesLoader.js';
import { elInfo } from '../elements.js';
import { defaultMassForZ } from '../physics/sim-model.js';
import { getCellParameters, buildCellFromParameters } from '../util/pbc.js';
import { OMOL25_ELEMENTS } from '../data/periodicTable.js';

const PERIOD_ROWS = [
  ['H', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'He'],
  ['Li', 'Be', '', '', '', '', '', '', '', '', '', '', 'B', 'C', 'N', 'O', 'F', 'Ne'],
  ['Na', 'Mg', '', '', '', '', '', '', '', '', '', '', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
  [
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'Ge',
    'As',
    'Se',
    'Br',
    'Kr',
  ],
  [
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
  ],
  [
    'Cs',
    'Ba',
    'La',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
  ],
  ['Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'],
];

const LANTH = [
  'La',
  'Ce',
  'Pr',
  'Nd',
  'Pm',
  'Sm',
  'Eu',
  'Gd',
  'Tb',
  'Dy',
  'Ho',
  'Er',
  'Tm',
  'Yb',
  'Lu',
];

const ACTIN = [
  'Ac',
  'Th',
  'Pa',
  'U',
  'Np',
  'Pu',
  'Am',
  'Cm',
  'Bk',
  'Cf',
  'Es',
  'Fm',
  'Md',
  'No',
  'Lr',
];

function createSection(id, title, { defaultOpen = false } = {}) {
  const section = document.createElement('div');
  section.className = 'panel-section';
  section.id = id;

  const header = document.createElement('button');
  header.className = 'panel-header';
  header.type = 'button';
  header.setAttribute('aria-expanded', defaultOpen ? 'true' : 'false');
  header.textContent = title;

  // Right-aligned +/- indicator
  const indicator = document.createElement('span');
  indicator.className = 'collapse-indicator';
  const setIndicator = () => {
    const open = header.getAttribute('aria-expanded') === 'true';
    indicator.textContent = open ? '−' : '+';
  };
  setIndicator();
  header.appendChild(indicator);

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
    setIndicator();
  });

  // Keep indicator synced if aria-expanded is changed programmatically
  try {
    const mo = new MutationObserver((muts) => {
      for (const m of muts) {
        if (m.attributeName === 'aria-expanded') {
          setIndicator();
          break;
        }
      }
    });
    mo.observe(header, { attributes: true, attributeFilter: ['aria-expanded'] });
  } catch {}

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
      /* Hide scrollbars but keep scrollable */
      #controlPanel { scrollbar-width: none; -ms-overflow-style: none; }
      #controlPanel::-webkit-scrollbar { display: none; width: 0; height: 0; }

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
        position: relative;
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
      /* Classic switch variant */
      .toggle.switch { position: relative; width: 46px; height: 24px; padding: 0; color: transparent; }
      .toggle.switch::after {
        content: '';
        position: absolute; top: 2px; left: 2px;
        width: 20px; height: 20px; border-radius: 50%;
        background: #ffffff; box-shadow: 0 1px 2px rgba(0,0,0,0.4);
        transition: transform .18s ease;
      }
      .toggle.switch[data-on="true"] { background: #3B82F6; border-color: #3B82F6; }
      .toggle.switch[data-on="true"]::after { transform: translateX(22px); }

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
      #energyYLabel { position: absolute; left: 8px; top: 50%; transform: translateY(-50%) rotate(-90deg); transform-origin: center; font-size: 10px; color: var(--muted); white-space: nowrap; }
      #energyXLabel { position: absolute; left: 50%; bottom: 4px; transform: translateX(-50%); font-size: 10px; color: var(--muted); white-space: nowrap; }

      /* Small helpers */
      .stat { font-size: 12px; color: var(--muted); }
      .value { font-size: 12px; }
      .divider { border-top: 1px solid var(--border); margin: 8px 0; }
  /* Inline controls inside headers */
  .inline-controls { position: absolute; right: 28px; top: 50%; transform: translateY(-50%); display: inline-flex; gap: 8px; align-items: center; }
  .inline-controls label { font-size: 12px; color: var(--muted); }
  .inline-controls .state-label { font-size: 12px; color: var(--muted); min-width: 24px; text-align: left; }
  .radio { display: inline-flex; align-items: center; gap: 4px; }
  .radio input[type="radio"] { accent-color: #3B82F6; cursor: pointer; }
  /* Collapse indicator at far right of header */
  .panel-header .collapse-indicator { position: absolute; right: 10px; top: 50%; transform: translateY(-50%); color: var(--muted); font-weight: 700; pointer-events: none; }
      /* Fixed XR widget top-right */
      #xrModeWidget {
        position: fixed;
        right: 12px;
        top: 12px;
        z-index: 1000;
        background: rgba(30,38,48,0.85);
        color: #d8e6f3;
        border: 1px solid rgba(255,255,255,0.1);
        border-radius: 10px;
        padding: 6px 10px;
        display: inline-flex; align-items: center; gap: 8px;
      }
      /* Mobile top bar – single rounded rectangle; tabs always visible */
      #mobileTopBar {
        position: absolute;
        top: 8px;
        left: 8px;
        right: 8px;
        display: flex;
        flex-direction: column;
        gap: 0;
        z-index: 5;
        background: rgba(15,18,24,0.9);
        border: 1px solid rgba(255,255,255,0.08);
        border-radius: 10px;
        padding: 6px;
      }
      #mobileTopBar .tabs-row { display: flex; gap: 8px; }
      #mobileTopBar .tab {
        flex: 1 1 0;
        background: #2A2A2A;
        color: var(--text);
        border: 1px solid var(--border);
        border-radius: 8px;
        padding: 8px 10px;
        font-size: 12px;
        cursor: pointer;
      }
      #mobileTopBar .tab[data-active="true"] { background: var(--accent); border-color: var(--accent); }
      /* In-bar content area that expands the same box vertically */
      #mobileTopBar .mobile-sheet { display: none; padding: 8px; margin-top: 6px; background: var(--bg-panel); border: 1px solid rgba(255,255,255,0.08); border-radius: 8px; }
      #mobileTopBar[data-open="true"] .mobile-sheet { display: block; }
      /* When open previously blended with panel; now content lives inside the bar */
      /* #mobileTopBar[data-open="true"] { border-bottom-left-radius: 0; border-bottom-right-radius: 0; border-bottom: none; } */
      /* Mobile sheet behavior */
  #controlPanel[data-mobile-open="true"] { left: 8px; right: 8px; width: auto; top: 50px; border-top-left-radius: 0; border-top-right-radius: 0; border-top: none; }
      #controlPanel[data-mobile-open="true"] .panel-section { display: none; }
      #controlPanel[data-mobile-open="true"] .panel-section[data-mobile-active="true"] { display: block; }
    `;
    document.head.appendChild(style);
  }

  // Live Metrics (open)
  const live = createSection('section-live-stats', 'Live Metrics', { defaultOpen: true });
  {
    const statsRow = document.createElement('div');
    statsRow.className = 'row';
    // Energy first, then temperature, then RPS
    const instEnergy = document.createElement('span');
    instEnergy.id = 'instEnergy';
    instEnergy.className = 'value mono';
    instEnergy.textContent = 'E: —';
    const instTemp = document.createElement('span');
    instTemp.id = 'instTemp';
    instTemp.className = 'value mono';
    instTemp.textContent = 'T: — K';
    const rps = document.createElement('span');
    rps.id = 'rpsLabel';
    rps.className = 'value mono';
    rps.textContent = 'RPS: —';
    statsRow.append(instEnergy, instTemp, rps);
    live.content.appendChild(statsRow);

    const plot = document.createElement('div');
    plot.id = 'energyPlot';
    const canvas = document.createElement('canvas');
    canvas.id = 'energyCanvas';
    canvas.style.position = 'absolute';
    canvas.style.left = '0';
    canvas.style.top = '0';
    const yLabel = document.createElement('div');
    yLabel.id = 'energyYLabel';
    yLabel.textContent = 'Energy';
    const xLabel = document.createElement('div');
    xLabel.id = 'energyXLabel';
    xLabel.textContent = 'Time(step)';
    plot.append(canvas, yLabel, xLabel);
    live.content.appendChild(plot);
  }

  // Selection (collapsed by default) — shows selected atom/bond info and mini periodic table
  // Note: Periodic table cells are intentionally non-clickable for v1; this is a purely informative view.
  // VR/AR parity plan (future): Mirror this info in XR by rendering a small Babylon GUI (AdvancedDynamicTexture)
  // panel pinned to the controller ray or wrist anchor. Subscribe to the same state bus events (selectionChanged,
  // positionsChanged) and render simple GUI controls (text blocks + colored circles) for the spheres and labels.
  // For element highlighting, draw a compact grid of TextBlocks/images; apply highlight style (border/emissive) when selected.
  // Wire via viewerApi.vr.* helpers and ensure the panel is disabled in desktop to avoid duplication.
  const selSec = createSection('section-selection', 'Selection', { defaultOpen: false });
  {
    // Track auto-expand state within session
    let hasAutoExpandedSelection = false;
    let userCollapsedSelection = false;
    // Detect user toggles of this section to avoid auto re-open after manual collapse
    // Hook the section header click to toggle userCollapsedSelection when closing
    setTimeout(() => {
      try {
        const hdr = selSec.section.querySelector('.panel-header');
        const content = selSec.content;
        if (hdr && content) {
          hdr.addEventListener('click', () => {
            const willCollapse =
              content.getAttribute('data-collapsed') !== 'true' &&
              hdr.getAttribute('aria-expanded') === 'true';
            // If the user is collapsing after an auto-expand, mark it so we don't auto-open again
            if (willCollapse) userCollapsedSelection = true;
          });
        }
      } catch {}
    }, 0);
    // Styles specific to selection section (only once)
    const styleId = 'selectionSectionStyles';
    if (!document.getElementById(styleId)) {
      const s = document.createElement('style');
      s.id = styleId;
      s.textContent = `
        /* Selection section */
        #section-selection .panel-content { padding-top:6px; padding-bottom:8px; }
  /* Selection header: keep title left-aligned; spheres on the right */
  #section-selection .panel-header { padding-left: 10px; }
        #section-selection .panel-header .header-spheres { position:absolute; right:10px; top:50%; transform: translateY(-50%); display:flex; gap:6px; pointer-events:none; }
        /* Centered selection name in header */
        #section-selection .panel-header .sel-name {
          position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);
          max-width: calc(100% - 140px); /* leave room for left text and right spheres */
          text-align: center; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
          font-weight: 600; color: var(--text); pointer-events: none;
        }
        #section-selection .panel-header .header-spheres .atom-sphere { width:20px; height:20px; }
        #section-selection .atom-sphere { width:30px; height:30px; border-radius:50%; border:1px solid rgba(255,255,255,0.15); box-shadow: inset 0 0 6px rgba(0,0,0,0.5); }
  #section-selection .info { font-size:12px; color: var(--text); line-height:1.4; }
  #section-selection .info #selPosition { font-size: 11px; color: var(--text); }
        #section-selection .info .label { color: var(--muted); margin-right:6px; }
  /* Compact rotate buttons (smaller button, same icon size) */
  #section-selection #bondRotateWrap .btn { width: 28px; height: 28px; padding: 2px; display:inline-flex; align-items:center; justify-content:center; font-size: 18px; line-height: 1; }
        #section-selection .grid { margin-top:8px; display:grid; grid-template-columns: repeat(10, 1fr); gap:4px; }
        #section-selection .pt-el { text-align:center; padding:0; border-radius:2px; border:1px solid var(--border); font-size:10px; color: var(--text); background:#171717; user-select:none; box-sizing:border-box; }
        #section-selection .pt-el.highlight { border-color: var(--accent); box-shadow: 0 0 0 1px var(--accent) inset; }
        #section-selection .pt-el.omol25 { background: linear-gradient(135deg, rgba(134, 174, 252, 0.32), rgba(79, 139, 255, 0.12)), #1d222b; color:#f2f6ff; border-color: rgba(134, 174, 252, 0.7); box-shadow: inset 0 0 0 1px rgba(134, 174, 252, 0.45); }
      `;
      document.head.appendChild(s);
    }

    // Color palette for common elements (CSS colors; approximate Babylon colors)
    const EL_COLORS = {
      H: '#f2f7ff',
      C: '#383d47',
      N: '#6290f7',
      O: '#f25454',
      S: '#f7d64a',
      Cl: '#7bd47b',
      F: '#8ef06f',
      P: '#ffa94d',
      Br: '#a3552e',
      I: '#8a2be2',
      Si: '#f7786b',
      Fe: '#c44e4e',
      Co: '#cc6699',
      Ni: '#6fbf73',
      Cu: '#b87333',
      Zn: '#9bb5c1',
      Ag: '#c0c0c0',
      Au: '#ffd700',
      Na: '#abdee6',
      K: '#ffc8dd',
      Ca: '#ffe1a8',
      Li: '#cdeac0',
      He: '#c7d2fe',
      Ne: '#a7f3d0',
    };
    const EL_NAMES = {
      H: 'Hydrogen',
      C: 'Carbon',
      N: 'Nitrogen',
      O: 'Oxygen',
      S: 'Sulfur',
      Cl: 'Chlorine',
      F: 'Fluorine',
      P: 'Phosphorus',
      Br: 'Bromine',
      I: 'Iodine',
      Si: 'Silicon',
      Fe: 'Iron',
      Co: 'Cobalt',
      Ni: 'Nickel',
      Cu: 'Copper',
      Zn: 'Zinc',
      Ag: 'Silver',
      Au: 'Gold',
      Na: 'Sodium',
      K: 'Potassium',
      Ca: 'Calcium',
      Li: 'Lithium',
      He: 'Helium',
      Ne: 'Neon',
    };
    const SYMBOL_TO_Z = {
      H: 1,
      He: 2,
      Li: 3,
      Be: 4,
      B: 5,
      C: 6,
      N: 7,
      O: 8,
      F: 9,
      Ne: 10,
      Na: 11,
      Mg: 12,
      Al: 13,
      Si: 14,
      P: 15,
      S: 16,
      Cl: 17,
      Ar: 18,
      K: 19,
      Ca: 20,
    };
    const Z_TO_SYMBOL = Object.fromEntries(Object.entries(SYMBOL_TO_Z).map(([k, v]) => [v, k]));
    function getSymbol(el) {
      if (!el) return 'X';
      if (typeof el === 'string') return el;
      if (typeof el === 'number') return Z_TO_SYMBOL[el] || 'X';
      return el.symbol || el.sym || el.S || Z_TO_SYMBOL[el.Z || el.z || el.atomicNumber] || 'X';
    }
    function fmtPos(p) {
      if (!p) return '(-,-,-)';
      return `(${p.x.toFixed(1)}, ${p.y.toFixed(1)}, ${p.z.toFixed(1)})`;
    }
    function dist(a, b) {
      const dx = a.x - b.x,
        dy = a.y - b.y,
        dz = a.z - b.z;
      return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    // UI structure — spheres moved into header (smaller) and not shown in content
    const selHeader = selSec.section.querySelector('.panel-header');
    const headerSpheres = document.createElement('div');
    headerSpheres.className = 'header-spheres';
    const headerName = document.createElement('div');
    headerName.id = 'selElementName';
    headerName.className = 'sel-name';
    headerName.textContent = '—';
    const sphereA = document.createElement('div');
    sphereA.className = 'atom-sphere';
    sphereA.style.display = 'none';
    sphereA.id = 'selSphereA';
    const sphereB = document.createElement('div');
    sphereB.className = 'atom-sphere';
    sphereB.style.display = 'none';
    sphereB.id = 'selSphereB';
    headerSpheres.append(sphereA, sphereB);
    if (selHeader) {
      selHeader.appendChild(headerName);
      selHeader.appendChild(headerSpheres);
    }
    // Helper: shrink header name text to fit available width
    function fitHeaderName(node) {
      try {
        if (!node) return;
        // Reset to base size (inherit) before measuring
        node.style.fontSize = '';
        const maxLoops = 8;
        let loops = 0;
        let size = parseFloat(getComputedStyle(node).fontSize) || 13;
        while (node.scrollWidth > node.clientWidth && size > 10 && loops < maxLoops) {
          size -= 1;
          loops += 1;
          node.style.fontSize = size + 'px';
        }
      } catch {}
    }

    const info = document.createElement('div');
    info.className = 'info';
    info.innerHTML = `
      <div><span class="label">Position:</span><span id="selPosition">(-,-,-)</span></div>
      <div><span class="label">Atomic weight (amu):</span><span id="selAtomicWeight">—</span></div>
      <div><span class="label">vdW radius (Å):</span><span id="selVdw">—</span></div>
      <div id="bondRow" class="row">
        <span class="label">Bond:</span>
        <span id="bondLength" class="mono" style="display:inline-block; width:8ch; text-align:right; background:#171717; border:1px solid var(--border); border-radius:8px; padding:4px 6px;">N/A</span>
        <span id="bondRotateWrap" style="display:none; gap:6px; align-items:center;">
          <button id="bondRotMinus" class="btn" title="Rotate counter-clockwise">↺</button>
          <button id="bondRotPlus" class="btn" title="Rotate clockwise">↻</button>
        </span>
      </div>
    `;

    // Legacy test alias: some older tests expect an element with id "rotateBtns" to toggle visibility
    // alongside the new #bondRotateWrap. Create a lightweight sibling that mirrors display state.
    try {
      const bondRowEl = info.querySelector('#bondRow');
      const rotWrapEl = info.querySelector('#bondRotateWrap');
      if (bondRowEl && rotWrapEl && !info.querySelector('#rotateBtns')) {
        const legacy = document.createElement('span');
        legacy.id = 'rotateBtns';
        legacy.style.display = 'none';
        legacy.style.gap = '6px';
        legacy.style.alignItems = 'center';
        bondRowEl.insertBefore(legacy, rotWrapEl.nextSibling);
      }
    } catch {}

    // Full periodic table layout as a fixed HTML table (18 columns).
    // Lanthanides/Actinides shown as separate indented rows. Using a table ensures consistent cell sizing and
    // visible grid lines for rows/columns.
    const mini = document.createElement('table');
    mini.id = 'miniPeriodic';
    mini.setAttribute('role', 'grid');
    const miniBody = document.createElement('tbody');
    mini.appendChild(miniBody);
    // Table structure styles
    const ptStyleId = 'miniPeriodicStyles';
    if (!document.getElementById(ptStyleId)) {
      const st = document.createElement('style');
      st.id = ptStyleId;
      st.textContent = `
        #miniPeriodic { width:100%; table-layout: fixed; border-collapse: collapse; margin-top:6px; background:#232831; border:1px solid rgba(255,255,255,0.08); }
        #miniPeriodic tr { height:18px; }
        #miniPeriodic td.pt-el { text-align:center; padding:0; height:18px; line-height:18px; font-size:10px; color: var(--text); background:#1b2027; border:1px solid rgba(255,255,255,0.08); box-sizing:border-box; cursor:pointer; user-select:none; }
        #miniPeriodic td.pt-el.empty { background:#232831; border:none; }
        #miniPeriodic td.pt-el.highlight { border-color: var(--accent); box-shadow: inset 0 0 0 1px var(--accent); }
        #miniPeriodic td.pt-el.omol25 {
          background: linear-gradient(135deg, rgba(134, 174, 252, 0.28), rgba(79, 139, 255, 0.08)) , #1f2531;
          color: #f2f6ff;
          border-color: rgba(134, 174, 252, 0.7);
          box-shadow: inset 0 0 0 1px rgba(134, 174, 252, 0.45);
        }
      `;
      document.head.appendChild(st);
    }

    const OMOL25_SET = new Set(OMOL25_ELEMENTS);
    const dragAddState = {
      active: false,
      pointerId: null,
      symbol: null,
      startX: 0,
      startY: 0,
      currentX: 0,
      currentY: 0,
      dragging: false,
      sourceCell: null,
    };
    const previewId = 'periodicElementPreview';
    let previewEl = typeof document !== 'undefined' ? document.getElementById(previewId) : null;

    function ensurePreview() {
      if (previewEl) return previewEl;
      if (typeof document === 'undefined') return null;
      const el = document.createElement('div');
      el.id = previewId;
      el.style.display = 'none';
      el.style.position = 'fixed';
      el.style.pointerEvents = 'none';
      el.style.zIndex = '2000';
      el.style.transform = 'translate(-50%, -50%)';
      el.style.padding = '6px 10px';
      el.style.borderRadius = '999px';
      el.style.background = 'rgba(37,44,54,0.9)';
      el.style.color = '#f0f6ff';
      el.style.fontSize = '12px';
      el.style.border = '1px solid rgba(160,190,255,0.4)';
      el.style.boxShadow = '0 4px 10px rgba(0,0,0,0.3)';
      document.body.appendChild(el);
      previewEl = el;
      return previewEl;
    }

    function hidePreview() {
      if (previewEl) previewEl.style.display = 'none';
    }

    function showPreview(sym, x, y) {
      const el = ensurePreview();
      if (!el) return;
      el.textContent = sym;
      el.style.left = `${x}px`;
      el.style.top = `${y}px`;
      el.style.display = 'block';
    }

    function updatePreviewPosition(x, y) {
      if (!previewEl) return;
      previewEl.style.left = `${x}px`;
      previewEl.style.top = `${y}px`;
    }

    function resetDragAddState() {
      dragAddState.active = false;
      dragAddState.pointerId = null;
      dragAddState.symbol = null;
      dragAddState.startX = 0;
      dragAddState.startY = 0;
      dragAddState.currentX = 0;
      dragAddState.currentY = 0;
      dragAddState.dragging = false;
      dragAddState.sourceCell = null;
    }

    function projectToScene(clientX, clientY) {
      try {
        const api = window.viewerApi || window._viewer || null;
        if (!api || !api.scene || !api.camera) return null;
        const engine =
          api.scene.getEngine && typeof api.scene.getEngine === 'function'
            ? api.scene.getEngine()
            : null;
        const canvas = engine && typeof engine.getRenderingCanvas === 'function'
          ? engine.getRenderingCanvas()
          : null;
        if (!canvas) return null;
        const rect = canvas.getBoundingClientRect();
        const inside =
          clientX >= rect.left &&
          clientX <= rect.right &&
          clientY >= rect.top &&
          clientY <= rect.bottom;
        const px = clientX - rect.left;
        const py = clientY - rect.top;
        const camera = api.camera;
        const scene = api.scene;
        const ray = scene.createPickingRay(px, py, BABYLON.Matrix.Identity(), camera);
        const camPos = camera.position
          ? camera.position
          : camera.getFrontPosition
            ? camera.getFrontPosition(0)
            : { x: 0, y: 0, z: -10 };
        const target = camera.target || { x: 0, y: 0, z: 0 };
        const normal = new BABYLON.Vector3(
          target.x - camPos.x,
          target.y - camPos.y,
          target.z - camPos.z
        );
        if (normal.lengthSquared() < 1e-6) {
          normal.x = 0;
          normal.y = 1;
          normal.z = 0;
        }
        normal.normalize();
        const planePoint = new BABYLON.Vector3(target.x, target.y, target.z);
        const denom = BABYLON.Vector3.Dot(normal, ray.direction);
        let point = planePoint;
        if (Math.abs(denom) >= 1e-6) {
          const t = BABYLON.Vector3.Dot(planePoint.subtract(ray.origin), normal) / denom;
          if (Number.isFinite(t)) {
            point = ray.origin.add(ray.direction.scale(t));
          }
        }
        return {
          inside,
          point: { x: point.x, y: point.y, z: point.z },
        };
      } catch {
        return null;
      }
    }

    function addAtomAtOrigin(sym) {
      try {
        const api = window.viewerApi || window._viewer || null;
        if (!api) return;
        if (typeof api.addAtomAtOrigin === 'function') {
          Promise.resolve(api.addAtomAtOrigin(sym)).catch(() => {});
        } else if (typeof api.addAtom === 'function') {
          Promise.resolve(api.addAtom({ element: sym })).catch(() => {});
        }
      } catch {}
    }

    function addAtomAtPoint(sym, point) {
      try {
        const api = window.viewerApi || window._viewer || null;
        if (!api) return;
        if (typeof api.addAtomAtPosition === 'function') {
          Promise.resolve(api.addAtomAtPosition(sym, point)).catch(() => {});
        } else {
          addAtomAtOrigin(sym);
        }
      } catch {}
    }

    function handleAddMove(e) {
      if (!dragAddState.active || dragAddState.pointerId !== e.pointerId) return;
      dragAddState.currentX = e.clientX;
      dragAddState.currentY = e.clientY;
      const dx = dragAddState.currentX - dragAddState.startX;
      const dy = dragAddState.currentY - dragAddState.startY;
      if (!dragAddState.dragging) {
        const dist = Math.hypot(dx, dy);
        if (dist >= 6) {
          dragAddState.dragging = true;
          dragAddState.sourceCell?.classList?.add('drag-source');
          showPreview(dragAddState.symbol, dragAddState.currentX, dragAddState.currentY);
        }
      } else {
        updatePreviewPosition(dragAddState.currentX, dragAddState.currentY);
      }
      if (typeof e?.preventDefault === 'function') e.preventDefault();
    }

    function endAddSession(e) {
      if (!dragAddState.active || dragAddState.pointerId !== e.pointerId) return;
      window.removeEventListener('pointermove', handleAddMove, true);
      window.removeEventListener('pointerup', endAddSession, true);
      window.removeEventListener('pointercancel', cancelAddSession, true);
      dragAddState.sourceCell?.classList?.remove('drag-source');
      hidePreview();
      const sym = dragAddState.symbol;
      const dragging = dragAddState.dragging;
      const currentX = dragAddState.currentX;
      const currentY = dragAddState.currentY;
      resetDragAddState();
      if (!sym) return;
      if (!dragging) {
        addAtomAtOrigin(sym);
        return;
      }
      const projection = projectToScene(currentX, currentY);
      if (projection && projection.inside && projection.point) {
        addAtomAtPoint(sym, projection.point);
      } else {
        addAtomAtOrigin(sym);
      }
    }

    function cancelAddSession(e) {
      if (dragAddState.pointerId !== e.pointerId) return;
      window.removeEventListener('pointermove', handleAddMove, true);
      window.removeEventListener('pointerup', endAddSession, true);
      window.removeEventListener('pointercancel', cancelAddSession, true);
      dragAddState.sourceCell?.classList?.remove('drag-source');
      hidePreview();
      resetDragAddState();
    }

    function beginAddSession(e, sym, cell) {
      if (!OMOL25_SET.has(sym)) return;
      ensurePreview();
      dragAddState.active = true;
      dragAddState.pointerId = e.pointerId;
      dragAddState.symbol = sym;
      dragAddState.startX = e.clientX;
      dragAddState.startY = e.clientY;
      dragAddState.currentX = e.clientX;
      dragAddState.currentY = e.clientY;
      dragAddState.dragging = false;
      dragAddState.sourceCell = cell;
      window.addEventListener('pointermove', handleAddMove, true);
      window.addEventListener('pointerup', endAddSession, true);
      window.addEventListener('pointercancel', cancelAddSession, true);
      if (typeof e?.preventDefault === 'function') e.preventDefault();
    }
    // Periods 1..7 layout defined at module scope (PERIOD_ROWS, LANTH, ACTIN)
    // Symbol->full name map (118). For brevity, keep names concise.
    const SYMBOL_TO_NAME = {
      H: 'Hydrogen',
      He: 'Helium',
      Li: 'Lithium',
      Be: 'Beryllium',
      B: 'Boron',
      C: 'Carbon',
      N: 'Nitrogen',
      O: 'Oxygen',
      F: 'Fluorine',
      Ne: 'Neon',
      Na: 'Sodium',
      Mg: 'Magnesium',
      Al: 'Aluminum',
      Si: 'Silicon',
      P: 'Phosphorus',
      S: 'Sulfur',
      Cl: 'Chlorine',
      Ar: 'Argon',
      K: 'Potassium',
      Ca: 'Calcium',
      Sc: 'Scandium',
      Ti: 'Titanium',
      V: 'Vanadium',
      Cr: 'Chromium',
      Mn: 'Manganese',
      Fe: 'Iron',
      Co: 'Cobalt',
      Ni: 'Nickel',
      Cu: 'Copper',
      Zn: 'Zinc',
      Ga: 'Gallium',
      Ge: 'Germanium',
      As: 'Arsenic',
      Se: 'Selenium',
      Br: 'Bromine',
      Kr: 'Krypton',
      Rb: 'Rubidium',
      Sr: 'Strontium',
      Y: 'Yttrium',
      Zr: 'Zirconium',
      Nb: 'Niobium',
      Mo: 'Molybdenum',
      Tc: 'Technetium',
      Ru: 'Ruthenium',
      Rh: 'Rhodium',
      Pd: 'Palladium',
      Ag: 'Silver',
      Cd: 'Cadmium',
      In: 'Indium',
      Sn: 'Tin',
      Sb: 'Antimony',
      Te: 'Tellurium',
      I: 'Iodine',
      Xe: 'Xenon',
      Cs: 'Cesium',
      Ba: 'Barium',
      La: 'Lanthanum',
      Ce: 'Cerium',
      Pr: 'Praseodymium',
      Nd: 'Neodymium',
      Pm: 'Promethium',
      Sm: 'Samarium',
      Eu: 'Europium',
      Gd: 'Gadolinium',
      Tb: 'Terbium',
      Dy: 'Dysprosium',
      Ho: 'Holmium',
      Er: 'Erbium',
      Tm: 'Thulium',
      Yb: 'Ytterbium',
      Lu: 'Lutetium',
      Hf: 'Hafnium',
      Ta: 'Tantalum',
      W: 'Tungsten',
      Re: 'Rhenium',
      Os: 'Osmium',
      Ir: 'Iridium',
      Pt: 'Platinum',
      Au: 'Gold',
      Hg: 'Mercury',
      Tl: 'Thallium',
      Pb: 'Lead',
      Bi: 'Bismuth',
      Po: 'Polonium',
      At: 'Astatine',
      Rn: 'Radon',
      Fr: 'Francium',
      Ra: 'Radium',
      Ac: 'Actinium',
      Th: 'Thorium',
      Pa: 'Protactinium',
      U: 'Uranium',
      Np: 'Neptunium',
      Pu: 'Plutonium',
      Am: 'Americium',
      Cm: 'Curium',
      Bk: 'Berkelium',
      Cf: 'Californium',
      Es: 'Einsteinium',
      Fm: 'Fermium',
      Md: 'Mendelevium',
      No: 'Nobelium',
      Lr: 'Lawrencium',
      Rf: 'Rutherfordium',
      Db: 'Dubnium',
      Sg: 'Seaborgium',
      Bh: 'Bohrium',
      Hs: 'Hassium',
      Mt: 'Meitnerium',
      Ds: 'Darmstadtium',
      Rg: 'Roentgenium',
      Cn: 'Copernicium',
      Nh: 'Nihonium',
      Fl: 'Flerovium',
      Mc: 'Moscovium',
      Lv: 'Livermorium',
      Ts: 'Tennessine',
      Og: 'Oganesson',
    };

    // Approximate standard atomic weights (amu) keyed by symbol (for UI display).
    // Values are representative; for synthetic/unstable elements, common isotopic masses are used.
    const MASS_BY_SYMBOL = {
      H: 1.008,
      He: 4.0026,
      Li: 6.94,
      Be: 9.0122,
      B: 10.81,
      C: 12.011,
      N: 14.007,
      O: 15.999,
      F: 18.998,
      Ne: 20.18,
      Na: 22.99,
      Mg: 24.305,
      Al: 26.982,
      Si: 28.085,
      P: 30.974,
      S: 32.06,
      Cl: 35.45,
      Ar: 39.948,
      K: 39.0983,
      Ca: 40.078,
      Sc: 44.9559,
      Ti: 47.867,
      V: 50.9415,
      Cr: 51.9961,
      Mn: 54.938,
      Fe: 55.845,
      Co: 58.933,
      Ni: 58.693,
      Cu: 63.546,
      Zn: 65.38,
      Ga: 69.723,
      Ge: 72.63,
      As: 74.922,
      Se: 78.971,
      Br: 79.904,
      Kr: 83.798,
      Rb: 85.468,
      Sr: 87.62,
      Y: 88.906,
      Zr: 91.224,
      Nb: 92.906,
      Mo: 95.95,
      Tc: 98,
      Ru: 101.07,
      Rh: 102.905,
      Pd: 106.42,
      Ag: 107.868,
      Cd: 112.414,
      In: 114.818,
      Sn: 118.71,
      Sb: 121.76,
      Te: 127.6,
      I: 126.904,
      Xe: 131.293,
      Cs: 132.905,
      Ba: 137.327,
      La: 138.905,
      Ce: 140.116,
      Pr: 140.908,
      Nd: 144.242,
      Pm: 145,
      Sm: 150.36,
      Eu: 151.964,
      Gd: 157.25,
      Tb: 158.925,
      Dy: 162.5,
      Ho: 164.93,
      Er: 167.259,
      Tm: 168.934,
      Yb: 173.045,
      Lu: 174.967,
      Hf: 178.49,
      Ta: 180.948,
      W: 183.84,
      Re: 186.207,
      Os: 190.23,
      Ir: 192.217,
      Pt: 195.084,
      Au: 196.96657,
      Hg: 200.592,
      Tl: 204.38,
      Pb: 207.2,
      Bi: 208.98,
      Po: 209,
      At: 210,
      Rn: 222,
      Fr: 223,
      Ra: 226,
      Ac: 227,
      Th: 232.0377,
      Pa: 231.0359,
      U: 238.0289,
      Np: 237,
      Pu: 244,
      Am: 243,
      Cm: 247,
      Bk: 247,
      Cf: 251,
      Es: 252,
      Fm: 257,
      Md: 258,
      No: 259,
      Lr: 266,
      Rf: 267,
      Db: 268,
      Sg: 269,
      Bh: 270,
      Hs: 277,
      Mt: 278,
      Ds: 281,
      Rg: 282,
      Cn: 285,
      Nh: 286,
      Fl: 289,
      Mc: 290,
      Lv: 293,
      Ts: 294,
      Og: 294,
    };

    // Build DOM rows
    function makeRow(symbols) {
      const tr = document.createElement('tr');
      for (let c = 0; c < 18; c++) {
        const sym = symbols[c] || '';
        const td = document.createElement('td');
        td.className = 'pt-el' + (sym ? '' : ' empty');
        if (sym) {
          td.setAttribute('data-symbol', sym);
          td.textContent = sym;
          if (OMOL25_SET.has(sym)) {
            td.classList.add('omol25');
            td.setAttribute('title', `${sym} — available in OMol25`);
          }
          td.addEventListener('click', () => {
            try {
              const api = getViewer();
              api && api.selection && api.selection.clear && api.selection.clear();
            } catch {}
            __overrideElementSym = sym;
            updateFromSelection();
          });
          td.addEventListener('pointerdown', (e) => {
            const isPrimary =
              e.isPrimary !== false &&
              (e.button === 0 || e.button === undefined || e.button === -1);
            if (!isPrimary) return;
            beginAddSession(e, sym, td);
          });
        } else {
          td.textContent = '';
        }
        tr.appendChild(td);
      }
      return tr;
    }
    PERIOD_ROWS.forEach((r) => miniBody.appendChild(makeRow(r)));
    const lanthRow = new Array(18).fill('');
    LANTH.forEach((sym, idx) => {
      const col = idx + 2;
      if (col < lanthRow.length) lanthRow[col] = sym;
    });
    const actinRow = new Array(18).fill('');
    ACTIN.forEach((sym, idx) => {
      const col = idx + 2;
      if (col < actinRow.length) actinRow[col] = sym;
    });
    miniBody.appendChild(makeRow(lanthRow));
    miniBody.appendChild(makeRow(actinRow));

    selSec.content.append(info, mini);

    // Live updater bound to selection changes
    function getViewer() {
      try {
        return window.viewerApi || window._viewer;
      } catch {
        return null;
      }
    }
    function highlight(symbols) {
      const cells = mini.querySelectorAll('.pt-el');
      cells.forEach((c) => c.classList.remove('highlight'));
      const set = new Set((symbols || []).filter(Boolean));
      set.forEach((sym) => {
        const el = mini.querySelector(`.pt-el[data-symbol="${sym}"]`);
        if (el) el.classList.add('highlight');
      });
    }
    function setSphere(el, sym) {
      if (!sym) {
        el.style.display = 'none';
        return;
      }
      el.style.display = 'block';
      const color = EL_COLORS[sym] || '#9aa0a6';
      el.style.background = `radial-gradient( circle at 30% 30%, #ffffff66, transparent 40%), ${color}`;
    }
    // If a periodic table element was clicked, we set an override symbol and clear it when a real selection happens.
    let __overrideElementSym = null;
    function symbolToName(sym) {
      return EL_NAMES[sym] || SYMBOL_TO_NAME[sym] || sym;
    }
    function ensureSelectionOpenOnce() {
      try {
        // Only auto-open once per session and only if the user hasn't manually collapsed it
        if (!hasAutoExpandedSelection && !userCollapsedSelection) {
          // Expand section
          selSec.content.removeAttribute('data-collapsed');
          const hdr = selSec.section.querySelector('.panel-header');
          if (hdr) hdr.setAttribute('aria-expanded', 'true');
          hasAutoExpandedSelection = true;
        }
      } catch {}
    }

    function updateFromSelection() {
      const api = getViewer();
      if (!api) return;
      const st = api.state;
      const sel = (st && st.selection) || { kind: null };
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_UI)
          console.log('[panel] updateFromSelection sel=', sel);
      } catch {}
      const elNameNode = selSec.section.querySelector('#selElementName');
      const posNode = selSec.content.querySelector('#selPosition');
      const weightNode = selSec.content.querySelector('#selAtomicWeight');
      const vdwNode = selSec.content.querySelector('#selVdw');
      const bondRow = selSec.content.querySelector('#bondRow');
      const bondLenNode = selSec.content.querySelector('#bondLength');
      const rotWrap = selSec.content.querySelector('#bondRotateWrap');
      const rotWrapLegacy = selSec.content.querySelector('#rotateBtns');
      const rotMinus = selSec.content.querySelector('#bondRotMinus');
      const rotPlus = selSec.content.querySelector('#bondRotPlus');
      // Bind handlers once
      if (!rotMinus.__bound) {
        const EPS = 0.1;
        rotMinus.addEventListener('click', () => {
          try {
            const v = getViewer();
            v && v.manipulation && v.manipulation.rotateBond && v.manipulation.rotateBond(-EPS);
          } catch {}
        });
        rotMinus.__bound = true;
        rotPlus.addEventListener('click', () => {
          try {
            const v = getViewer();
            v && v.manipulation && v.manipulation.rotateBond && v.manipulation.rotateBond(EPS);
          } catch {}
        });
        rotPlus.__bound = true;
      }
      function weightForSym(sym) {
        if (MASS_BY_SYMBOL[sym] != null) return MASS_BY_SYMBOL[sym];
        const Z = SYMBOL_TO_Z[sym];
        if (!Z) return null;
        const m = defaultMassForZ(Z);
        return typeof m === 'number' && isFinite(m) ? m : null;
      }
      function vdwForSym(sym) {
        const info = elInfo(sym);
        return info && typeof info.vdw === 'number' ? info.vdw : null;
      }
      if (__overrideElementSym) {
        ensureSelectionOpenOnce();
        // Show virtual element selection
        setSphere(sphereA, __overrideElementSym);
        setSphere(sphereB, null);
        elNameNode.textContent = symbolToName(__overrideElementSym);
        fitHeaderName(elNameNode);
        posNode.textContent = '(-,-,-)';
        const mw = weightForSym(__overrideElementSym);
        const rv = vdwForSym(__overrideElementSym);
        weightNode.textContent = mw != null ? mw.toFixed(3) : '—';
        vdwNode.textContent = rv != null ? rv.toFixed(2) : '—';
        bondRow.style.display = 'block';
        bondLenNode.textContent = 'N/A';
        rotWrap.style.display = 'none';
        if (rotWrapLegacy) rotWrapLegacy.style.display = 'none';
        highlight([__overrideElementSym]);
      } else if (sel.kind === 'atom') {
        ensureSelectionOpenOnce();
        const idx = sel.data.index;
        const sym = getSymbol(st.elements[idx]);
        const name = symbolToName(sym);
        const pos = st.positions[idx];
        setSphere(sphereA, sym);
        setSphere(sphereB, null);
        elNameNode.textContent = name;
        fitHeaderName(elNameNode);
        posNode.textContent = fmtPos(pos);
        const mw = weightForSym(sym);
        const rv = vdwForSym(sym);
        weightNode.textContent = mw != null ? mw.toFixed(3) : '—';
        vdwNode.textContent = rv != null ? rv.toFixed(2) : '—';
        bondRow.style.display = 'block';
        bondLenNode.textContent = 'N/A';
        rotWrap.style.display = 'none';
        if (rotWrapLegacy) rotWrapLegacy.style.display = 'none';
        highlight([sym]);
      } else if (sel.kind === 'bond') {
        ensureSelectionOpenOnce();
        const i = sel.data.i,
          j = sel.data.j;
        const symA = getSymbol(st.elements[i]);
        const symB = getSymbol(st.elements[j]);
        const nameA = symbolToName(symA);
        const nameB = symbolToName(symB);
        const posA = st.positions[i];
        const posB = st.positions[j];
        setSphere(sphereA, symA);
        setSphere(sphereB, symB);
        elNameNode.textContent = `${nameA} – ${nameB}`;
        fitHeaderName(elNameNode);
        posNode.textContent = `${fmtPos(posA)} – ${fmtPos(posB)}`;
        const mwA = weightForSym(symA),
          mwB = weightForSym(symB);
        const rvA = vdwForSym(symA),
          rvB = vdwForSym(symB);
        weightNode.textContent = `${mwA != null ? mwA.toFixed(3) : '—'} – ${mwB != null ? mwB.toFixed(3) : '—'}`;
        vdwNode.textContent = `${rvA != null ? rvA.toFixed(2) : '—'} – ${rvB != null ? rvB.toFixed(2) : '—'}`;
        const L = dist(posA, posB);
        bondRow.style.display = 'block';
        bondLenNode.textContent = `${L.toFixed(2)} Å`;
        rotWrap.style.display = 'inline-flex';
        if (rotWrapLegacy) rotWrapLegacy.style.display = 'inline-flex';
        try {
          if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_UI)
            console.log('[panel] showing rotate buttons for bond', i, j, 'L=', L.toFixed(3));
        } catch {}
        highlight([symA, symB]);
      } else {
        setSphere(sphereA, null);
        setSphere(sphereB, null);
        elNameNode.textContent = '—';
        fitHeaderName(elNameNode);
        posNode.textContent = '(-,-,-)';
        weightNode.textContent = '—';
        vdwNode.textContent = '—';
        bondRow.style.display = 'block';
        bondLenNode.textContent = 'N/A';
        rotWrap.style.display = 'none';
        highlight([]);
      }
    }
    // Initial draw and subscriptions
    updateFromSelection();
    // Keep header name fitting on resize
    try {
      window.addEventListener(
        'resize',
        () => {
          const n = selSec.section.querySelector('#selElementName');
          if (n) {
            // Give layout a moment to settle
            setTimeout(() => fitHeaderName(n), 0);
          }
        },
        { passive: true }
      );
    } catch {}
    try {
      const api = getViewer();
      if (api && api.state && api.state.bus) {
        api.state.bus.on('selectionChanged', () => {
          __overrideElementSym = null;
          updateFromSelection();
        });
        // Also keep positions updated in case of drag while selected
        api.state.bus.on('positionsChanged', updateFromSelection);
      }
    } catch {}
  }

  // Simulation (open) — Top utility toggles + Relax/MD toggles (side-by-side)
  const sim = createSection('section-simulation', 'Simulation', { defaultOpen: false });
  {
    // Top row: Show Forces
    const utilRow = document.createElement('div');
    utilRow.className = 'row';
    // Forces toggle
    const forcesT = makeToggle({
      id: 'toggleForces',
      labelOn: 'Show Forces: On',
      labelOff: 'Show Forces: Off',
      title: 'Render per-atom force vectors',
    });
    forcesT.addEventListener('click', () => {
      try {
        const api = window.viewerApi || window._viewer;
        if (!api) return;
        const on = forcesT.getAttribute('data-on') === 'true';
        api.setForceVectorsEnabled(on);
        const legacy = document.getElementById('btnToggleForces');
        if (legacy) legacy.textContent = on ? 'on' : 'off';
      } catch {}
    });
    utilRow.append(forcesT);

    // Relaxation toggle (mutually exclusive with MD)
    const togglesRow = document.createElement('div');
    togglesRow.className = 'row';
    const relaxToggle = makeToggle({
      id: 'toggleRelax',
      labelOn: 'Relaxation: On',
      labelOff: 'Relaxation: Off',
      title: 'Toggle geometry optimization',
    });
    // MD toggle
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
    function getApi() {
      try {
        return window.viewerApi || window._viewer;
      } catch {
        return null;
      }
    }
    function setBtnState(btn, on) {
      btn.setAttribute('data-on', String(!!on));
      btn.setAttribute('aria-checked', String(!!on));
      btn.textContent = on
        ? btn.id === 'toggleMD'
          ? 'Molecular Dynamics: On'
          : 'Relaxation: On'
        : btn.id === 'toggleMD'
          ? 'Molecular Dynamics: Off'
          : 'Relaxation: Off';
    }
    // Relaxation toggle wiring
    relaxToggle.addEventListener('click', (e) => {
      const api = getApi();
      if (!api) return;
      const on = relaxToggle.getAttribute('data-on') === 'true';
      // If turning on, stop MD and start relax run (single step semantics via run loop)
      if (on) {
        try {
          api.stopSimulation();
        } catch {}
        try {
          const shortRuns = typeof window !== 'undefined' && !!window.__MLIPVIEW_PANEL_SHORT_RUNS;
          const args = shortRuns ? { maxSteps: 2 } : {};
          api.startRelaxContinuous(args);
        } catch {}
      } else {
        try {
          api.stopSimulation();
        } catch {}
      }
      // Mirror legacy button text convention for compatibility
      try {
        const legacy = document.getElementById('btnRelaxRun');
        if (legacy) legacy.textContent = on ? 'stop' : 'run';
      } catch {}
    });
    // MD toggle wiring
    mdToggle.addEventListener('click', (e) => {
      const api = getApi();
      if (!api) return;
      const on = mdToggle.getAttribute('data-on') === 'true';
      if (on) {
        try {
          api.stopSimulation();
        } catch {}
        try {
          const shortRuns = typeof window !== 'undefined' && !!window.__MLIPVIEW_PANEL_SHORT_RUNS;
          const args = shortRuns ? { steps: 2 } : {};
          api.startMDContinuous(args);
        } catch {}
      } else {
        try {
          api.stopSimulation();
        } catch {}
      }
      try {
        const legacy = document.getElementById('btnMDRun');
        if (legacy) legacy.textContent = on ? 'stop' : 'run';
      } catch {}
    });

    togglesRow.append(relaxToggle, mdToggle);

    // Reflect simulation state from API into toggles (handles auto-start MD after load)
    function syncFromApiOnce() {
      try {
        const api = getApi();
        if (!api || !api.getMetrics) return;
        // Always query the current DOM toggles so this works across re-renders/tests
        const mdToggleEl = document.getElementById('toggleMD');
        const relaxToggleEl = document.getElementById('toggleRelax');
        if (!mdToggleEl || !relaxToggleEl) return;
        const m = api.getMetrics() || {};
        const kind = m.running || null;
        if (kind === 'md') {
          // MD running: set MD=On, Relax=Off
          if (mdToggleEl.getAttribute('data-on') !== 'true') setBtnState(mdToggleEl, true);
          if (relaxToggleEl.getAttribute('data-on') !== 'false') setBtnState(relaxToggleEl, false);
          try {
            const legacy = document.getElementById('btnMDRun');
            if (legacy) legacy.textContent = 'stop';
          } catch {}
          // Status element removed
        } else if (kind === 'relax') {
          // Relax running: set Relax=On, MD=Off
          if (relaxToggleEl.getAttribute('data-on') !== 'true') setBtnState(relaxToggleEl, true);
          if (mdToggleEl.getAttribute('data-on') !== 'false') setBtnState(mdToggleEl, false);
          try {
            const legacy = document.getElementById('btnRelaxRun');
            if (legacy) legacy.textContent = 'stop';
          } catch {}
          // Status element removed
        } else {
          // Neither running
          if (mdToggleEl.getAttribute('data-on') !== 'false') setBtnState(mdToggleEl, false);
          if (relaxToggleEl.getAttribute('data-on') !== 'false') setBtnState(relaxToggleEl, false);
          try {
            const legacyMD = document.getElementById('btnMDRun');
            if (legacyMD) legacyMD.textContent = 'run';
          } catch {}
          try {
            const legacyRX = document.getElementById('btnRelaxRun');
            if (legacyRX) legacyRX.textContent = 'run';
          } catch {}
          // Status element removed
        }
      } catch {}
    }
    // Initial sync now, then keep a lightweight periodic sync so programmatic runs update UI
    try {
      syncFromApiOnce();
    } catch {}
    try {
      if (!window.__MLIPVIEW_PANEL_SYNC_INTERVAL) {
        window.__MLIPVIEW_PANEL_SYNC_INTERVAL = setInterval(syncFromApiOnce, 250);
      }
    } catch {}

    sim.content.append(utilRow, togglesRow);

    // Sliders: temperature + friction (legacy tests expect both sliders to exist)
    initTemperatureSlider({ hudEl: sim.content, getViewer: () => null });
    try {
      initFrictionSlider({ hudEl: sim.content });
    } catch {}
  }

  // System (collapsed) — keep presets and SMILES, keep backend select
  const sys = createSection('section-system', 'System', { defaultOpen: false });
  {
    // Preset dropdown (existing helper)
    const presetsRow = document.createElement('div');
    presetsRow.className = 'row';
    const label = document.createElement('span');
    label.className = 'form-label';
    label.textContent = 'Preset:';
    presetsRow.appendChild(label);
    installMoleculeSelector({ hudEl: presetsRow, windowRef: window, documentRef: document });
    sys.content.appendChild(presetsRow);

    // SMILES input (single-line) + Generate button; Enter submits
    const smilesBlock = document.createElement('div');
    smilesBlock.className = 'block';
    const smilesLabel = document.createElement('div');
    smilesLabel.className = 'form-label';
    smilesLabel.textContent = 'SMILES';
    const inputRow = document.createElement('div');
    inputRow.className = 'input-row';
    const smilesInput = document.createElement('input');
    smilesInput.type = 'text';
    smilesInput.id = 'smilesInput';
    smilesInput.placeholder = 'Paste or type SMILES here';
    const smilesBtn = document.createElement('button');
    smilesBtn.id = 'smilesGoBtn';
    smilesBtn.className = 'btn';
    smilesBtn.textContent = 'Generate';
    smilesBtn.disabled = true;

    smilesInput.addEventListener('input', () => {
      smilesBtn.disabled = !smilesInput.value.trim();
    });
    smilesInput.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' && !smilesBtn.disabled) {
        smilesBtn.click();
      }
    });

    // Wire SMILES Generate: validate basic string and reload page with ?smiles=
    smilesBtn.addEventListener('click', () => {
      const s = (smilesInput.value || '').trim();
      if (!s) return;
      if (!isLikelySmiles(s)) {
        showErrorBanner('Invalid SMILES format');
        return;
      }
      try {
        const href = buildReloadUrlWithParam(window.location.href, 'smiles', s);
        window.location.assign ? window.location.assign(href) : (window.location.href = href);
      } catch (e) {
        showErrorBanner('SMILES navigation failed');
      }
    });

    inputRow.append(smilesInput, smilesBtn);
    smilesBlock.append(smilesLabel, inputRow);
    sys.content.appendChild(smilesBlock);

    // Upload XYZ button + hidden input
    const uploadRow = document.createElement('div');
    uploadRow.className = 'row';
    const uploadBtn = document.createElement('button');
    uploadBtn.id = 'uploadXyzBtn';
    uploadBtn.className = 'btn';
    uploadBtn.textContent = 'Upload XYZ';
    const fileInput = document.createElement('input');
    fileInput.type = 'file';
    fileInput.accept = '.xyz,text/plain';
    fileInput.style.display = 'none';
    fileInput.id = 'xyzFileInput';
    uploadBtn.addEventListener('click', () => fileInput.click());
    fileInput.addEventListener('change', async (ev) => {
      try {
        const f = ev.target.files && ev.target.files[0];
        if (!f) return;
        const text = await f.text();
        let parsed;
        try {
          parsed = parseXYZ(text);
        } catch (e) {
          showErrorBanner('XYZ parse failed');
          console.warn?.('[uploadXYZ] parse failed', e);
          return;
        }
        const v = validateParsedXYZ(parsed);
        if (!v.ok) {
          showErrorBanner(`XYZ invalid: ${v.error}`);
          console.warn?.('[uploadXYZ] validation failed', v.error);
          return;
        }
        // Encode to base64 (UTF-8 safe) and reload with ?molxyz
        const b64 = base64EncodeUtf8(text);
        const href = buildReloadUrlWithParam(window.location.href, 'molxyz', b64);
        window.location.assign ? window.location.assign(href) : (window.location.href = href);
      } catch (e) {
        showErrorBanner('XYZ upload failed');
        console.warn?.('[uploadXYZ] unexpected error', e);
      } finally {
        try {
          ev.target.value = '';
        } catch {}
      }
    });
    uploadRow.append(uploadBtn, fileInput);
    sys.content.appendChild(uploadRow);

    const sessionRow = document.createElement('div');
    sessionRow.className = 'row';
    const loadSessionBtn = document.createElement('button');
    loadSessionBtn.id = 'loadSessionBtn';
    loadSessionBtn.dataset.testid = 'session-load';
    loadSessionBtn.className = 'btn';
    loadSessionBtn.textContent = 'Load Session';
    const saveSessionBtn = document.createElement('button');
    saveSessionBtn.id = 'saveSessionBtn';
    saveSessionBtn.dataset.testid = 'session-save';
    saveSessionBtn.className = 'btn';
    saveSessionBtn.textContent = 'Save Session';
    const sessionInput = document.createElement('input');
    sessionInput.type = 'file';
    sessionInput.accept = '.json,application/json';
    sessionInput.style.display = 'none';
    loadSessionBtn.addEventListener('click', () => sessionInput.click());
    sessionInput.addEventListener('change', async (ev) => {
      try {
        const f = ev.target.files && ev.target.files[0];
        if (!f) return;
        const api = (typeof window !== 'undefined') ? window.viewerApi || window._viewer || null : null;
        if (!api?.session?.loadFromFile) {
          showErrorBanner('Session loading not available');
          return;
        }
        await api.session.loadFromFile(f);
      } catch (e) {
        showErrorBanner('Session load failed');
        console.warn?.('[sessionLoad] failed', e);
      } finally {
        try { ev.target.value = ''; } catch {}
      }
    });
    saveSessionBtn.addEventListener('click', () => {
      try {
        const api = (typeof window !== 'undefined') ? window.viewerApi || window._viewer || null : null;
        if (!api?.session?.saveToFile) {
          showErrorBanner('Session saving not available');
          return;
        }
        api.session.saveToFile();
      } catch (e) {
        showErrorBanner('Session save failed');
        console.warn?.('[sessionSave] failed', e);
      }
    });
    sessionRow.append(loadSessionBtn, saveSessionBtn, sessionInput);
    sys.content.appendChild(sessionRow);

    const libraryRow = document.createElement('div');
    libraryRow.className = 'row';
    const libraryLabel = document.createElement('span');
    libraryLabel.className = 'form-label';
    libraryLabel.textContent = 'Library:';
    const librarySelect = document.createElement('select');
    librarySelect.id = 'sessionLibrarySel';
    librarySelect.dataset.testid = 'session-library';
    librarySelect.disabled = true;
    librarySelect.innerHTML = '<option value="">Loading…</option>';
    libraryRow.append(libraryLabel, librarySelect);
    sys.content.appendChild(libraryRow);

    const loadLibraryEntry = async (id) => {
      if (!id) return;
      const api = (typeof window !== 'undefined') ? window.viewerApi || window._viewer || null : null;
      if (!api?.session?.loadFromLibrary) {
        showErrorBanner('Session library not available');
        return;
      }
      librarySelect.disabled = true;
      try {
        await api.session.loadFromLibrary(id);
      } catch (e) {
        showErrorBanner('Library load failed');
        console.warn?.('[sessionLibrary] load failed', e);
      } finally {
        librarySelect.disabled = false;
        librarySelect.value = '';
      }
    };

    librarySelect.addEventListener('change', (ev) => {
      const id = ev.target.value;
      if (!id) return;
      loadLibraryEntry(id);
    });

    (async () => {
      try {
        const api = (typeof window !== 'undefined') ? window.viewerApi || window._viewer || null : null;
        if (!api?.session?.getLibraryManifest) {
          librarySelect.innerHTML = '<option value="">Unavailable</option>';
          return;
        }
        const entries = await api.session.getLibraryManifest();
        if (!Array.isArray(entries) || entries.length === 0) {
          librarySelect.innerHTML = '<option value="">No entries</option>';
          librarySelect.disabled = true;
          return;
        }
        librarySelect.innerHTML = '<option value="">Select session…</option>';
        for (const entry of entries) {
          const opt = document.createElement('option');
          opt.value = entry.id;
          opt.textContent = entry.label || entry.id;
          if (entry.description) opt.title = entry.description;
          librarySelect.appendChild(opt);
        }
        librarySelect.disabled = false;
      } catch (e) {
        librarySelect.innerHTML = '<option value="">Load failed</option>';
        librarySelect.disabled = true;
        console.warn?.('[sessionLibrary] manifest load failed', e);
      }
    })();

    // Backend select only (PBC moved)
    const backendRow = document.createElement('div');
    backendRow.className = 'row';
    const backendSel = document.createElement('select');
    backendSel.id = 'forceProviderSel';
    backendSel.innerHTML = '<option value="fairchem">UMA (Remote)</option>';
    backendRow.append(backendSel);
    sys.content.appendChild(backendRow);
  }
  // Rendering section removed per request

  // Periodic (collapsed) — PBC + full unit cell parameters, with PBC toggle in header
  const periodic = createSection('section-periodic', 'Periodic', { defaultOpen: false });
  {
    // Insert PBC toggle control into the header (right-aligned)
    const hdr = periodic.section.querySelector('.panel-header');
    const pbcCtrls = document.createElement('span');
    pbcCtrls.className = 'inline-controls';
    pbcCtrls.id = 'pbcHeaderControls';
    const pbcToggle = makeToggle({
      id: 'togglePBCHeader',
      labelOn: 'On',
      labelOff: 'Off',
      title: 'Periodic Boundary Conditions',
    });
    const pbcStateLabel = document.createElement('span');
    pbcStateLabel.id = 'pbcStateLabel';
    pbcStateLabel.className = 'state-label';
    pbcStateLabel.textContent = 'Off';
    pbcToggle.classList.add('switch');
    // Prevent header toggle when interacting with toggle
    ['click', 'pointerdown', 'change'].forEach((ev) =>
      pbcToggle.addEventListener(ev, (e) => e.stopPropagation())
    );
    pbcCtrls.append(pbcStateLabel, pbcToggle);
    hdr.appendChild(pbcCtrls);

    // Full parameters: a, b, c (Å), alpha, beta, gamma (°)
    const monoBlock = document.createElement('div');
    monoBlock.className = 'block';
    const monoLabel = document.createElement('div');
    monoLabel.className = 'form-label';
    monoLabel.textContent = 'Unit cell (a,b,c,α,β,γ)';
    const makeValueDisplay = (id, placeholder) => {
      const span = document.createElement('span');
      span.id = id;
      span.title = placeholder || '';
      span.classList.add('mono');
      span.style.display = 'inline-block';
      span.style.width = '8ch';
      span.style.textAlign = 'right';
      span.style.background = '#171717';
      span.style.border = '1px solid var(--border)';
      span.style.borderRadius = '8px';
      span.style.padding = '6px 8px';
      return span;
    };
    const aIn = makeValueDisplay('cellA', 'a (Å)');
    const bIn = makeValueDisplay('cellB', 'b (Å)');
    const cIn = makeValueDisplay('cellC', 'c (Å)');
    const alphaIn = makeValueDisplay('cellAlpha', 'α (°)');
    const betaIn = makeValueDisplay('cellBeta', 'β (°)');
    const gammaIn = makeValueDisplay('cellGamma', 'γ (°)');
    // Nudge buttons (+/-) with press-and-hold auto-repeat
    function makeNudgers(valueEl, { idPrefix, step = 0.1, min = 0.01, max = 999, onApply }) {
      const minus = document.createElement('button');
      minus.className = 'btn';
      minus.textContent = '-';
      minus.title = 'Decrease';
      minus.id = idPrefix + 'Minus';
      const plus = document.createElement('button');
      plus.className = 'btn';
      plus.textContent = '+';
      plus.title = 'Increase';
      plus.id = idPrefix + 'Plus';
      const DECIMALS = step < 1 ? 2 : 0;
      const clamp = (v) => Math.max(min, Math.min(max, v));
      const setVal = (v) => {
        const nv = clamp(v);
        valueEl.textContent = nv.toFixed(DECIMALS);
        if (typeof onApply === 'function') onApply();
      };
      const bump = (dir) => {
        if (plus.disabled) return;
        const cur = parseFloat(valueEl.textContent);
        const base = Number.isFinite(cur) ? cur : dir > 0 ? step : step;
        setVal(base + dir * step);
      };
      // Click does one step
      minus.addEventListener('click', (e) => {
        e.preventDefault();
        bump(-1);
      });
      plus.addEventListener('click', (e) => {
        e.preventDefault();
        bump(+1);
      });
      // Press-and-hold: pointerdown starts, pointerup/cancel ends
      function attachHold(btn, dir) {
        let tid = null;
        let started = false;
        let first = true;
        const start = (e) => {
          e.preventDefault();
          if (plus.disabled) return;
          if (started) return;
          started = true;
          bump(dir);
          // after initial, repeat at interval
          tid = setInterval(() => bump(dir), 80);
          window.addEventListener('pointerup', stop, { once: true });
          window.addEventListener('pointercancel', stop, { once: true });
          btn.addEventListener('pointerleave', stop, { once: true });
        };
        const stop = () => {
          if (tid) {
            clearInterval(tid);
            tid = null;
          }
          started = false;
        };
        btn.addEventListener('pointerdown', start);
      }
      attachHold(minus, -1);
      attachHold(plus, +1);
      return { minus, plus };
    }
    // No Apply button: values are controlled only via +/- and auto-apply
    // Arrange inputs with their nudgers as one row per parameter: Label: [value] (+) (-)
    // Apply full cell parameters whenever any of a,b,c,alpha,beta,gamma changes
    const aN = makeNudgers(aIn, {
      idPrefix: 'cellA',
      step: 0.1,
      min: 0.01,
      max: 999,
      onApply: () => applyCellParams(),
    });
    const bN = makeNudgers(bIn, {
      idPrefix: 'cellB',
      step: 0.1,
      min: 0.01,
      max: 999,
      onApply: () => applyCellParams(),
    });
    const cN = makeNudgers(cIn, {
      idPrefix: 'cellC',
      step: 0.1,
      min: 0.01,
      max: 999,
      onApply: () => applyCellParams(),
    });
    const alphaN = makeNudgers(alphaIn, {
      idPrefix: 'cellAlpha',
      step: 1,
      min: 1,
      max: 179,
      onApply: () => applyCellParams(),
    });
    const betaN = makeNudgers(betaIn, {
      idPrefix: 'cellBeta',
      step: 1,
      min: 1,
      max: 179,
      onApply: () => applyCellParams(),
    });
    const gammaN = makeNudgers(gammaIn, {
      idPrefix: 'cellGamma',
      step: 1,
      min: 1,
      max: 179,
      onApply: () => applyCellParams(),
    });
    function makeParamRow(labelText, input, nudgers) {
      const row = document.createElement('div');
      row.className = 'row';
      const lbl = document.createElement('span');
      lbl.className = 'form-label';
      lbl.textContent = labelText;
      // Fix label width for alignment
      lbl.style.minWidth = '20px';
      lbl.style.display = 'inline-block';
      // Order: value then (+) then (-)
      row.append(lbl, input, nudgers.plus, nudgers.minus);
      return row;
    }
    const aRow = makeParamRow('A:', aIn, aN);
    const bRow = makeParamRow('B:', bIn, bN);
    const cRow = makeParamRow('C:', cIn, cN);
    const alphaRow = makeParamRow('α:', alphaIn, alphaN);
    const betaRow = makeParamRow('β:', betaIn, betaN);
    const gammaRow = makeParamRow('γ:', gammaIn, gammaN);
    // Append rows and apply button (its own row)
    monoBlock.append(monoLabel, aRow, bRow, cRow, alphaRow, betaRow, gammaRow);
    periodic.content.appendChild(monoBlock);

    function getApi() {
      try {
        return window.viewerApi || window._viewer;
      } catch {
        return null;
      }
    }
    function len(v) {
      return Math.hypot(v?.x || 0, v?.y || 0, v?.z || 0);
    }
    function dot(u, v) {
      return (u?.x || 0) * (v?.x || 0) + (u?.y || 0) * (v?.y || 0) + (u?.z || 0) * (v?.z || 0);
    }
    function clamp(v, lo, hi) {
      return Math.max(lo, Math.min(hi, v));
    }
    function rad2deg(r) {
      return (r * 180) / Math.PI;
    }
    function deg2rad(d) {
      return (d * Math.PI) / 180;
    }
    function extractParams(cell) {
      if (!cell) return { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 };
      const p = getCellParameters(cell) || { a: 1, b: 1, c: 1, alpha: 90, beta: 90, gamma: 90 };
      return p;
    }
    function prefill() {
      const api = getApi();
      const st = api && api.state;
      const c = st && st.cell;
      const p = extractParams(c);
      aIn.textContent = p.a.toFixed(2);
      bIn.textContent = p.b.toFixed(2);
      cIn.textContent = p.c.toFixed(2);
      alphaIn.textContent = p.alpha.toFixed(2);
      betaIn.textContent = p.beta.toFixed(2);
      gammaIn.textContent = p.gamma.toFixed(2);
      const enabled = !!(st && st.showCell && c && c.enabled);
      for (const el of [
        aN.minus,
        aN.plus,
        bN.minus,
        bN.plus,
        cN.minus,
        cN.plus,
        alphaN.minus,
        alphaN.plus,
        betaN.minus,
        betaN.plus,
        gammaN.minus,
        gammaN.plus,
      ])
        el.disabled = !enabled;
      // Reflect toggle state
      const on = !!(st && st.showCell && c && c.enabled);
      try {
        pbcToggle.setAttribute('data-on', String(on));
        pbcToggle.setAttribute('aria-checked', String(on));
        pbcToggle.textContent = on ? 'On' : 'Off';
        if (pbcStateLabel) pbcStateLabel.textContent = on ? 'On' : 'Off';
      } catch {}
    }
    function applyCellParams() {
      const api = getApi();
      const st = api && api.state;
      if (!st) return;
      let aL = parseFloat(aIn.textContent),
        bL = parseFloat(bIn.textContent),
        cL = parseFloat(cIn.textContent);
      let alphaD = parseFloat(alphaIn.textContent),
        betaD = parseFloat(betaIn.textContent),
        gammaD = parseFloat(gammaIn.textContent);
      if (!Number.isFinite(aL) || aL <= 0) aL = 1;
      if (!Number.isFinite(bL) || bL <= 0) bL = 1;
      if (!Number.isFinite(cL) || cL <= 0) cL = 1;
      if (!Number.isFinite(alphaD)) alphaD = 90;
      if (!Number.isFinite(betaD)) betaD = 90;
      if (!Number.isFinite(gammaD)) gammaD = 90;
      alphaD = clamp(alphaD, 1, 179);
      betaD = clamp(betaD, 1, 179);
      gammaD = clamp(gammaD, 1, 179);
      const cell = buildCellFromParameters({
        a: aL,
        b: bL,
        c: cL,
        alpha: alphaD,
        beta: betaD,
        gamma: gammaD,
      });
      const originOffset =
        st.cell && st.cell.originOffset ? { ...st.cell.originOffset } : { x: 0, y: 0, z: 0 };
      st.cell = { ...cell, originOffset, enabled: true };
      st.markCellChanged && st.markCellChanged();
    }
    // No Apply or Enter binding; auto-apply occurs on each nudge

    // Toggle behavior to set PBC On/Off
    function setPBC(wantOn) {
      try {
        const api = getApi();
        if (!api) return;
        const st = api.state;
        if (!st) return;
        const cur = !!st.showCell;
        if (cur === wantOn) return;
        const fn = st.toggleCellVisibilityEnhanced || st.toggleCellVisibility;
        if (typeof fn === 'function') fn.call(st);
        // Keep ghosts in sync with cell visibility
        const wantGhosts = !!st.showCell;
        if (!!st.showGhostCells !== wantGhosts && typeof st.toggleGhostCells === 'function') {
          st.toggleGhostCells();
        }
        const legacy = document.getElementById('btnCell');
        if (legacy) legacy.textContent = st.showCell || st.showGhostCells ? 'on' : 'off';
        prefill();
      } catch {}
    }
    pbcToggle.addEventListener('click', () => {
      const on = pbcToggle.getAttribute('data-on') === 'true';
      // makeToggle already flipped it; we want to apply the new state
      setPBC(on);
      try {
        if (pbcStateLabel) pbcStateLabel.textContent = on ? 'On' : 'Off';
      } catch {}
    });
    // Initial fill and react to cell changes (including external)
    try {
      prefill();
    } catch {}
    try {
      const api = getApi();
      api && api.state && api.state.bus && api.state.bus.on('cellChanged', prefill);
    } catch {}

  }

  // XR fixed widget (top-right)
  {
    if (!document.getElementById('xrModeWidget')) {
      const box = document.createElement('div');
      box.id = 'xrModeWidget';
      const label = document.createElement('span');
      label.className = 'form-label';
      label.textContent = 'XR:';
      const sel = document.createElement('select');
      sel.id = 'xrModeSelect';
      sel.innerHTML =
        '<option value="none">Off</option><option value="vr">VR</option><option value="ar">AR</option>';
      function getViewer() {
        try {
          return window.viewerApi || window._viewer;
        } catch {
          return null;
        }
      }
      sel.addEventListener('change', async (ev) => {
        const v = (ev && ev.target && ev.target.value) || sel.value;
        const api = getViewer();
        const vr = api && api.vr;
        const doSwitch = async (mode) => {
          try {
            if (vr && typeof vr.switchXR === 'function') return await vr.switchXR(mode);
            if (mode === 'none') {
              if (vr?.exitXR) return await vr.exitXR();
              if (vr?.exitVR) return await vr.exitVR();
              return false;
            }
            if (mode === 'vr') {
              if (vr?.enterVR) return await vr.enterVR();
              return false;
            }
            if (mode === 'ar') {
              if (vr?.enterAR) return await vr.enterAR();
              return false;
            }
            return false;
          } catch {
            return false;
          }
        };
        if (!vr) {
          try {
            sel.value = 'none';
          } catch {}
          return;
        }
        await doSwitch(v);
      });
      box.append(label, sel);
      document.body.appendChild(box);
    }
  }

  // Append all sections (Selection under Live Metrics)
  // Build panel: sections only (XR lives in a fixed widget)
  panel.append(live.section, selSec.section, sim.section, periodic.section, sys.section);
  host.appendChild(panel);

  // Add a persistent reset button in bottom-right that resets the viewer state without reloading the page
  try {
    if (!document.getElementById('resetAllBtn')) {
      const resetBtn = document.createElement('button');
      resetBtn.id = 'resetAllBtn';
      resetBtn.textContent = 'Reset';
      resetBtn.title = 'Reset simulation to initial state (no page reload)';
      Object.assign(resetBtn.style, {
        position: 'fixed',
        right: '12px',
        bottom: '12px',
        zIndex: 1000,
        background: 'rgba(30,38,48,0.85)',
        color: '#d8e6f3',
        border: '1px solid rgba(255,255,255,0.1)',
        borderRadius: '10px',
        padding: '10px 14px',
        cursor: 'pointer',
      });
      resetBtn.addEventListener('click', () => {
        try {
          // Prefer in-app reset if available
          const api = typeof window !== 'undefined' ? window.viewerApi || null : null;
          if (api && typeof api.resetToInitialPositions === 'function') {
            // Prevent double clicks during reset
            resetBtn.disabled = true;
            const prevText = resetBtn.textContent;
            resetBtn.textContent = 'Resetting…';
            Promise.resolve(api.resetToInitialPositions()).finally(() => {
              resetBtn.disabled = false;
              resetBtn.textContent = prevText;
            });
          } else {
            // Fallback to page reload if API missing
            const href = window.location.href; // preserve current query (mol, smiles, etc.)
            if (typeof window.location.assign === 'function') window.location.assign(href);
            else window.location.href = href;
          }
        } catch {}
      });
      document.body.appendChild(resetBtn);
    }
  } catch {}

  // --- Mobile top bar with three tabs (Live Metrics, Simulation, System) ---
  const topBar = document.createElement('div');
  topBar.id = 'mobileTopBar';
  const tabsRow = document.createElement('div');
  tabsRow.className = 'tabs-row';
  const tabs = [
    { id: 'live', label: 'Live Metrics', sectionId: 'section-live-stats' },
    { id: 'simulation', label: 'Simulation', sectionId: 'section-simulation' },
    { id: 'periodic', label: 'Periodic', sectionId: 'section-periodic' },
    { id: 'system', label: 'System', sectionId: 'section-system' },
  ];
  const tabEls = new Map();
  tabs.forEach((t) => {
    const b = document.createElement('button');
    b.className = 'tab';
    b.type = 'button';
    b.id = `mobileTab-${t.id}`;
    b.textContent = t.label;
    b.setAttribute('data-active', 'false');
    b.addEventListener('click', () => toggleMobileTab(t.sectionId, b));
    tabEls.set(t.sectionId, b);
    tabsRow.appendChild(b);
  });
  const mobileSheet = document.createElement('div');
  mobileSheet.className = 'mobile-sheet';
  topBar.appendChild(tabsRow);
  topBar.appendChild(mobileSheet);
  host.appendChild(topBar);

  // Map section id -> content node (stable reference even if reparented)
  const sectionContentMap = new Map([
    ['section-live-stats', live.content],
    ['section-selection', selSec.content],
    ['section-simulation', sim.content],
    ['section-periodic', periodic.content],
    ['section-system', sys.content],
  ]);

  // Track original parent/nextSibling to restore reparented content
  const restoreInfo = new Map(); // node -> { parent, nextSibling, collapsed, headerExpanded }

  function restoreNode(node) {
    try {
      const info = restoreInfo.get(node);
      if (!info) return;
      if (node.parentElement === mobileSheet) {
        if (info.nextSibling && info.parent === info.nextSibling.parentNode) {
          info.parent.insertBefore(node, info.nextSibling);
        } else if (info.parent) {
          info.parent.appendChild(node);
        }
      }
      if (info.collapsed) node.setAttribute('data-collapsed', 'true');
      else node.removeAttribute('data-collapsed');
      const secEl = node.parentElement;
      const hdr = secEl && secEl.querySelector && secEl.querySelector('.panel-header');
      if (hdr) hdr.setAttribute('aria-expanded', info.collapsed ? 'false' : 'true');
      restoreInfo.delete(node);
    } catch {}
  }

  function mountToBar(sectionId) {
    try {
      const node = sectionContentMap.get(sectionId);
      if (!node) return;
      // Save original placement once
      if (!restoreInfo.has(node)) {
        restoreInfo.set(node, {
          parent: node.parentElement,
          nextSibling: node.nextSibling,
          collapsed: node.getAttribute('data-collapsed') === 'true',
          headerExpanded:
            (node.previousElementSibling &&
              node.previousElementSibling.getAttribute('aria-expanded')) ||
            'false',
        });
      }
      // Ensure expanded state for visibility in the bar
      try {
        node.removeAttribute('data-collapsed');
      } catch {}
      const hdr = document.querySelector(`#${sectionId} .panel-header`);
      if (hdr) hdr.setAttribute('aria-expanded', 'true');
      // Move node into the mobile sheet
      if (node.parentElement !== mobileSheet) mobileSheet.appendChild(node);
    } catch {}
  }

  function restoreAllToPanel() {
    try {
      restoreInfo.forEach((info, node) => {
        if (!info || !node) return;
        if (node.parentElement === mobileSheet) {
          if (info.nextSibling && info.parent === info.nextSibling.parentNode) {
            info.parent.insertBefore(node, info.nextSibling);
          } else if (info.parent) {
            info.parent.appendChild(node);
          }
        }
        // Restore collapsed + header state
        if (info.collapsed) node.setAttribute('data-collapsed', 'true');
        else node.removeAttribute('data-collapsed');
        const secEl = node.parentElement; // should now be the section
        const hdr = secEl && secEl.querySelector && secEl.querySelector('.panel-header');
        if (hdr) hdr.setAttribute('aria-expanded', info.collapsed ? 'false' : 'true');
      });
      restoreInfo.clear();
    } catch {}
  }

  function makeExclusiveInSheet(sectionId) {
    try {
      const want = sectionContentMap.get(sectionId) || null;
      const children = Array.from(mobileSheet.children || []);
      for (const ch of children) {
        if (ch !== want) restoreNode(ch);
      }
      mountToBar(sectionId);
    } catch {}
  }

  function updateTopBarVisibility(open) {
    try {
      if (open) {
        topBar.setAttribute('data-open', 'true');
      } else {
        topBar.removeAttribute('data-open');
      }
      // Tabs remain visible at all times; ensure no inline hiding remains
      const buttons = topBar.querySelectorAll('button.tab');
      if (buttons && buttons.length) {
        buttons.forEach((b) => {
          b.style.display = '';
        });
      }
      // Control the in-bar sheet visibility
      if (mobileSheet) mobileSheet.style.display = open ? 'block' : 'none';
    } catch {}
  }

  // --- Responsive behavior (desktop vs mobile) ---
  function isMobileViewport() {
    try {
      const w = typeof window !== 'undefined' ? window.innerWidth || 0 : 0;
      const mqNarrow =
        typeof window !== 'undefined' &&
        window.matchMedia &&
        window.matchMedia('(max-width: 800px)').matches;
      const mqCoarse =
        typeof window !== 'undefined' &&
        window.matchMedia &&
        window.matchMedia('(pointer: coarse)').matches;
      return (w > 0 && w <= 800) || mqNarrow || mqCoarse;
    } catch {
      return false;
    }
  }
  function applyLayout() {
    const mobile = isMobileViewport();
    if (mobile) {
      // Mobile: hide left panel; show the top bar and render content inside it
      try {
        panel.setAttribute('data-mode', 'mobile');
      } catch {}
      try {
        topBar.style.display = 'flex';
      } catch {}
      const open = panel.getAttribute('data-mobile-open') === 'true';
      // Keep the left panel hidden on mobile; content lives in the bar
      panel.style.display = 'none';
      // If open, ensure only the active section content is mounted into the bar
      if (open) {
        const current = panel.querySelector('.panel-section[data-mobile-active="true"]');
        if (current) makeExclusiveInSheet(current.id);
      } else {
        restoreAllToPanel();
      }
      // Reflect state on the top bar
      updateTopBarVisibility(open);
    } else {
      // Desktop: always show left panel with all sections; hide the top bar
      try {
        panel.setAttribute('data-mode', 'desktop');
      } catch {}
      try {
        topBar.style.display = 'none';
      } catch {}
      panel.style.display = 'block';
      // Ensure mobile-only attributes are reset so desktop shows all sections
      panel.removeAttribute('data-mobile-open');
      const secs = panel.querySelectorAll('.panel-section');
      secs.forEach((sec) => sec.removeAttribute('data-mobile-active'));
      // Ensure any reparented content is restored to original sections
      restoreAllToPanel();
      clearActiveTabs();
      updateTopBarVisibility(false);
    }
  }

  function setOnlySectionActive(sectionId) {
    const secs = panel.querySelectorAll('.panel-section');
    secs.forEach((sec) => {
      if (sec.id === sectionId) {
        sec.setAttribute('data-mobile-active', 'true');
        const content = sectionContentMap.get(sectionId) || sec.querySelector('.panel-content');
        const header = sec.querySelector('.panel-header');
        if (content) content.removeAttribute('data-collapsed');
        if (header) header.setAttribute('aria-expanded', 'true');
      } else {
        sec.removeAttribute('data-mobile-active');
      }
    });
  }
  function clearActiveTabs() {
    tabEls.forEach((btn) => btn.setAttribute('data-active', 'false'));
  }
  function toggleMobileTab(sectionId, btn) {
    const isOpen = panel.getAttribute('data-mobile-open') === 'true';
    const current = panel.querySelector('.panel-section[data-mobile-active="true"]');
    const isSame = !!(current && current.id === sectionId);
    if (isOpen && isSame) {
      panel.removeAttribute('data-mobile-open');
      const secs = panel.querySelectorAll('.panel-section');
      secs.forEach((sec) => sec.removeAttribute('data-mobile-active'));
      clearActiveTabs();
      // On close: restore content to panel and show only tabs
      restoreAllToPanel();
      updateTopBarVisibility(false);
      return;
    }
    panel.setAttribute('data-mobile-open', 'true');
    setOnlySectionActive(sectionId);
    clearActiveTabs();
    if (btn) btn.setAttribute('data-active', 'true');
    // Render content inside the bar; keep the side panel hidden on mobile
    makeExclusiveInSheet(sectionId);
    updateTopBarVisibility(true);
  }

  // Legacy hidden controls for test compatibility (text reflects on/off)
  const legacyIds = [
    'btnRelax',
    'btnRelaxRun',
    'btnMD',
    'btnMDRun',
    'btnToggleForces',
    'toggleEnergyPlot',
    'btnCell',
  ];
  const legacyBox = document.createElement('div');
  legacyBox.style.display = 'none';
  legacyBox.id = 'legacyControlsHidden';
  for (const id of legacyIds) {
    const b = document.createElement('button');
    b.id = id;
    b.textContent = 'run';
    legacyBox.appendChild(b);
  }
  panel.appendChild(legacyBox);

  // Wire legacy "toggleEnergyPlot" to show/hide the energyPlot for backward-compat tests
  try {
    const plot = panel.querySelector('#energyPlot');
    const btn = legacyBox.querySelector('#toggleEnergyPlot');
    if (plot && btn) {
      // default visible
      if (!plot.style.display) plot.style.display = 'block';
      btn.addEventListener('click', () => {
        const cur = plot.style.display;
        plot.style.display = cur === 'none' ? 'block' : 'none';
      });
    }
  } catch {}

  // Initial layout and resize handling
  try {
    applyLayout();
  } catch {}
  try {
    if (typeof window !== 'undefined') {
      window.addEventListener('resize', applyLayout, { passive: true });
    }
  } catch {}

  return {
    panelEl: panel,
    containers: {
      live: live.content,
      selection: selSec.content,
      simulation: sim.content,
      system: sys.content,
    },
  };
}
