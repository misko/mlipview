// Builds the left-side desktop UI panel with collapsible sections.
// Sections: Live Metrics, Simulation, System, Rendering, XR
// Returns { panelEl, containers: { live, simulation, system, rendering, xr } }

import { initTemperatureSlider } from './temperatureSlider.js';
import { installMoleculeSelector } from './moleculeSelect.js';
import { elInfo } from '../elements.js';
import { defaultMassForZ } from '../physics/sim-model.js';

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
      #energyYLabel { position: absolute; left: 8px; top: 50%; transform: translateY(-50%) rotate(-90deg); transform-origin: center; font-size: 10px; color: var(--muted); white-space: nowrap; }
      #energyXLabel { position: absolute; left: 50%; bottom: 4px; transform: translateX(-50%); font-size: 10px; color: var(--muted); white-space: nowrap; }

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
    const yLabel = document.createElement('div'); yLabel.id = 'energyYLabel'; yLabel.textContent = 'Energy';
    const xLabel = document.createElement('div'); xLabel.id = 'energyXLabel'; xLabel.textContent = 'Time(step)';
    plot.append(canvas, yLabel, xLabel);
    live.content.appendChild(plot);
  }

  // Selection (open) — shows selected atom/bond info and mini periodic table
  // Note: Periodic table cells are intentionally non-clickable for v1; this is a purely informative view.
  // VR/AR parity plan (future): Mirror this info in XR by rendering a small Babylon GUI (AdvancedDynamicTexture)
  // panel pinned to the controller ray or wrist anchor. Subscribe to the same state bus events (selectionChanged,
  // positionsChanged) and render simple GUI controls (text blocks + colored circles) for the spheres and labels.
  // For element highlighting, draw a compact grid of TextBlocks/images; apply highlight style (border/emissive) when selected.
  // Wire via viewerApi.vr.* helpers and ensure the panel is disabled in desktop to avoid duplication.
  const selSec = createSection('section-selection', 'Selection', { defaultOpen: true });
  {
    // Styles specific to selection section (only once)
    const styleId = 'selectionSectionStyles';
    if (!document.getElementById(styleId)) {
      const s = document.createElement('style');
      s.id = styleId;
      s.textContent = `
        /* Selection section */
        #section-selection .panel-content { padding-top:6px; padding-bottom:8px; }
        #section-selection .mini-spheres { display:flex; align-items:center; gap:8px; height:46px; }
        #section-selection .atom-sphere { width:30px; height:30px; border-radius:50%; border:1px solid rgba(255,255,255,0.15); box-shadow: inset 0 0 6px rgba(0,0,0,0.5); }
  #section-selection .info { font-size:12px; color: var(--text); line-height:1.4; }
  #section-selection .info #selPosition { font-size: 11px; color: var(--text); }
        #section-selection .info .label { color: var(--muted); margin-right:6px; }
        #section-selection .grid { margin-top:8px; display:grid; grid-template-columns: repeat(10, 1fr); gap:4px; }
        #section-selection .pt-el { text-align:center; padding:0; border-radius:2px; border:1px solid var(--border); font-size:10px; color: var(--text); background:#171717; user-select:none; box-sizing:border-box; }
        #section-selection .pt-el.highlight { border-color: var(--accent); box-shadow: 0 0 0 1px var(--accent) inset; }
      `;
      document.head.appendChild(s);
    }

    // Color palette for common elements (CSS colors; approximate Babylon colors)
    const EL_COLORS = {
      H: '#f2f7ff', C:'#383d47', N:'#6290f7', O:'#f25454', S:'#f7d64a', Cl:'#7bd47b',
      F:'#8ef06f', P:'#ffa94d', Br:'#a3552e', I:'#8a2be2', Si:'#f7786b', Fe:'#c44e4e',
      Co:'#cc6699', Ni:'#6fbf73', Cu:'#b87333', Zn:'#9bb5c1', Ag:'#c0c0c0', Au:'#ffd700',
      Na:'#abdee6', K:'#ffc8dd', Ca:'#ffe1a8', Li:'#cdeac0', He:'#c7d2fe', Ne:'#a7f3d0'
    };
    const EL_NAMES = {
      H:'Hydrogen', C:'Carbon', N:'Nitrogen', O:'Oxygen', S:'Sulfur', Cl:'Chlorine',
      F:'Fluorine', P:'Phosphorus', Br:'Bromine', I:'Iodine', Si:'Silicon', Fe:'Iron',
      Co:'Cobalt', Ni:'Nickel', Cu:'Copper', Zn:'Zinc', Ag:'Silver', Au:'Gold',
      Na:'Sodium', K:'Potassium', Ca:'Calcium', Li:'Lithium', He:'Helium', Ne:'Neon'
    };
    const SYMBOL_TO_Z = { H:1, He:2, Li:3, Be:4, B:5, C:6, N:7, O:8, F:9, Ne:10, Na:11, Mg:12, Al:13, Si:14, P:15, S:16, Cl:17, Ar:18, K:19, Ca:20 };
    const Z_TO_SYMBOL = Object.fromEntries(Object.entries(SYMBOL_TO_Z).map(([k,v])=>[v,k]));
    function getSymbol(el){
      if (!el) return 'X';
      if (typeof el === 'string') return el;
      if (typeof el === 'number') return Z_TO_SYMBOL[el] || 'X';
      return el.symbol || el.sym || el.S || (Z_TO_SYMBOL[el.Z||el.z||el.atomicNumber] || 'X');
    }
  function fmtPos(p){ if(!p) return '(-,-,-)'; return `(${p.x.toFixed(1)}, ${p.y.toFixed(1)}, ${p.z.toFixed(1)})`; }
    function dist(a,b){ const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return Math.sqrt(dx*dx+dy*dy+dz*dz); }

    // UI structure
    const spheres = document.createElement('div'); spheres.className = 'mini-spheres';
    const sphereA = document.createElement('div'); sphereA.className = 'atom-sphere'; sphereA.style.display='none'; sphereA.id='selSphereA';
    const sphereB = document.createElement('div'); sphereB.className = 'atom-sphere'; sphereB.style.display='none'; sphereB.id='selSphereB';
    spheres.append(sphereA, sphereB);

    const info = document.createElement('div'); info.className = 'info';
    info.innerHTML = `
      <div><span class="label">Element:</span><span id="selElementName">—</span></div>
      <div><span class="label">Position:</span><span id="selPosition">(-,-,-)</span></div>
      <div><span class="label">Atomic weight (amu):</span><span id="selAtomicWeight">—</span></div>
      <div><span class="label">vdW radius (Å):</span><span id="selVdw">—</span></div>
      <div id="bondInfoRow"><span class="label">Bond length:</span><span id="bondLength">N/A</span></div>
      <div id="rotateRow"><span class="label">Rotate bond:</span>
        <span id="rotateNA">N/A</span>
        <span id="rotateBtns" style="display:none; gap:6px;">
          <button id="bondRotMinus" class="btn" title="Rotate -">-</button>
          <button id="bondRotPlus" class="btn" title="Rotate +">+</button>
        </span>
      </div>
    `;

    // Full periodic table layout as a fixed HTML table (18 columns).
    // Lanthanides/Actinides shown as separate indented rows. Using a table ensures consistent cell sizing and
    // visible grid lines for rows/columns.
    const mini = document.createElement('table'); mini.id = 'miniPeriodic'; mini.setAttribute('role','grid');
    const miniBody = document.createElement('tbody'); mini.appendChild(miniBody);
    // Table structure styles
    const ptStyleId = 'miniPeriodicStyles';
    if (!document.getElementById(ptStyleId)) {
      const st = document.createElement('style'); st.id = ptStyleId; st.textContent = `
        #miniPeriodic { width:100%; table-layout: fixed; border-collapse: collapse; margin-top:6px; background:#232831; border:1px solid rgba(255,255,255,0.08); }
        #miniPeriodic tr { height:18px; }
        #miniPeriodic td.pt-el { text-align:center; padding:0; height:18px; line-height:18px; font-size:10px; color: var(--text); background:#1b2027; border:1px solid rgba(255,255,255,0.08); box-sizing:border-box; }
        #miniPeriodic td.pt-el.empty { background:#232831; border:none; }
        #miniPeriodic td.pt-el.highlight { border-color: var(--accent); box-shadow: inset 0 0 0 1px var(--accent); }
      `; document.head.appendChild(st);
    }
    // Periods 1..7
    const PERIOD_ROWS = [
      // 1
      ['H','','','','','','','','','','','','','','','','','He'],
      // 2
      ['Li','Be','','','','','','','','','','','B','C','N','O','F','Ne'],
      // 3
      ['Na','Mg','','','','','','','','','','','Al','Si','P','S','Cl','Ar'],
      // 4
      ['K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'],
      // 5
      ['Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe'],
      // 6 (La is shown in group 3; rest of lanthanides moved to separate row)
      ['Cs','Ba','La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'],
      // 7 (Ac similarly)
      ['Fr','Ra','Ac','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']
    ];
    const LANTH = ['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'];
    const ACTIN = ['Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'];
    // Symbol->full name map (118). For brevity, keep names concise.
    const SYMBOL_TO_NAME = {
      H:'Hydrogen', He:'Helium', Li:'Lithium', Be:'Beryllium', B:'Boron', C:'Carbon', N:'Nitrogen', O:'Oxygen', F:'Fluorine', Ne:'Neon',
      Na:'Sodium', Mg:'Magnesium', Al:'Aluminum', Si:'Silicon', P:'Phosphorus', S:'Sulfur', Cl:'Chlorine', Ar:'Argon',
      K:'Potassium', Ca:'Calcium', Sc:'Scandium', Ti:'Titanium', V:'Vanadium', Cr:'Chromium', Mn:'Manganese', Fe:'Iron', Co:'Cobalt', Ni:'Nickel', Cu:'Copper', Zn:'Zinc', Ga:'Gallium', Ge:'Germanium', As:'Arsenic', Se:'Selenium', Br:'Bromine', Kr:'Krypton',
      Rb:'Rubidium', Sr:'Strontium', Y:'Yttrium', Zr:'Zirconium', Nb:'Niobium', Mo:'Molybdenum', Tc:'Technetium', Ru:'Ruthenium', Rh:'Rhodium', Pd:'Palladium', Ag:'Silver', Cd:'Cadmium', In:'Indium', Sn:'Tin', Sb:'Antimony', Te:'Tellurium', I:'Iodine', Xe:'Xenon',
      Cs:'Cesium', Ba:'Barium', La:'Lanthanum', Ce:'Cerium', Pr:'Praseodymium', Nd:'Neodymium', Pm:'Promethium', Sm:'Samarium', Eu:'Europium', Gd:'Gadolinium', Tb:'Terbium', Dy:'Dysprosium', Ho:'Holmium', Er:'Erbium', Tm:'Thulium', Yb:'Ytterbium', Lu:'Lutetium',
      Hf:'Hafnium', Ta:'Tantalum', W:'Tungsten', Re:'Rhenium', Os:'Osmium', Ir:'Iridium', Pt:'Platinum', Au:'Gold', Hg:'Mercury', Tl:'Thallium', Pb:'Lead', Bi:'Bismuth', Po:'Polonium', At:'Astatine', Rn:'Radon',
      Fr:'Francium', Ra:'Radium', Ac:'Actinium', Th:'Thorium', Pa:'Protactinium', U:'Uranium', Np:'Neptunium', Pu:'Plutonium', Am:'Americium', Cm:'Curium', Bk:'Berkelium', Cf:'Californium', Es:'Einsteinium', Fm:'Fermium', Md:'Mendelevium', No:'Nobelium', Lr:'Lawrencium',
      Rf:'Rutherfordium', Db:'Dubnium', Sg:'Seaborgium', Bh:'Bohrium', Hs:'Hassium', Mt:'Meitnerium', Ds:'Darmstadtium', Rg:'Roentgenium', Cn:'Copernicium', Nh:'Nihonium', Fl:'Flerovium', Mc:'Moscovium', Lv:'Livermorium', Ts:'Tennessine', Og:'Oganesson'
    };

    // Approximate standard atomic weights (amu) keyed by symbol (for UI display).
    // Values are representative; for synthetic/unstable elements, common isotopic masses are used.
    const MASS_BY_SYMBOL = {
      H:1.008, He:4.0026, Li:6.94, Be:9.0122, B:10.81, C:12.011, N:14.007, O:15.999, F:18.998, Ne:20.180,
      Na:22.990, Mg:24.305, Al:26.982, Si:28.085, P:30.974, S:32.06, Cl:35.45, Ar:39.948,
      K:39.0983, Ca:40.078, Sc:44.9559, Ti:47.867, V:50.9415, Cr:51.9961, Mn:54.938, Fe:55.845, Co:58.933, Ni:58.693, Cu:63.546, Zn:65.38, Ga:69.723, Ge:72.630, As:74.922, Se:78.971, Br:79.904, Kr:83.798,
      Rb:85.468, Sr:87.62, Y:88.906, Zr:91.224, Nb:92.906, Mo:95.95, Tc:98, Ru:101.07, Rh:102.905, Pd:106.42, Ag:107.868, Cd:112.414, In:114.818, Sn:118.710, Sb:121.760, Te:127.60, I:126.904, Xe:131.293,
      Cs:132.905, Ba:137.327, La:138.905, Ce:140.116, Pr:140.908, Nd:144.242, Pm:145, Sm:150.36, Eu:151.964, Gd:157.25, Tb:158.925, Dy:162.500, Ho:164.930, Er:167.259, Tm:168.934, Yb:173.045, Lu:174.967,
      Hf:178.49, Ta:180.948, W:183.84, Re:186.207, Os:190.23, Ir:192.217, Pt:195.084, Au:196.96657, Hg:200.592, Tl:204.38, Pb:207.2, Bi:208.980, Po:209, At:210, Rn:222,
      Fr:223, Ra:226, Ac:227, Th:232.0377, Pa:231.0359, U:238.0289, Np:237, Pu:244, Am:243, Cm:247, Bk:247, Cf:251, Es:252, Fm:257, Md:258, No:259, Lr:266,
      Rf:267, Db:268, Sg:269, Bh:270, Hs:277, Mt:278, Ds:281, Rg:282, Cn:285, Nh:286, Fl:289, Mc:290, Lv:293, Ts:294, Og:294
    };

    // Build DOM rows
    function makeRow(symbols){
      const tr = document.createElement('tr');
      for (let c=0; c<18; c++){
        const sym = symbols[c] || '';
        const td = document.createElement('td');
        td.className = 'pt-el' + (sym? '' : ' empty');
        if (sym){
          td.setAttribute('data-symbol', sym);
          td.textContent = sym;
          td.addEventListener('click', ()=>{
            try { const api = getViewer(); api && api.selection && api.selection.clear && api.selection.clear(); } catch {}
            __overrideElementSym = sym; updateFromSelection();
          });
        } else {
          td.textContent = '';
        }
        tr.appendChild(td);
      }
      return tr;
    }
  PERIOD_ROWS.forEach(r => miniBody.appendChild(makeRow(r)));

    selSec.content.append(spheres, info, mini);

    // Live updater bound to selection changes
    function getViewer(){ try { return window.viewerApi || window._viewer; } catch { return null; } }
    function highlight(symbols){
      const cells = mini.querySelectorAll('.pt-el');
      cells.forEach(c => c.classList.remove('highlight'));
      const set = new Set((symbols||[]).filter(Boolean));
      set.forEach(sym => {
        const el = mini.querySelector(`.pt-el[data-symbol="${sym}"]`);
        if (el) el.classList.add('highlight');
      });
    }
    function setSphere(el, sym){
      if(!sym){ el.style.display='none'; return; }
      el.style.display='block';
      const color = EL_COLORS[sym] || '#9aa0a6';
      el.style.background = `radial-gradient( circle at 30% 30%, #ffffff66, transparent 40%), ${color}`;
    }
    // If a periodic table element was clicked, we set an override symbol and clear it when a real selection happens.
    let __overrideElementSym = null;
    function symbolToName(sym){ return EL_NAMES[sym] || SYMBOL_TO_NAME[sym] || sym; }
    function updateFromSelection(){
      const api = getViewer(); if(!api) return;
      const st = api.state;
      const sel = st && st.selection || { kind:null };
      const elNameNode = selSec.content.querySelector('#selElementName');
      const posNode = selSec.content.querySelector('#selPosition');
      const weightNode = selSec.content.querySelector('#selAtomicWeight');
      const vdwNode = selSec.content.querySelector('#selVdw');
      const bondRow = selSec.content.querySelector('#bondInfoRow');
      const bondLenNode = selSec.content.querySelector('#bondLength');
      const rotNA = selSec.content.querySelector('#rotateNA');
      const rotBtns = selSec.content.querySelector('#rotateBtns');
      const rotMinus = selSec.content.querySelector('#bondRotMinus');
      const rotPlus = selSec.content.querySelector('#bondRotPlus');
      // Bind handlers once
      if (!rotMinus.__bound) {
        const EPS = 0.1;
        rotMinus.addEventListener('click', ()=>{ try { const v=getViewer(); v && v.manipulation && v.manipulation.rotateBond && v.manipulation.rotateBond(-EPS); } catch{} }); rotMinus.__bound=true;
        rotPlus.addEventListener('click', ()=>{ try { const v=getViewer(); v && v.manipulation && v.manipulation.rotateBond && v.manipulation.rotateBond(EPS); } catch{} }); rotPlus.__bound=true;
      }
      function weightForSym(sym){
        if (MASS_BY_SYMBOL[sym] != null) return MASS_BY_SYMBOL[sym];
        const Z = SYMBOL_TO_Z[sym];
        if(!Z) return null;
        const m = defaultMassForZ(Z);
        return (typeof m === 'number' && isFinite(m)) ? m : null;
      }
      function vdwForSym(sym){
        const info = elInfo(sym);
        return (info && typeof info.vdw === 'number') ? info.vdw : null;
      }
      if (__overrideElementSym) {
        // Show virtual element selection
        setSphere(sphereA, __overrideElementSym); setSphere(sphereB, null);
        elNameNode.textContent = symbolToName(__overrideElementSym);
        posNode.textContent = '(-,-,-)';
        const mw = weightForSym(__overrideElementSym);
        const rv = vdwForSym(__overrideElementSym);
        weightNode.textContent = (mw!=null) ? mw.toFixed(3) : '—';
        vdwNode.textContent = (rv!=null) ? rv.toFixed(2) : '—';
        bondRow.style.display = 'block'; bondLenNode.textContent = 'N/A';
        rotNA.style.display = 'inline'; rotBtns.style.display = 'none';
        highlight([__overrideElementSym]);
      } else if (sel.kind === 'atom') {
        const idx = sel.data.index;
        const sym = getSymbol(st.elements[idx]);
        const name = symbolToName(sym);
        const pos = st.positions[idx];
        setSphere(sphereA, sym); setSphere(sphereB, null);
        elNameNode.textContent = name;
        posNode.textContent = fmtPos(pos);
        const mw = weightForSym(sym);
        const rv = vdwForSym(sym);
        weightNode.textContent = (mw!=null) ? mw.toFixed(3) : '—';
        vdwNode.textContent = (rv!=null) ? rv.toFixed(2) : '—';
        bondRow.style.display = 'block'; bondLenNode.textContent = 'N/A';
        rotNA.style.display = 'inline'; rotBtns.style.display = 'none';
        highlight([sym]);
      } else if (sel.kind === 'bond') {
        const i = sel.data.i, j = sel.data.j;
        const symA = getSymbol(st.elements[i]);
        const symB = getSymbol(st.elements[j]);
        const nameA = symbolToName(symA);
        const nameB = symbolToName(symB);
        const posA = st.positions[i];
        const posB = st.positions[j];
        setSphere(sphereA, symA); setSphere(sphereB, symB);
        elNameNode.textContent = `${nameA} – ${nameB}`;
        posNode.textContent = `${fmtPos(posA)} – ${fmtPos(posB)}`;
        const mwA = weightForSym(symA), mwB = weightForSym(symB);
        const rvA = vdwForSym(symA), rvB = vdwForSym(symB);
        weightNode.textContent = `${mwA!=null?mwA.toFixed(3):'—'} – ${mwB!=null?mwB.toFixed(3):'—'}`;
        vdwNode.textContent = `${rvA!=null?rvA.toFixed(2):'—'} – ${rvB!=null?rvB.toFixed(2):'—'}`;
        const L = dist(posA, posB);
        bondRow.style.display = 'block';
        bondLenNode.textContent = `${L.toFixed(3)} Å`;
        rotNA.style.display = 'none'; rotBtns.style.display = 'inline-flex';
        highlight([symA, symB]);
      } else {
        setSphere(sphereA, null); setSphere(sphereB, null);
        elNameNode.textContent = '—';
        posNode.textContent = '(-,-,-)';
        weightNode.textContent = '—';
        vdwNode.textContent = '—';
        bondRow.style.display = 'block'; bondLenNode.textContent = 'N/A';
        rotNA.style.display = 'inline'; rotBtns.style.display = 'none';
        highlight([]);
      }
    }
    // Initial draw and subscriptions
    updateFromSelection();
    try {
      const api = getViewer();
      if (api && api.state && api.state.bus) {
        api.state.bus.on('selectionChanged', ()=>{ __overrideElementSym = null; updateFromSelection(); });
        // Also keep positions updated in case of drag while selected
        api.state.bus.on('positionsChanged', updateFromSelection);
      }
    } catch {}
  }

  // Simulation (open) — Top utility toggles + Relax/MD toggles (side-by-side)
  const sim = createSection('section-simulation', 'Simulation', { defaultOpen: true });
  {
    // Top row: Show Forces + PBC
    const utilRow = document.createElement('div'); utilRow.className = 'row';
    // Forces toggle
    const forcesT = makeToggle({
      id: 'toggleForces',
      labelOn: 'Show Forces: On',
      labelOff: 'Show Forces: Off',
      title: 'Render per-atom force vectors',
    });
    forcesT.addEventListener('click', ()=>{
      try {
        const api = (window.viewerApi||window._viewer); if(!api) return;
        const on = forcesT.getAttribute('data-on')==='true';
        api.setForceVectorsEnabled(on);
        const legacy = document.getElementById('btnToggleForces'); if(legacy) legacy.textContent = on? 'on':'off';
      } catch{}
    });
    // PBC toggle (moved from System)
    const pbcToggle = makeToggle({
      id: 'togglePBC',
      labelOn: 'PBC: On',
      labelOff: 'PBC: Off',
      title: 'Periodic Boundary Conditions',
    });
    pbcToggle.addEventListener('click', ()=>{
      try {
        const api = (window.viewerApi||window._viewer); if(!api) return;
        const st = api.state; if(!st) return;
        const fn = (st.toggleCellVisibilityEnhanced||st.toggleCellVisibility);
        if (typeof fn === 'function') fn.call(st);
        // Keep ghost cells in sync with cell visibility
        const wantGhosts = !!st.showCell;
        if (!!st.showGhostCells !== wantGhosts && typeof st.toggleGhostCells === 'function') {
          st.toggleGhostCells();
        }
        const legacy = document.getElementById('btnCell'); if(legacy) legacy.textContent = (st.showCell || st.showGhostCells) ? 'on' : 'off';
        const status = document.getElementById('status'); if(status) status.textContent = (st.showCell || st.showGhostCells) ? 'PBC ON' : 'PBC OFF';
      } catch {}
    });
    utilRow.append(forcesT, pbcToggle);

    // Relaxation toggle (mutually exclusive with MD)
    const togglesRow = document.createElement('div'); togglesRow.className = 'row';
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
    function getApi(){ try { return window.viewerApi || window._viewer; } catch { return null; } }
    function setBtnState(btn, on){ btn.setAttribute('data-on', String(!!on)); btn.setAttribute('aria-checked', String(!!on)); btn.textContent = on ? (btn.id==='toggleMD'?'Molecular Dynamics: On':'Relaxation: On') : (btn.id==='toggleMD'?'Molecular Dynamics: Off':'Relaxation: Off'); }
    // Relaxation toggle wiring
    relaxToggle.addEventListener('click', (e)=>{
      const api = getApi(); if(!api) return;
      const on = relaxToggle.getAttribute('data-on')==='true';
      // If turning on, stop MD and start relax run (single step semantics via run loop)
      if(on){
        try { api.stopSimulation(); } catch{}
        try {
          const shortRuns = (typeof window !== 'undefined') && !!window.__MLIPVIEW_PANEL_SHORT_RUNS;
          const args = shortRuns ? { maxSteps: 2 } : {};
          api.startRelaxContinuous(args);
        } catch{}
      } else { try { api.stopSimulation(); } catch{} }
      // Mirror legacy button text convention for compatibility
      try { const legacy = document.getElementById('btnRelaxRun'); if(legacy) legacy.textContent = on? 'stop':'run'; } catch{}
    });
    // MD toggle wiring
    mdToggle.addEventListener('click', (e)=>{
      const api = getApi(); if(!api) return;
      const on = mdToggle.getAttribute('data-on')==='true';
      if(on){
        try { api.stopSimulation(); } catch{}
        try {
          const shortRuns = (typeof window !== 'undefined') && !!window.__MLIPVIEW_PANEL_SHORT_RUNS;
          const args = shortRuns ? { steps: 2 } : {};
          api.startMDContinuous(args);
        } catch{}
      } else { try { api.stopSimulation(); } catch{} }
      try {
        const legacy = document.getElementById('btnMDRun'); if(legacy) legacy.textContent = on? 'stop':'run';
        const status = document.getElementById('status'); if(status) status.textContent = on? 'MD running':'MD stopped';
      } catch{}
    });

  togglesRow.append(relaxToggle, mdToggle);

    // Reflect simulation state from API into toggles (handles auto-start MD after load)
    function syncFromApiOnce(){
      try {
        const api = getApi(); if(!api || !api.getMetrics) return;
        // Always query the current DOM toggles so this works across re-renders/tests
        const mdToggleEl = document.getElementById('toggleMD');
        const relaxToggleEl = document.getElementById('toggleRelax');
        if (!mdToggleEl || !relaxToggleEl) return;
        const m = api.getMetrics() || {}; const kind = m.running || null;
        if (kind === 'md') {
          // MD running: set MD=On, Relax=Off
          if (mdToggleEl.getAttribute('data-on') !== 'true') setBtnState(mdToggleEl, true);
          if (relaxToggleEl.getAttribute('data-on') !== 'false') setBtnState(relaxToggleEl, false);
          try { const legacy = document.getElementById('btnMDRun'); if(legacy) legacy.textContent = 'stop'; } catch {}
          try { const status = document.getElementById('status'); if(status) status.textContent = 'MD running'; } catch {}
        } else if (kind === 'relax') {
          // Relax running: set Relax=On, MD=Off
          if (relaxToggleEl.getAttribute('data-on') !== 'true') setBtnState(relaxToggleEl, true);
          if (mdToggleEl.getAttribute('data-on') !== 'false') setBtnState(mdToggleEl, false);
          try { const legacy = document.getElementById('btnRelaxRun'); if(legacy) legacy.textContent = 'stop'; } catch {}
          try { const status = document.getElementById('status'); if(status) status.textContent = 'Relax running'; } catch {}
        } else {
          // Neither running
          if (mdToggleEl.getAttribute('data-on') !== 'false') setBtnState(mdToggleEl, false);
          if (relaxToggleEl.getAttribute('data-on') !== 'false') setBtnState(relaxToggleEl, false);
          try { const legacyMD = document.getElementById('btnMDRun'); if(legacyMD) legacyMD.textContent = 'run'; } catch {}
          try { const legacyRX = document.getElementById('btnRelaxRun'); if(legacyRX) legacyRX.textContent = 'run'; } catch {}
          try { const status = document.getElementById('status'); if(status && !/Ready|stopped/i.test(status.textContent||'')) status.textContent = 'Ready'; } catch {}
        }
      } catch {}
    }
    // Initial sync now, then keep a lightweight periodic sync so programmatic runs update UI
    try { syncFromApiOnce(); } catch {}
    try {
      if (!window.__MLIPVIEW_PANEL_SYNC_INTERVAL) {
        window.__MLIPVIEW_PANEL_SYNC_INTERVAL = setInterval(syncFromApiOnce, 250);
      }
    } catch {}

    sim.content.append(utilRow, togglesRow);

    // Sliders: keep temperature only, friction fixed at default (0.5)
    initTemperatureSlider({ hudEl: sim.content, getViewer: () => null });
  }

  // System (collapsed) — PBC moved to Simulation; keep presets and SMILES, keep backend select
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

    // Backend select only (PBC moved)
    const backendRow = document.createElement('div'); backendRow.className = 'row';
    const backendSel = document.createElement('select'); backendSel.id = 'forceProviderSel';
    backendSel.innerHTML = '<option value="fairchem">UMA (Remote)</option>';
    backendRow.append(backendSel);
    sys.content.appendChild(backendRow);
  }
  // Rendering section removed (energy plot toggle removed; forces moved to Simulation)

  // XR (collapsed)
  const xr = createSection('section-xr', 'XR', { defaultOpen: false });
  {
    const row = document.createElement('div'); row.className = 'row';
    const label = document.createElement('span'); label.className = 'form-label'; label.textContent = 'XR mode:';
    const sel = document.createElement('select'); sel.id = 'xrModeSelect';
    sel.innerHTML = '<option value="none">Off</option><option value="vr">VR</option><option value="ar">AR</option>';
    // Wire dropdown to XR entry/exit
    function getViewer(){ try { return window.viewerApi || window._viewer; } catch { return null; } }
    sel.addEventListener('change', async (ev)=>{
      const v = (ev && ev.target && ev.target.value) || sel.value;
      const api = getViewer();
      const vr = api && api.vr;
      // Fallbacks if switchXR is not present (older API)
      const doSwitch = async (mode)=>{
        try {
          if (vr && typeof vr.switchXR === 'function') {
            return await vr.switchXR(mode);
          }
          if (mode === 'none') {
            if (vr && typeof vr.exitXR === 'function') return await vr.exitXR();
            if (vr && typeof vr.exitVR === 'function') return await vr.exitVR();
            return false;
          }
          if (mode === 'vr') {
            if (vr && typeof vr.enterVR === 'function') return await vr.enterVR();
            return false;
          }
          if (mode === 'ar') {
            if (vr && typeof vr.enterAR === 'function') return await vr.enterAR();
            return false;
          }
          return false;
        } catch(e){ return false; }
      };
      if (!vr) {
        // No VR subsystem yet; reset selection and warn softly
        try { sel.value = 'none'; } catch {}
        try { const s=document.getElementById('status'); if(s) s.textContent='XR unavailable'; } catch{}
        return;
      }
      const ok = await doSwitch(v);
      // Reflect status text minimally (optional UX)
      try {
        const s = document.getElementById('status');
        if (s) {
          if (v === 'none') s.textContent = ok? 'XR exited' : 'XR exit failed';
          else s.textContent = ok? (v.toUpperCase()+' entered') : (v.toUpperCase()+' failed');
        }
      } catch {}
    });
    row.append(label, sel);
    xr.content.appendChild(row);
  }

  // Append all sections (Selection under Live Metrics)
  panel.append(live.section, selSec.section, sim.section, sys.section, xr.section);
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
        selection: selSec.content,
      simulation: sim.content,
      system: sys.content,
      xr: xr.content,
    }
  };
}
