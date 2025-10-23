/** @jest-environment jsdom */

describe('desktop left panel UI', () => {
  beforeEach(() => {
    document.body.innerHTML = '<div id="app" style="position:relative"></div>';
    window.__MLIP_DEV_MODE = true;
  });

  test('builds collapsible sections with defaults and controls', async () => {
    // Import builders and build panel
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');
    const { initFrictionSlider } = await import('../public/ui/frictionSlider.js');

    const res = buildDesktopPanel({ attachTo: document.getElementById('app') });
    expect(res).toBeTruthy();
    const panel = document.getElementById('controlPanel');
    expect(panel).toBeTruthy();

    // Sections exist
    const live = document.getElementById('section-live-stats');
    const sim = document.getElementById('section-simulation');
    const sys = document.getElementById('section-system');
    // Rendering section removed; XR lives as a fixed widget instead of a panel section
    expect(live && sim && sys).toBeTruthy();

    const isCollapsed = (sec) =>
      sec.querySelector('.panel-content').getAttribute('data-collapsed') === 'true';
    // Defaults: live open; simulation + system collapsed
    expect(isCollapsed(live)).toBe(false);
    expect(isCollapsed(sim)).toBe(true);
    expect(isCollapsed(sys)).toBe(true);

    // Live stats controls exist
    // Status element removed from Live Metrics
    expect(document.getElementById('instTemp')).toBeTruthy();
    expect(document.getElementById('rpsLabel')).toBeTruthy();
    expect(document.getElementById('energyPlot')).toBeTruthy();
    expect(document.getElementById('energyCanvas')).toBeTruthy();

    // Simulation controls exist
    expect(document.getElementById('btnRelax')).toBeTruthy();
    expect(document.getElementById('btnRelaxRun')).toBeTruthy();
    expect(document.getElementById('btnMD')).toBeTruthy();
    expect(document.getElementById('btnMDRun')).toBeTruthy();

    // Temperature and friction sliders (initialized during build)
    expect(document.getElementById('mdTempSlider')).toBeTruthy();
    expect(document.getElementById('mdFrictionSlider')).toBeTruthy();

    // System controls
    expect(document.getElementById('moleculeSelect')).toBeTruthy();
    expect(document.getElementById('smilesInput')).toBeTruthy();
    expect(document.getElementById('smilesGoBtn')).toBeTruthy();
    expect(document.getElementById('btnCell')).toBeTruthy();

    // Rendering controls
    expect(document.getElementById('btnToggleForces')).toBeTruthy();

    // XR controls (fixed widget)
    const xrSel = document.getElementById('xrModeSelect');
    expect(xrSel).toBeTruthy();
    expect(Array.from(xrSel.options).map((o) => o.value)).toEqual(['none', 'vr', 'ar']);
  });
});
