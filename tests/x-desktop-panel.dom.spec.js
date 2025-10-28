/** @jest-environment jsdom */

beforeEach(() => {
  document.body.innerHTML = '<div id="mount"></div>';
});

describe('x-desktop-panel', () => {
  test('renders expected sections and controls', async () => {
    window.__MLIP_DEV_MODE = true;
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');

    const api = buildDesktopPanel({ attachTo: document.getElementById('mount') });
    expect(api).toBeTruthy();

    const panel = document.getElementById('controlPanel');
    expect(panel).toBeTruthy();

    const live = document.getElementById('section-live-stats');
    const sim = document.getElementById('section-simulation');
    const sys = document.getElementById('section-system');
    expect(live && sim && sys).toBeTruthy();

    const isCollapsed = (id) =>
      document
        .getElementById(id)
        .querySelector('.panel-content')
        .getAttribute('data-collapsed') === 'true';
    expect(isCollapsed('section-live-stats')).toBe(false);
    expect(isCollapsed('section-simulation')).toBe(true);
    expect(isCollapsed('section-system')).toBe(true);

    expect(document.getElementById('btnRelax')).toBeTruthy();
    expect(document.getElementById('btnMD')).toBeTruthy();
    expect(document.getElementById('mdTempSlider')).toBeTruthy();
    expect(document.getElementById('mdFrictionSlider')).toBeTruthy();

    expect(document.getElementById('moleculeSelect')).toBeTruthy();
    expect(document.getElementById('btnCell')).toBeTruthy();
    expect(document.getElementById('btnToggleForces')).toBeTruthy();

    const xrSelect = document.getElementById('xrModeSelect');
    expect(xrSelect).toBeTruthy();
    expect(Array.from(xrSelect.options).map((o) => o.value)).toEqual(['none', 'vr', 'ar']);
  });
});
