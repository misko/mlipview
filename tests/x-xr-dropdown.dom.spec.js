/** @jest-environment jsdom */

describe('x-xr-dropdown', () => {
  beforeEach(() => {
    document.body.innerHTML = `
      <div id="app"></div>
      <canvas id="viewer" width="120" height="120"></canvas>
    `;

    global.window = global.window || {};
    window.viewerApi = {
      vr: {
        switchXR: jest.fn(async (mode) => mode === 'vr' || mode === 'ar' || mode === 'none'),
        enterVR: jest.fn(async () => true),
        enterAR: jest.fn(async () => true),
        exitXR: jest.fn(async () => true),
      },
    };
  });

  test('routes select changes through switchXR', async () => {
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');

    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const select = document.getElementById('xrModeSelect');
    expect(select).toBeTruthy();

    select.value = 'vr';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('vr');

    select.value = 'ar';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('ar');

    select.value = 'none';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('none');
  });

  test('falls back to enter/exit helpers when switchXR missing', async () => {
    window.viewerApi.vr = {
      enterVR: jest.fn(async () => true),
      enterAR: jest.fn(async () => true),
      exitXR: jest.fn(async () => true),
    };

    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const select = document.getElementById('xrModeSelect');

    select.value = 'vr';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.enterVR).toHaveBeenCalled();

    select.value = 'ar';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.enterAR).toHaveBeenCalled();

    select.value = 'none';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.exitXR).toHaveBeenCalled();
  });

  test('gracefully handles missing viewerApi.vr', async () => {
    delete window.viewerApi.vr;

    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const select = document.getElementById('xrModeSelect');
    select.value = 'vr';
    select.dispatchEvent(new Event('change'));
    await Promise.resolve();

    expect(select.value).toBe('none');
  });
});

