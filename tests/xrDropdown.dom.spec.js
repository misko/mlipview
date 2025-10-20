/** @jest-environment jsdom */

// Verify that changing the XR mode dropdown invokes the VR API hooks
// and updates status text accordingly.

describe('XR dropdown -> VR API wiring (jsdom)', () => {
  beforeEach(() => {
    // Minimal DOM shell expected by desktopPanel
    document.body.innerHTML = `
      <div id="app"></div>
      <canvas id="viewer" width="100" height="100"></canvas>
    `;
    // Ensure a window.viewerApi.vr mock is present
    global.window = global.window || {};
    window.viewerApi = {
      vr: {
        switchXR: jest.fn(async (mode) => {
          // simulate success for vr/ar/none
          return mode === 'vr' || mode === 'ar' || mode === 'none';
        }),
        enterVR: jest.fn(async () => true),
        enterAR: jest.fn(async () => true),
        exitXR: jest.fn(async () => true),
      },
    };
  });

  test('changing select calls switchXR with vr/ar/none', async () => {
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const sel = document.getElementById('xrModeSelect');
    expect(sel).toBeTruthy();

    // Select VR
    sel.value = 'vr';
    sel.dispatchEvent(new Event('change'));
    // Allow microtask
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('vr');

    // Select AR
    sel.value = 'ar';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('ar');

    // Select Off
    sel.value = 'none';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.switchXR).toHaveBeenCalledWith('none');
  });

  test('fallback to enterVR/enterAR/exitXR when switchXR missing', async () => {
    // Replace vr object without switchXR
    window.viewerApi.vr = {
      enterVR: jest.fn(async () => true),
      enterAR: jest.fn(async () => true),
      exitXR: jest.fn(async () => true),
    };
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const sel = document.getElementById('xrModeSelect');
    // VR
    sel.value = 'vr';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.enterVR).toHaveBeenCalled();
    // AR
    sel.value = 'ar';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.enterAR).toHaveBeenCalled();
    // None
    sel.value = 'none';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    expect(window.viewerApi.vr.exitXR).toHaveBeenCalled();
  });

  test('graceful when viewerApi.vr missing', async () => {
    delete window.viewerApi.vr;
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const sel = document.getElementById('xrModeSelect');
    sel.value = 'vr';
    sel.dispatchEvent(new Event('change'));
    await Promise.resolve();
    // No throws; selection auto-resets to none in handler
    expect(sel.value).toBe('none');
  });
});
