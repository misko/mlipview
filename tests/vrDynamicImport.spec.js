// Dynamic import test: emulate browser-style import of the VR module.
// We avoid executing Babylon specifics; just ensure the module loads and exports the function.

describe('VR dynamic import', () => {
  test('import("./public/vr/main-vr.js") resolves with initVRApp', async () => {
    const mod = await import('../public/vr/main-vr.js');
    expect(typeof mod.initVRApp).toBe('function');
  });
});
