describe('x-vr dynamic import', () => {
  test('main-vr module exports initVRApp', async () => {
    const mod = await import('../public/vr/main-vr.js');
    expect(typeof mod.initVRApp).toBe('function');
  });
});
