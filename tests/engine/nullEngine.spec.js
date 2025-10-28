const { createNullScene, MeshBuilder } = require('./nullEngineHarness');

describe('NullEngine harness basic', () => {
  test('can create a sphere mesh placeholder', () => {
    const { scene } = createNullScene();
    const sphere = MeshBuilder.CreateSphere('s', { diameter: 1 });
    // In fallback mock, meshes might not auto-register; tolerate missing push logic
    if (scene.meshes && !scene.meshes.includes(sphere)) scene.meshes.push(sphere);
    expect(sphere.name).toBe('s');
  });
});
