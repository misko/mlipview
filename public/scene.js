export function createBasicScene(engine, canvas) {
  const scene = new BABYLON.Scene(engine);
  scene.clearColor = new BABYLON.Color4(0.05, 0.06, 0.08, 1.0);

  // Camera (pulled back a bit to be safe)
  const camera = new BABYLON.ArcRotateCamera(
    "cam",
    -Math.PI / 2, Math.PI / 2.35,
    16,
    BABYLON.Vector3.Zero(),
    scene
  );
  camera.wheelPrecision = 50;
  camera.panningSensibility = 800;
  camera.attachControl(canvas, true);

  // Lights
  const hemi = new BABYLON.HemisphericLight("hemi", new BABYLON.Vector3(0,1,0), scene);
  hemi.intensity = 1.3;
  const dir = new BABYLON.DirectionalLight("dir", new BABYLON.Vector3(-1,-2,-1), scene);
  dir.intensity = 0.7;

  // Ground grid for reference
  const ground = BABYLON.MeshBuilder.CreateGround("g", { width: 30, height: 30, subdivisions: 2 }, scene);
  const gmat = new BABYLON.StandardMaterial("gm", scene);
  gmat.diffuseColor = new BABYLON.Color3(0.12,0.14,0.18);
  gmat.specularColor = new BABYLON.Color3(0,0,0);
  ground.material = gmat;
  ground.position.y = -1.5;

  return scene;
}
