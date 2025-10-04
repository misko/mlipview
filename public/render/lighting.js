// Viewer-consistent lighting helper.
// Creates a low ambient hemispheric fill and a camera-aligned directional "headlight".
// Directional light direction & position are synced every frame to the (active) camera
// so shading appears constant as the user pitches/yaws the view.
// Exported function returns a disposable handle.

export function attachConsistentLighting(scene, camera, opts = {}) {
  if (!scene || !camera || !scene.onBeforeRenderObservable) return null;
  const ambientIntensity = opts.ambientIntensity ?? 0.2;
  const directionalIntensity = opts.directionalIntensity ?? 0.9;
  const useActiveCamera = opts.useActiveCamera !== false; // default true

  // Remove any default hemi lights we previously created by name.
  try {
    const names = ['hemi','vrHemi'];
    for (const n of names) {
      const l = scene.getLightByName && scene.getLightByName(n);
      if (l && l.dispose) l.dispose();
    }
  } catch {}

  let ambient=null, head=null;
  try {
    ambient = new BABYLON.HemisphericLight('ambientFill', new BABYLON.Vector3(0,1,0), scene);
    ambient.intensity = ambientIntensity;
  } catch {}
  try {
    head = new BABYLON.DirectionalLight('viewerHeadlight', new BABYLON.Vector3(0,0,-1), scene);
    head.intensity = directionalIntensity;
  } catch {}

  const tmp = new BABYLON.Vector3(0,0,0);
  function updateDirection(){
    const cam = useActiveCamera && scene.activeCamera ? scene.activeCamera : camera;
    if(!cam) return;
    try {
      let fwd = null;
      if (typeof cam.getFrontPosition === 'function' && cam.position) {
        const front = cam.getFrontPosition(1);
        if (front) { tmp.copyFrom(front).subtractInPlace(cam.position).normalize(); fwd = tmp; }
      }
      if (!fwd && typeof cam.alpha === 'number' && typeof cam.beta === 'number') {
        const sinB = Math.sin(cam.beta);
        tmp.set(Math.cos(cam.alpha)*sinB, Math.cos(cam.beta), Math.sin(cam.alpha)*sinB).normalize();
        fwd = tmp;
      }
      if (fwd && head && head.direction) {
        head.direction.x = fwd.x; head.direction.y = fwd.y; head.direction.z = fwd.z;
      }
      if (head && head.position && cam.position) {
        head.position.x = cam.position.x; head.position.y = cam.position.y; head.position.z = cam.position.z;
      }
    } catch {}
  }
  const obs = scene.onBeforeRenderObservable.add(updateDirection);
  updateDirection();
  return { ambient, head, update: updateDirection, dispose(){ try { scene.onBeforeRenderObservable.remove(obs); } catch {}; try { ambient?.dispose?.(); } catch {}; try { head?.dispose?.(); } catch {}; } };
}
