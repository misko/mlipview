import { attachConsistentLighting } from './lighting.js';

export async function createScene(canvas) {
  const engine = new BABYLON.Engine(canvas, true, { preserveDrawingBuffer:true, stencil:true });
  const scene = new BABYLON.Scene(engine);
  // Ensure pointer observable exists (some Babylon builds include it, but guard for safety)
  if (!scene.onPointerObservable) {
    scene.onPointerObservable = { _l:[], add(fn){ this._l.push(fn); }, notify(ev){ this._l.forEach(f=>f(ev)); }, notifyObservers(ev){ this._l.forEach(f=>f(ev)); } };
  } else {
    // Normalize API: ensure both notify & notifyObservers exist for our usage
    const o = scene.onPointerObservable;
    if (o.notify && !o.notifyObservers) o.notifyObservers = o.notify.bind(o);
    if (o.notifyObservers && !o.notify) o.notify = o.notifyObservers.bind(o);
  }
  // Mirror engine pointer coordinates onto scene for picking service convenience
  if (engine) {
    scene.pointerX = 0; scene.pointerY = 0;
    canvas.addEventListener('pointermove', e => { const rect = canvas.getBoundingClientRect(); scene.pointerX = e.clientX - rect.left; scene.pointerY = e.clientY - rect.top; });
    // Only synthesize pointerdown if Babylon core doesn't already populate it (heuristic: check for internal observable class shape)
    // If Babylon's native onPointerObservable exists, it will call observers itself. We mark ours if we create it manually.
    const createdManually = !scene.onPointerObservable._isNative;
    if (createdManually) {
      canvas.addEventListener('pointerdown', e => { scene.onPointerObservable.notify({ type: BABYLON.PointerEventTypes ? BABYLON.PointerEventTypes.POINTERDOWN : 1, event:e }); });
    }
  }
  // Set background to white (desktop + baseline for XR background override)
  scene.clearColor = new BABYLON.Color4(1,1,1,1);
  // Move camera closer so objects appear ~50% larger (radius scaled by 1/1.5)
  const camera = new BABYLON.ArcRotateCamera('cam', Math.PI/4, Math.PI/3, 25/1.5, new BABYLON.Vector3(0,0,0), scene);
  camera.attachControl(canvas, true);
  attachConsistentLighting(scene, camera, { ambientIntensity:0.2, directionalIntensity:0.9 });
  return { engine, scene, camera };
}
