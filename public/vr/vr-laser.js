// VR Laser management: creates and updates controller laser beams with dynamic thickness.
// Extracted from setup.js to keep core VR setup lean.
// Public factory: createVRLaserManager({ scene, xrHelper, controllerState, getControllerNode, pickAtomOrBondWithRay })
// Returned API: { ensureForController(c), updateFrame(inputs), getInfo() }

export function createVRLaserManager({ scene, xrHelper, controllerState, getControllerNode, pickAtomOrBondWithRay }){
  const cfg = {
    baseDiam: (typeof window!=='undefined' && window.vrLaserBaseDiam) || 0.010,
    hoverDiam: (typeof window!=='undefined' && window.vrLaserHoverDiam) || 0.020,
  // Nominal maximum visual length. Treated as effectively "infinite" unless hit found sooner.
  maxLength: (typeof window!=='undefined' && window.vrLaserMaxLength) || 60.0,
    lerp: (typeof window!=='undefined' && window.vrLaserLerp) || 0.25,
    baseColor: { r:0.35,g:0.8,b:1.0 },
    hoverColor: { r:0.9,g:0.95,b:1.0 }
  };
  const tmpForward = { x:0,y:0,z:1 };
  function ensureForController(c){
    let st = controllerState.get(c);
    if(!st){ st={ pressed:false,lastYaw:0,lastPitch:0,downTime:0,tapHandled:false }; controllerState.set(c,st); }
    if(typeof st.laserDiameterCurrent !== 'number') st.laserDiameterCurrent = cfg.baseDiam;
    if(!st.laserMesh && typeof BABYLON!=='undefined'){
      try {
        const node = getControllerNode(c);
        if(node){
          const mesh = BABYLON.MeshBuilder.CreateCylinder('xrLaser_'+(c.uniqueId||'u'), { height:cfg.maxLength, diameter:cfg.baseDiam, tessellation:12, subdivisions:1 }, scene);
          mesh.isPickable=false;
          mesh.alwaysSelectAsActiveMesh=false;
          mesh.parent=node;
          mesh.position = new BABYLON.Vector3(0,0,cfg.maxLength/2); // along local +Z (Babylon forward)
          // Orient cylinder so its axis aligns with +Z: Babylon cylinders point up (Y) by default; rotate -90° about X.
          mesh.rotationQuaternion = BABYLON.Quaternion.RotationAxis(BABYLON.Axis.X, -Math.PI/2);
          const mat = new BABYLON.StandardMaterial('xrLaserMat_'+(c.uniqueId||'u'), scene);
          mat.disableLighting=true;
          mat.emissiveColor = new BABYLON.Color3(cfg.baseColor.r,cfg.baseColor.g,cfg.baseColor.b);
          mat.alpha=0.95;
          mesh.material=mat;
          st.laserMesh=mesh;
          st.hoverKind=null;
        }
      } catch {/* ignore */}
    }
  }
  function updateController(c){
    const st = controllerState.get(c); if(!st) return;
    // Build forward pick ray for hover detection
    let hoverKind=null;
    let hitDistance=null;
    if(typeof BABYLON!=='undefined'){
      try {
        const node = getControllerNode(c);
        if(node?.getDirection && node.getAbsolutePosition){
          const origin=node.getAbsolutePosition();
            const dir=node.getDirection(BABYLON.Vector3.Forward()).normalize();
          const ray=new BABYLON.Ray(origin,dir,cfg.maxLength);
          const res = pickAtomOrBondWithRay(ray);
          if(res){
            if(res.kind==='atom') hoverKind='atom'; else if(res.kind==='bond') hoverKind='bond';
            if(typeof res.dist==='number') hitDistance=res.dist;
          }
        }
      } catch{}
    }
    st.hoverKind=hoverKind;
    const desired = hoverKind? cfg.hoverDiam : cfg.baseDiam;
    st.laserDiameterCurrent += (desired - st.laserDiameterCurrent)*cfg.lerp;
    // Apply scale & color
    const mesh=st.laserMesh;
    if(mesh){
      const scaleFactor = st.laserDiameterCurrent / cfg.baseDiam;
      mesh.scaling.x = scaleFactor; // radial axes after rotation (X/Z after rotation represent diameter axes)
      mesh.scaling.z = scaleFactor;
      // Adjust visible length: Babylon cylinder originally along Y (before rotation). After rotation -90° about X, its local Z length corresponds to original Y scale.
      // We keep base mesh at maxLength and shorten by scaling Y plus repositioning so base stays at controller.
      let lengthScale=1.0;
      if(typeof hitDistance==='number' && hitDistance>0){
        const clamped = Math.min(hitDistance, cfg.maxLength);
        lengthScale = clamped / cfg.maxLength;
      }
      mesh.scaling.y = lengthScale;
      // Reposition so tip ends at hit (mesh centered at half-length along +Z after rotation+position). Original position.z = maxLength/2. After scaling Y, effective half-length becomes (maxLength*lengthScale)/2.
      if(mesh.position){
        mesh.position.z = (cfg.maxLength * lengthScale)/2;
      }
      const mat = mesh.material;
      if(mat?.emissiveColor){
        if(hoverKind){
          mat.emissiveColor.r=cfg.hoverColor.r; mat.emissiveColor.g=cfg.hoverColor.g; mat.emissiveColor.b=cfg.hoverColor.b; mat.alpha=1.0;
        } else {
          mat.emissiveColor.r=cfg.baseColor.r; mat.emissiveColor.g=cfg.baseColor.g; mat.emissiveColor.b=cfg.baseColor.b; mat.alpha=0.95;
        }
      }
    }
  }
  function updateFrame(inputs){
    for(const c of inputs){ ensureForController(c); updateController(c); }
  }
  function getInfo(){
    const out=[];
    try {
      const inputs = (xrHelper && xrHelper.input && xrHelper.input.controllers) ? xrHelper.input.controllers : [];
      for (const c of inputs) {
        const st = controllerState.get(c);
        if (!st) { out.push({ id: c?.uniqueId || 'unknown', hoverKind: null, diameter: undefined }); continue; }
        out.push({ id: c?.uniqueId || 'unknown', hoverKind: st.hoverKind || null, diameter: st.laserDiameterCurrent || null });
      }
    } catch {}
    return out;
  }
  return { ensureForController, updateFrame, getInfo, config:cfg };
}
