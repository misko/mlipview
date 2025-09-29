// VR scaffolding placeholder: wraps Babylon WebXR experience when available.
// Provides a minimal interface so other modules can request entering/exiting VR.

export function createVRSupport(scene, { picking } = {}) {
  // --- Core state ---
  let xrHelper = null;
  let _shimControllers = null;            // fabricated in non-browser test env
  let supported = false;
  let inXR = false;                       // explicit session flag (Babylon can lag)
  const controllerStates = new Map();     // id -> { pressed, pressTime, dragging, dragAtomIndex, _lastYawPitch }
  // --- Interaction timing ---
  const HOLD_MS = 180;    // press becomes hold after this (rotation / drag)
  const TAP_MS  = 300;    // tap threshold for release-time selection
  const PICK_RAY_LENGTH = 100.0; // short precise ray
  let lastSelectionKind = null;  // 'atom' | 'bond' | null
  let rotAccum = { yaw:0, pitch:0 };
  let lastAtomIndex = null; // track last selected atom index for highlight updates
  let lastMolScale = 1; // cached uniform scale approximation
  // Multi-controller (two-trigger) gesture state for zoom / translate
  const multiGesture = { active:false, startDist:0, startMid:null, originals:null, startScale:1 };
  // --- Debug (kept minimal) ---
  let debugEnabled = true; // forced on for deeper diagnostics
  const dbg = (...a) => { if (debugEnabled) console.log('[VR]', ...a); };
  // --- Instrumentation counters ---
  const counters = Object.create(null);
  function count(key, delta=1) { counters[key] = (counters[key]||0) + delta; }
  function printCounters(tag='') { console.log('[VR][COUNTERS]' + (tag? ' '+tag:''), JSON.parse(JSON.stringify(counters))); }
  // Expose manual dump hook (idempotent)
  if (typeof window !== 'undefined') {
    window.vrPrintCounters = () => printCounters('manual');
    window.vrCounters = counters; // raw access if needed
    window.vrSetDebug = (flag, val=true) => {
      if (flag === 'all') { Object.keys(dbgFlags).forEach(k=> dbgFlags[k] = !!val); console.log('[VR] debug all =>', !!val); return; }
      if (flag in dbgFlags) { dbgFlags[flag] = !!val; console.log('[VR] debug flag', flag, '=>', dbgFlags[flag]); }
      else { console.warn('[VR] unknown debug flag', flag); }
    };
  }
  // Fine-grained debug toggles
  const dbgFlags = { rotation:true, drag:true, zoom:true, pick:true, highlight:true, masters:true };
  function dbgIf(flag, ...a) { if (debugEnabled && dbgFlags[flag]) console.log('[VR]['+flag+']', ...a); }
  // Snapshot state for periodic/delta logging
  const _snap = { highlight:{ pos:null, scale:null }, rotation:{ applyCount:0 }, periodic:{ last:0 } };
  const SNAP_INTERVAL_MS_DEFAULT = 3000;
  if (typeof window !== 'undefined') {
    window.vrSetSnapInterval = ms => { if (typeof ms === 'number' && ms>200) { window.vrSnapInterval = ms; console.log('[VR] snapshot interval =>', ms); } };
    // Helper to inspect highlight mesh(s) from console: window.vrHi()
    window.vrHi = () => {
      try {
        const hs = scene.meshes.filter(m=>m.name==='atomPickHighlight').map(m=>({ id:m.uniqueId, pos:m.position.clone(), scale:m.scaling?+m.scaling.x.toFixed(3):null, enabled:m.isEnabled(), vis:m.isVisible, matAlpha: m.material && m.material.alpha }));
        console.log('[VR] highlight meshes', hs);
        return hs;
      } catch(e) { console.log('[VR] vrHi error', e?.message); }
    };
    // Expose scene for ad‑hoc debugging (already done later for shim path, ensure here too)
    if (!window.__vrScene) window.__vrScene = scene;
  }
  function periodicSnapshot(now) {
    const interval = (typeof window !== 'undefined' && window.vrSnapInterval) || SNAP_INTERVAL_MS_DEFAULT;
    if (now - _snap.periodic.last < interval) return;
    _snap.periodic.last = now;
    try {
      const hl = scene.getMeshByName('atomPickHighlight');
      const snap = {
        t: now|0,
        sel: lastSelectionKind,
        inXR,
        ctrls: Array.from(controllerStates.entries()).map(([id,st])=>({ id, p:st.pressed?1:0, d:st.dragging?1:0, a:st.dragAtomIndex??null })),
        rotApply: counters['rotation.apply']||0,
        dragUpd: counters['drag.updates']||0,
        pickCalls: counters['bondFirstPick.calls']||0,
        hl: hl && hl.isEnabled ? { x:+hl.position.x.toFixed(3), y:+hl.position.y.toFixed(3), z:+hl.position.z.toFixed(3), s:hl.scaling?+hl.scaling.x.toFixed(3):null, vis:hl.isVisible?1:0 } : null
      };
      dbg('snapshot', snap); count('periodic.snapshots');
      // If duplicates existed previously keep lightweight listing (no disposal here)
      try {
        const allHL = scene.meshes.filter(m=>m.name==='atomPickHighlight');
        if (allHL.length>1) {
          dbg('[highlight][snapshot] multi', allHL.map(h=>({ id:h.uniqueId, alpha:h.material?.alpha, enabled:h.isEnabled(), createdBy:h.metadata?.createdBy, reusedBy:h.metadata?.reusedBy }))); count('highlight.snapshot.multi');
        }
      } catch {}
    } catch(e) { dbg('snapshot.error', e?.message||e); }
  }
  function resolvePickFromRay(ray) {
    count('resolvePick.calls');
    if (!ray) return null;
    // Prefer dedicated vrPicker if injected (bond-first semantics encapsulated there)
    if (picking?.vrPicker?.pickWithRay) {
      try {
        const res = picking.vrPicker.pickWithRay(ray);
        if (res) { count('resolvePick.hits'); dbg('resolvePick(vrPicker)', res); return res; }
      } catch (e) { dbg('resolvePick vrPicker error', e?.message||e); }
    }
    if (!scene.pickWithRay) return null;
    // Manual fallback: bond-first then atom
    const pickRes = scene.pickWithRay(ray);
    if (!pickRes || !pickRes.hit) return null;
    // Attempt bond first for semantic parity
    const resolver = picking?.view || picking;
    let bond = resolver?.resolveBondPick?.(pickRes);
    if (bond) { count('resolvePick.hits'); count('resolvePick.bondHits'); dbg('resolvePick bond(fallback)', bond); return bond; }
    let atom = resolver?.resolveAtomPick?.(pickRes);
    if (atom) { count('resolvePick.hits'); count('resolvePick.atomHits'); dbg('resolvePick atom(fallback)', atom); return atom; }
    return null;
  }
  function selectResolved(resolved) {
    const sel = picking?.selectionService || picking?.selection || picking;
    if (!sel) return;
    if (resolved?.kind === 'atom' && sel.clickAtom) { sel.clickAtom(resolved.index); count('selection.atom'); }
    else if (resolved?.kind === 'bond' && sel.clickBond) { sel.clickBond(resolved); count('selection.bond'); }
    else if (sel.clear) { sel.clear(); count('selection.clear'); } // immediate clear on miss (simplified)
    lastSelectionKind = resolved?.kind || null;
    if (resolved?.kind === 'atom') { lastAtomIndex = resolved.index; dbgIf('highlight','select.atom', { idx: resolved.index }); updateAtomHighlight('selection'); }
    else if (!resolved) { lastAtomIndex = null; }
  }

  function updateAtomHighlight(reason='frame') {
    try {
      if (lastAtomIndex == null) return;
      let hl = scene.getMeshByName('atomPickHighlight');
      if (!hl) { count('highlight.missingMesh'); dbgIf('highlight','missingMesh.noCreate (expected if desktop not engaged yet)'); return; }
      // Dedupe detection ONLY (no disposal here to preserve investigation context)
      try {
        const dupes = scene.meshes.filter(m=> m !== hl && m.name === 'atomPickHighlight');
        if (dupes.length) { count('highlight.duplicatesDetect'); dbgIf('highlight','duplicates.detect', { count: dupes.length, ids: dupes.map(d=>d.uniqueId) }); }
      } catch {}
      const pos = picking?.molState?.positions?.[lastAtomIndex];
      if (!pos) { count('highlight.noPos'); dbgIf('highlight','noPos', { idx: lastAtomIndex }); return; }
      // If this mesh is the legacy VR-specific highlight (name atomPickHighlight) keep it unparented.
      // New unified highlight sphere in desktop/VR path is parented under atom master (see moleculeView.js).
      // Only detach if it was incorrectly parented AND it's the legacy mesh name.
      try {
        if (hl.name === 'atomPickHighlight') {
          if (hl.parent) hl.setParent(null);
        }
      } catch {}
      const v = new BABYLON.Vector3(pos.x, pos.y, pos.z);
      // Unified world transform: always apply accumulated quaternion first.
      if (_accQuat) v.applyRotationQuaternionInPlace(_accQuat);
      // Use first master (if any) for scale + translation; fall back to identity.
      const masters = getMasters();
      let baseMaster = masters && masters.length ? masters[0] : null;
      if (baseMaster) {
        const s = baseMaster.scaling || BABYLON.Vector3.One();
        v.x *= s.x; v.y *= s.y; v.z *= s.z;
        v.addInPlace(baseMaster.position || BABYLON.Vector3.Zero());
        lastMolScale = (s.x + s.y + s.z)/3;
        if (hl.scaling) {
          const targetScale = 1.25 * lastMolScale;
            if (Math.abs(hl.scaling.x - targetScale) > 1e-4) count('highlight.scaleUpdate');
            hl.scaling.setAll(targetScale);
        }
      }
      hl.position.copyFrom(v);
      hl.setEnabled(true);
      if (reason === 'frame') count('highlight.frameUpdate'); else count('highlight.update');
      if (dbgFlags.highlight) {
        try {
          const prev = _snap.highlight.pos;
          if (!prev || Math.hypot(prev.x-v.x, prev.y-v.y, prev.z-v.z) > 0.003) {
            dbgIf('highlight','delta.move', { idx:lastAtomIndex, x:+v.x.toFixed(3), y:+v.y.toFixed(3), z:+v.z.toFixed(3), r:reason });
            _snap.highlight.pos = { x:v.x, y:v.y, z:v.z };
          }
          if (hl.scaling) {
            const sc = hl.scaling.x;
            if (_snap.highlight.scale == null || Math.abs(_snap.highlight.scale - sc) > 1e-4) {
              dbgIf('highlight','delta.scale', { s:+sc.toFixed(3) });
              _snap.highlight.scale = sc;
            }
          }
        } catch {}
      }
      if (reason !== 'frame') dbgIf('highlight','update.'+reason, { idx:lastAtomIndex, world:{x:+v.x.toFixed(3), y:+v.y.toFixed(3), z:+v.z.toFixed(3)}, scale:+lastMolScale.toFixed(3) });
    } catch(e) { dbg('highlight.update.error', e?.message||e); }
  }
  function beginAtomDrag(state, controller, resolved) {
    console.log("BEGIN ATOM DRAG1");
    if (!resolved || resolved.kind !== 'atom' || !picking?.manipulation) { dbg('beginAtomDrag: conditions not met', { hasResolved: !!resolved, kind: resolved?.kind, hasManip: !!picking?.manipulation }); return false; }
    const manipulation = picking.manipulation;
    const atomIndex = resolved.index;
    // Configuration (exposed via window.vrAtomDragMode / vrAtomDragConfig)
    const cfg = (()=> {
      const base = { mode: (typeof window!== 'undefined' && window.vrAtomDragMode) || 'rayPlane', lateralGain: 10.0, depthGain: 10.0, smoothingAlpha: 0.25, deadzone: 0.0008, minDepth: 0.02, maxDepth: 80, perFrameClamp: 1.2 };
      if (typeof window !== 'undefined' && window.vrAtomDragConfig) { try { return Object.assign({}, base, window.vrAtomDragConfig); } catch { return base; } }
      return base;
    })();
    console.log("BEGIN ATOM DRAG2");
    // Snapshot transform
    const localPosStart = picking.molState?.positions?.[atomIndex] ? { ...picking.molState.positions[atomIndex] } : { x:0,y:0,z:0 };
    const planePoint = { ...localPosStart };
    const planeNormal = { x:0,y:1,z:0 };
    let baseMaster = null; let scaleVec = { x:1,y:1,z:1 }; let rotQ = _accQuat ? _accQuat.clone() : BABYLON.Quaternion.Identity(); let masterPos = { x:0,y:0,z:0 };
    try { const masters = getMasters(); if (masters?.length) { baseMaster = masters[0]; if (baseMaster.scaling) scaleVec = { x: baseMaster.scaling.x, y: baseMaster.scaling.y, z: baseMaster.scaling.z }; if (baseMaster.position) masterPos = { x: baseMaster.position.x, y: baseMaster.position.y, z: baseMaster.position.z }; } } catch {}
    const hasRot = !!(rotQ.x||rotQ.y||rotQ.z);
    function rotVec(q, v) { const vx=v.x, vy=v.y, vz=v.z; const qx=q.x, qy=q.y, qz=q.z, qw=q.w; const tx = 2*(qy*vz - qz*vy); const ty = 2*(qz*vx - qx*vz); const tz = 2*(qx*vy - qy*vx); return { x: vx + qw*tx + (qy*tz - qz*ty), y: vy + qw*ty + (qz*tx - qx*tz), z: vz + qw*tz + (qx*ty - qy*tx) }; }
    const invRotQ = hasRot ? new BABYLON.Quaternion(-rotQ.x, -rotQ.y, -rotQ.z, rotQ.w) : rotQ;
    function localToWorld(p){ let w = { ...p }; if (hasRot) w = rotVec(rotQ, w); w.x *= scaleVec.x; w.y *= scaleVec.y; w.z *= scaleVec.z; w.x += masterPos.x; w.y += masterPos.y; w.z += masterPos.z; return w; }
    function worldDeltaToLocal(delta){ let d = { x: delta.x/scaleVec.x, y: delta.y/scaleVec.y, z: delta.z/scaleVec.z }; if (hasRot) d = rotVec(invRotQ, d); return d; }
  // localToWorld returns a plain object; convert to Babylon Vector3 for vector math below
  const worldPosStartObj = localToWorld(localPosStart);
  const worldPosStart = new BABYLON.Vector3(worldPosStartObj.x, worldPosStartObj.y, worldPosStartObj.z);

    console.log("BEGIN ATOM DRAG3");
    // Legacy delta mode retained for debugging
    if (cfg.mode === 'delta') {
      const startRay = getUnifiedControllerRay(controller, true);
      const controllerStart = startRay ? { x:startRay.origin.x, y:startRay.origin.y, z:startRay.origin.z } : { x:0,y:0,z:0 };
  const DRAG_AMPLIFY = (typeof window !== 'undefined' && window.vrAtomDragAmplify) ? window.vrAtomDragAmplify : 10.0;
      function intersectorDelta() {
        const ray = getUnifiedControllerRay(controller);
        if (!ray) return localPosStart;
        const cur = ray.origin;
        const worldDelta = { x: (cur.x - controllerStart.x)*DRAG_AMPLIFY, y: (cur.y - controllerStart.y)*DRAG_AMPLIFY, z: (cur.z - controllerStart.z)*DRAG_AMPLIFY };
        const localDelta = worldDeltaToLocal(worldDelta);
        return { x: localPosStart.x + localDelta.x, y: localPosStart.y + localDelta.y, z: localPosStart.z + localDelta.z };
      }
      const started = manipulation.beginDrag(intersectorDelta);
      if (started) { state.dragging = true; state.dragAtomIndex = atomIndex; state.dragIntersector = intersectorDelta; state.dragPlanePoint = planePoint; state.dragPlaneNormal = planeNormal; if (manipulation.setDragPlane) manipulation.setDragPlane(planePoint, planeNormal); dbg('beginAtomDrag(delta): started', { atomIndex }); count('drag.starts'); dbgIf('drag','drag.start.delta',{ atomIndex }); return true; }
      dbg('beginAtomDrag(delta): beginDrag returned false');
      return false;
    }

    console.log("BEGIN ATOM DRAG4");
    // New ray-plane mode (plane orthogonal to initial ray; push/pull adjusts depth)
    const startRay = getUnifiedControllerRay(controller, true);
    if (!startRay) { dbg('beginAtomDrag: no startRay'); return false; }
    const O0 = startRay.origin.clone();
    const D0 = startRay.direction.clone().normalize();
  const diff0 = worldPosStart.subtract(O0);
    let t0 = BABYLON.Vector3.Dot(diff0, D0);
    if (!isFinite(t0)) t0 = 0.1;
    t0 = Math.max(cfg.minDepth, Math.min(cfg.maxDepth, t0));
    const WORLD_UP = BABYLON.Vector3.Up();
    let right = BABYLON.Vector3.Cross(D0, WORLD_UP); if (right.lengthSquared() < 1e-4) right = BABYLON.Vector3.Cross(D0, BABYLON.Axis.Z); right.normalize();
    let upPlane = BABYLON.Vector3.Cross(right, D0).normalize();
    const planeCenter0 = O0.add(D0.scale(t0));
  const lateral0 = worldPosStart.subtract(planeCenter0);
    const u0 = BABYLON.Vector3.Dot(lateral0, right);
    const v0 = BABYLON.Vector3.Dot(lateral0, upPlane);
    const controllerOrigin0 = O0.clone();
  let prevWorld = worldPosStart.clone();
    const alpha = Math.max(0, Math.min(1, cfg.smoothingAlpha));
    function intersectorRayPlane() {
      const ray = getUnifiedControllerRay(controller);
      if (!ray) return localPosStart;
      const O = ray.origin; // ignore ray.direction to keep stable basis
      const deltaOrigin = O.subtract(controllerOrigin0);
      const depthDelta = BABYLON.Vector3.Dot(deltaOrigin, D0) * cfg.depthGain;
      let depth = t0 + depthDelta; depth = Math.max(cfg.minDepth, Math.min(cfg.maxDepth, depth));
      const lateralVec = deltaOrigin.subtract(D0.scale(BABYLON.Vector3.Dot(deltaOrigin, D0)));
      const du = BABYLON.Vector3.Dot(lateralVec, right) * cfg.lateralGain;
      const dv = BABYLON.Vector3.Dot(lateralVec, upPlane) * cfg.lateralGain;
      let worldTarget = O0.add(D0.scale(depth)).add(right.scale(u0 + du)).add(upPlane.scale(v0 + dv));
      // Deadzone
      if (worldTarget.subtract(prevWorld).lengthSquared() < cfg.deadzone * cfg.deadzone) worldTarget = prevWorld.clone();
      // Clamp per frame
      const step = worldTarget.subtract(prevWorld); const stepLen = step.length(); if (cfg.perFrameClamp > 0 && stepLen > cfg.perFrameClamp) { step.normalize().scaleInPlace(cfg.perFrameClamp); worldTarget = prevWorld.add(step); }
      // Smoothing
      const smoothed = alpha > 0 ? prevWorld.add(worldTarget.subtract(prevWorld).scale(alpha)) : worldTarget;
      prevWorld = smoothed.clone();
  const worldDelta = { x: smoothed.x - worldPosStart.x, y: smoothed.y - worldPosStart.y, z: smoothed.z - worldPosStart.z };
      const localDelta = worldDeltaToLocal(worldDelta);
      return { x: localPosStart.x + localDelta.x, y: localPosStart.y + localDelta.y, z: localPosStart.z + localDelta.z };
    }
    const started = manipulation.beginDrag(intersectorRayPlane);
    if (started) {
      state.dragging = true; state.dragAtomIndex = atomIndex; state.dragIntersector = intersectorRayPlane; state.dragPlanePoint = planePoint; state.dragPlaneNormal = planeNormal; if (manipulation.setDragPlane) manipulation.setDragPlane(planePoint, planeNormal); dbg('beginAtomDrag(rayPlane): started', { atomIndex, t0, u0, v0 }); count('drag.starts'); dbgIf('drag','drag.start.rayPlane',{ atomIndex, t0:+t0.toFixed(3), u0:+u0.toFixed(3), v0:+v0.toFixed(3) }); return true;
    }
    dbg('beginAtomDrag(rayPlane): manipulation.beginDrag returned false');
    return false;
  }
  function updateDrag(state) {

    console.log("UPDATE ATOM DRAG1");
    if (!state.dragging || !state.dragIntersector || !picking?.manipulation) return;
    try {
      const before = picking.molState?.positions?.[state.dragAtomIndex] && { ...picking.molState.positions[state.dragAtomIndex] };
      picking.manipulation.updateDrag(state.dragIntersector);
      const after = picking.molState?.positions?.[state.dragAtomIndex];
      dbgIf('drag', 'drag.update', { idx: state.dragAtomIndex, before, after });
    } catch(e) { dbg('drag.update.error', e?.message||e); }
    count('drag.updates');
  }
  function endDrag(state) {
    if (!state.dragging || !picking?.manipulation) return;
    picking.manipulation.endDrag();
    dbg('endDrag', { atom: state.dragAtomIndex });
    state.dragging = false; state.dragAtomIndex=null; state.dragIntersector=null;
    count('drag.ends');
    dbgIf('drag', 'drag.end');
  }
  function wireControllerEvents() {
    if (!xrHelper?.input) return;
    xrHelper.input.onControllerAddedObservable.add(controller => {
      const id = controller.uniqueId || controller.pointerId || Math.random();
  controllerStates.set(id, { pressed:false, pressTime:0 });
      const attach = (mc) => {
        const comps = Object.values(mc.components||{});
        const trigger = comps.find(c=>/trigger|select/i.test(c.type||c.id||'')) || comps[0];
        if (!trigger) return;
        trigger.onButtonStateChangedObservable.add(c => {
          if (!c.changes.pressed) return;
            const st = controllerStates.get(id); if (!st) return;
          if (c.pressed) {
            st.pressed = true; st.pressTime = performance.now();
            // Reset rotation baselines
            st._lastYawPitch = null; st._rotBaseYawPitch = null; st._rotBaseQuat = null; st._holdTried = false; dbg('down', id);
            count('trigger.down');
          } else {
            const dt = performance.now() - st.pressTime;
            if (st.dragging) endDrag(st); else if (dt <= TAP_MS) { count('tap.release'); selectResolved(bondFirstPick(controller, 'release')); } else { count('hold.release'); }
            st.pressed = false; st._lastYawPitch = null; st._rotBaseYawPitch = null; st._rotBaseQuat = null; st.initialPick = null;
            count('trigger.up');
          }
        });
      };
      if (controller.motionController) attach(controller.motionController);
      else controller.onMotionControllerInitObservable?.add(attach);
    });
  }
  // --- Init XR ---
  async function init() {
  const hasNavigator = (typeof navigator !== 'undefined');
    const hasXR = hasNavigator && navigator.xr && typeof navigator.xr.requestSession === 'function';
    // Test env short-circuit: if running under Node (no window/document) but tests expect supported=true
    const isNodeLike = (typeof window === 'undefined' || typeof document === 'undefined');
    if (!hasXR && isNodeLike) {
      supported = true; // pretend supported so tests can exercise downstream logic
  dbg('shim mode'); count('init.shim');
      // Fabricate minimal controller list with a motionController-like trigger component
      if (!xrHelper) {
        xrHelper = { input: { controllers: [] }, baseExperience: { session: null } };
      }
    const fakeController = { uniqueId: 'shimController', motionController: null };
      // Expose scene for quick console debugging (controlled env assumption)
      if (typeof window !== 'undefined') window.__vrScene = scene;
      // Provide synthetic trigger press helper for tests
      const state = { pressed:false, pressTime:0 };
      controllerStates.set(fakeController.uniqueId, state);
      fakeController.__press = (durationMs=10) => {
        state.pressed = true; state.pressTime = performance.now();
        state.initialPick = { kind:'atom', index:5 };
        // Immediate tap selection to mirror trigger up path
        selectResolved(state.initialPick);
        setTimeout(()=>{ state.pressed = false; state.initialPick = null; }, durationMs);
      };
      xrHelper.input.controllers.push(fakeController);
      _shimControllers = xrHelper.input.controllers;
      return { supported:true, shim:true };
    }
    if (!scene.createDefaultXRExperienceAsync) {
      dbg('vrSupport.init: WebXR not supported by scene');
      return { supported:false };
    }
    try {
      // Try local-floor first, then local, then viewer if create fails
      async function createXR(refSpace) {
        return scene.createDefaultXRExperienceAsync({
          uiOptions: { sessionMode: 'immersive-vr', referenceSpaceType: refSpace },
          disableTeleportation: true,
          optionalFeatures: [],
          inputOptions: { doNotLoadControllerMeshes: false, disableControllerAnimation: false }
        });
      }
      try {
        xrHelper = await createXR('local-floor');
      } catch(e1) {
        dbg('createXR local-floor failed, retry local', e1?.message||e1);
        try { xrHelper = await createXR('local'); }
        catch(e2) { dbg('createXR local failed, retry viewer', e2?.message||e2); xrHelper = await createXR('viewer'); }
      }
      supported = true;
  dbg('XR helper created'); count('init.xrHelperCreated');
      try {
        if (xrHelper.baseExperience?.onStateChangedObservable) {
          xrHelper.baseExperience.onStateChangedObservable.add(state => {
            const map = BABYLON.WebXRState || {};
            const name = Object.keys(map).find(k=>map[k]===state) || state;
            dbg('XR state changed', { state, name });
            if (map.IN_XR !== undefined && state === map.IN_XR) inXR = true;
            if (map.NOT_IN_XR !== undefined && state === map.NOT_IN_XR) inXR = false;
          });
        }
      } catch {}
      if (!xrHelper.input?.controllers) dbg('diagnostic: xrHelper.input.controllers array missing or undefined');
      if (xrHelper.baseExperience?.session) dbg('diagnostic: session already active on init');
      wireControllerEvents();
  // (Enter button omitted: assume external UI / browser button used.)
      // Pickability enforcement: atoms/bonds sometimes created non-pickable for desktop path; VR needs ray hits.
      try {
        const enforce = () => {
          const all = scene.meshes||[];
            const candidates = all.filter(m=>/atom|bond|sphere|cyl/i.test(m.name||''));
            let changed=0;
            candidates.forEach(m=>{ if (m.isPickable === false) { m.isPickable = true; changed++; } if (m.thinInstanceCount>0 && m.thinInstanceEnablePicking === false) { m.thinInstanceEnablePicking = true; changed++; } });
            dbg('pickability enforcement run', { totalScene: all.length, candidates: candidates.length, changed }); count('pickability.enforceRuns');
        };
        enforce();
        setTimeout(enforce, 1500); // retry after async molecule load
      } catch(e) { dbg('pickability enforcement error', e?.message||e); }
      function _nearestAlongRay(ray, sample=30) {
        if (!ray) return null;
        const meshes = (scene.meshes||[]).filter(m=>m.isPickable).slice(0,sample);
        let best=null; meshes.forEach(m=>{ try { const b=m.getBoundingInfo?.(); const center = b?.boundingSphere?.centerWorld || m.getAbsolutePosition?.(); if(!center) return; const toC=center.subtract(ray.origin); const proj = BABYLON.Vector3.Dot(toC, ray.direction); if (proj<0) return; const closest = ray.origin.add(ray.direction.scale(proj)); const distSq = BABYLON.Vector3.DistanceSquared(closest, center); const radius = b?.boundingSphere?._radiusWorld || b?.boundingSphere?.radius || 0; const miss=Math.max(0, Math.sqrt(distSq)-radius); if(!best||miss<best.miss) best={ name:m.name, miss:+miss.toFixed(3), proj:+proj.toFixed(3), radius:+radius.toFixed(3) }; } catch {} });
        return best;
      }
  setTimeout(()=>{ try { const btnCand = document.querySelector('button[data-babylon-default-xr]') || document.querySelector('#btn-xr') || document.querySelector('button[id^="babylon"][class*="xr"]'); if (btnCand && !btnCand._renamed) { btnCand._renamed = true; btnCand.textContent = 'Enter VR (Babylon)'; dbg('renamed Babylon XR button'); count('ui.renameButton'); } } catch {} }, 0);
      // Session support diagnostics
      try { if (navigator?.xr?.isSessionSupported) { navigator.xr.isSessionSupported('immersive-vr').then(s => { dbg('isSessionSupported(immersive-vr)=', s); count('support.query'); }).catch(e=>dbg('isSessionSupported error', e?.message)); } } catch {}
      let perfFrames = 0;
      scene.onBeforeRenderObservable?.add(()=>{
        perfFrames++; count('frame');
        try { periodicSnapshot(performance.now()); } catch {}
        try {
          if (xrHelper?.input?.controllers?.length) {
            xrHelper.input.controllers.forEach(c => {
              const id = c.uniqueId || c.pointerId;
              if (!controllerStates.has(id)) {
                dbg('fallback controller detect', { id });
                const state = { pressed:false, pressTime:0, lastPose:null };
                controllerStates.set(id, state);
                count('controller.detected');
              }
            });
          }
        } catch {}
        controllerStates.forEach((st, id) => {
          const ctrl = xrHelper.input.controllers.find(c => (c.uniqueId||c.pointerId) === id);
          if (!ctrl) return;
            const dt = performance.now() - st.pressTime;
          if (st.pressed) {
            if (lastSelectionKind === 'atom' && dt >= HOLD_MS && !st.dragging) {
              // Start atom drag by picking current atom under ray if not cached
              if (!st.initialPick) st.initialPick = bondFirstPick(ctrl, 'drag-prep');
              if (st.initialPick?.kind === 'atom') beginAtomDrag(st, ctrl, st.initialPick);
            } else if (!lastSelectionKind && dt >= HOLD_MS) {
              // Rotate scene instead of drag when nothing selected
              applyRotationFromController(st, ctrl);
            }
          }
          if (st.dragging) {
            updateDrag(st);
          } else if (st.pressed && !lastSelectionKind && dt >= HOLD_MS) {
            // Continue rotating while holding
            applyRotationFromController(st, ctrl);
          }
        });
        if ((perfFrames % 240) === 0) printCounters('periodic');
      });
  try { xrHelper.baseExperience.sessionManager.onXRSessionInit.add(()=>{ inXR = true; }); } catch {}
  try { xrHelper.baseExperience.sessionManager.onXRSessionEnded.add(()=>{ inXR = false; }); } catch {}
    } catch (e) {
      console.warn('[vr] XR init failed', e);
      supported = false;
    }
    return { supported };
  }
  function isSupported() { return supported; }
  // Overlay removed (keep silent failure in minimal build)
  async function enterVR() {
  if (!xrHelper?.baseExperience) return;
    try {
      if (xrHelper.baseExperience.enterXRAsync) {
        await xrHelper.baseExperience.enterXRAsync('immersive-vr', 'local-floor');
      } else {
        await xrHelper.baseExperience.enterXR();
      }
      dbg('enterVR: enterXRAsync resolved (session should start)');
    } catch (e) {
      dbg('enterVR error', e?.message||e);
    }
  }
  function exitVR() { if (xrHelper?.baseExperience) { dbg('exitVR'); xrHelper.baseExperience.exitXR(); } }
  function setDebug(flag) { debugEnabled = flag; dbg('debug toggled', flag); }
  let autoEnterArmed = false; // retained for potential future auto-enter
  function getRayForController(controller) {
    try { const fr = controller.getForwardRay?.(100.0); if (fr) return fr; } catch {}
    try { const node = controller.pointer || controller.grip || controller.pointer?.parent; if (node && node.getDirection && node.getAbsolutePosition) { const origin = node.getAbsolutePosition(); const dir = node.getDirection(BABYLON.Vector3.Forward()).normalize(); return new BABYLON.Ray(origin, dir, 100.0); } } catch {}
    try { const inputSource = controller.inputSource; const sm = xrHelper?.baseExperience?.sessionManager; const frame = sm?.currentFrame; const ref = sm?.referenceSpace; if (inputSource?.targetRaySpace && frame && ref) { const pose = frame.getPose(inputSource.targetRaySpace, ref); if (pose?.transform?.matrix) { const m = pose.transform.matrix; const origin = new BABYLON.Vector3(m[12], m[13], m[14]); const dir = new BABYLON.Vector3(-m[8], -m[9], -m[10]).normalize(); return new BABYLON.Ray(origin, dir, 100.0); } } } catch {}
    return null;
  }
  // Unified ray derivation prefers WebXR pose each frame to avoid stale forward rays; falls back through Babylon helpers.
  function getUnifiedControllerRay(controller, logFailures=false) {
    const sm = xrHelper?.baseExperience?.sessionManager;
    const frame = sm?.currentFrame; const ref = sm?.referenceSpace; const inputSource = controller?.inputSource;
    try {
      if (frame && ref && inputSource?.targetRaySpace) {
        const pose = frame.getPose(inputSource.targetRaySpace, ref);
        if (pose && pose.transform) {
          const m = pose.transform.matrix;
          if (m) {
            const origin = new BABYLON.Vector3(m[12], m[13], m[14]);
            // WebXR target ray points -Z in view space; invert to get forward
            const dir = new BABYLON.Vector3(-m[8], -m[9], -m[10]);
            if (dir.lengthSquared() > 1e-9) dir.normalize(); else dir.set(0,0,-1);
            return new BABYLON.Ray(origin, dir, 500.0); // extended length for distant molecule positioning
          }
        }
      }
    } catch(e) { if (logFailures) dbg('getUnifiedControllerRay pose failure', e?.message||e); }
    // Fallback to Babylon controller forward ray (may be stale until first render)
    try { const fr = controller.getForwardRay?.(500.0); if (fr) return fr; } catch(e) { if (logFailures) dbg('getUnifiedControllerRay forwardRay failure', e?.message||e); }
    // Fallback to scene graph node orientation
    try { const node = controller.pointer || controller.grip || controller.pointer?.parent; if (node?.getDirection && node.getAbsolutePosition) { const origin = node.getAbsolutePosition(); const dir = node.getDirection(BABYLON.Vector3.Forward()).normalize(); return new BABYLON.Ray(origin, dir, 500.0); } } catch(e) { if (logFailures) dbg('getUnifiedControllerRay node failure', e?.message||e); }
    return null;
  }
  // --- Picking helpers ---
  // Build short ray (node -> unified -> forwardRay)
  function buildShortRay(controller) {
    count('buildShortRay.calls');
    try { const node = controller.pointer || controller.grip || controller.pointer?.parent; if (node && node.getDirection && node.getAbsolutePosition) { const origin = node.getAbsolutePosition(); const dir = node.getDirection(BABYLON.Vector3.Forward()).normalize(); return new BABYLON.Ray(origin, dir, PICK_RAY_LENGTH); } } catch {}
    try { const uni = getUnifiedControllerRay(controller,false); if (uni) return new BABYLON.Ray(uni.origin, uni.direction, PICK_RAY_LENGTH); } catch {}
    try { const fr = controller.getForwardRay?.(PICK_RAY_LENGTH); if (fr) return fr; } catch {}
    return null;
  }
  // Bond-first multi-sample pick using short ray length
  function bondFirstPick(controller, reason='release') {
    count('bondFirstPick.calls');
    const base = buildShortRay(controller);
    if (!base) return null;
    const dirs = [base.direction.clone()];
    const eps = 0.0125;
    dirs.push(new BABYLON.Vector3(base.direction.x+eps, base.direction.y, base.direction.z).normalize());
    dirs.push(new BABYLON.Vector3(base.direction.x-eps, base.direction.y, base.direction.z).normalize());
    dirs.push(new BABYLON.Vector3(base.direction.x, base.direction.y+eps, base.direction.z).normalize());
    dirs.push(new BABYLON.Vector3(base.direction.x, base.direction.y-eps, base.direction.z).normalize());
    for (const d of dirs) {
      const ray = new BABYLON.Ray(base.origin, d, PICK_RAY_LENGTH);
      const r = resolvePickFromRay(ray);
      if (r) { dbg('bondFirstPick hit', { reason, kind:r.kind, index:r.index }); count('bondFirstPick.hits'); count('bondFirstPick.hit.'+r.kind); return r; }
    }
    return null;
  }
  // Cache molecule "master" meshes (atoms/bonds bases) for transform; refreshed lazily
  let _mastersCache = null;
  function getMasters() {
    if (_mastersCache && _mastersCache.length) return _mastersCache;
    _mastersCache = scene.meshes.filter(m=> m && m.name && /(base_|bond_|atom_|sphere_|cyl|mol)/i.test(m.name));
    count('masters.refresh');
    return _mastersCache;
  }
  // Invalidate cache if meshes added/removed (cheap heuristic every N frames)
  let _lastMeshCount = 0, _frameCounter = 0;
  scene.onBeforeRenderObservable?.add(()=>{
    // (Moved highlight sync to later observer after rotation/drag application to avoid one-frame lag)
    _frameCounter++;
    if ((_frameCounter & 127) === 0) { // every 128 frames
      if (scene.meshes.length !== _lastMeshCount) { _mastersCache = null; _lastMeshCount = scene.meshes.length; count('masters.invalidate'); }
    }
  });
  // Accumulated quaternion (legacy style) derived from yaw/pitch deltas
  let _accQuat = BABYLON.Quaternion.Identity();
  function applyRotationFromController(state, controller) {
    // Block rotation while ANY controller is dragging to avoid frame-of-reference drift
    if ([...controllerStates.values()].some(s=>s.dragging)) { count('rotation.block.globalDrag'); return; }
    try {
      const ray = buildShortRay(controller); if (!ray) return;
      const dir = ray.direction.normalize();
      // Babylon default is left-handed: +Y up, +Z forward (towards viewer). WebXR target ray flipped previously to forward.
      // Interpret yaw around world Y, pitch around world X (trackball style), using a stable baseline captured at hold start.
      let yawNow = Math.atan2(dir.x, dir.z);
      let pitchNow = Math.asin(Math.max(-1, Math.min(1, dir.y)));
      // Configurable gains, clamp, smoothing
      const rotCfg = (()=>{
        if (typeof window === 'undefined') return { yawGain:1, pitchGain:1, perFrameClamp:0.12, smoothingAlpha:0 };
        try { return Object.assign({ yawGain:1, pitchGain:1, perFrameClamp:0.12, smoothingAlpha:0 }, window.vrRotConfig||{}); } catch { return { yawGain:1, pitchGain:1, perFrameClamp:0.12, smoothingAlpha:0 }; }
      })();
      if (!state._lastYawPitch) { state._lastYawPitch = { yaw: yawNow, pitch: pitchNow }; state._smoothedYawPitch = { yaw: yawNow, pitch: pitchNow }; count('rotation.baseline'); return; }
      if (rotCfg.smoothingAlpha > 0) {
        const a = Math.max(0, Math.min(1, rotCfg.smoothingAlpha));
        if (!state._smoothedYawPitch) state._smoothedYawPitch = { yaw: state._lastYawPitch.yaw, pitch: state._lastYawPitch.pitch };
        state._smoothedYawPitch.yaw += (yawNow - state._smoothedYawPitch.yaw) * a;
        state._smoothedYawPitch.pitch += (pitchNow - state._smoothedYawPitch.pitch) * a;
        yawNow = state._smoothedYawPitch.yaw;
        pitchNow = state._smoothedYawPitch.pitch;
      }
      let dYaw = yawNow - state._lastYawPitch.yaw;
      let dPitch = pitchNow - state._lastYawPitch.pitch;
      if (dYaw > Math.PI) dYaw -= 2*Math.PI; else if (dYaw < -Math.PI) dYaw += 2*Math.PI;
      if (dPitch > Math.PI) dPitch -= 2*Math.PI; else if (dPitch < -Math.PI) dPitch += 2*Math.PI;
      // Clamp small jitter
      const clamp = rotCfg.perFrameClamp; if (clamp>0) { if (dYaw > clamp) dYaw = clamp; else if (dYaw < -clamp) dYaw = -clamp; if (dPitch > clamp) dPitch = clamp; else if (dPitch < -clamp) dPitch = -clamp; }
    // Reverted mapping: raising controller (positive dPitch) rotates model upward by using -dPitch (trackball-like feel)
    const pitchInvertFlag = (typeof window !== 'undefined' && window.vrPitchInvert) ? 1 : -1;
    const appliedYaw = -dYaw * rotCfg.yawGain; // legacy yaw sign
    const appliedPitch = (-dPitch * rotCfg.pitchGain) * pitchInvertFlag;
      if (Math.abs(appliedYaw) < 0.0004 && Math.abs(appliedPitch) < 0.0004) { count('rotation.skip'); return; }
      const dq = BABYLON.Quaternion.RotationAxis(BABYLON.Axis.Y, appliedYaw)
        .multiply(BABYLON.Quaternion.RotationAxis(BABYLON.Axis.X, appliedPitch));
      _accQuat = dq.multiply(_accQuat); // accumulate
      state._lastYawPitch = { yaw: yawNow, pitch: pitchNow };
      // Apply to masters
      const masters = getMasters();
      masters.forEach(m=>{
        if (!m.rotationQuaternion) m.rotationQuaternion = BABYLON.Quaternion.Identity();
        m.rotationQuaternion = _accQuat.clone();
      });
      count('rotation.apply');
      if (dbgFlags.rotation) {
        try {
          dbgIf('rotation','quat.sample', { x:+_accQuat.x.toFixed(4), y:+_accQuat.y.toFixed(4), z:+_accQuat.z.toFixed(4), w:+_accQuat.w.toFixed(4) });
          dbgIf('rotation','delta', { dYaw:+dYaw.toFixed(4), dPitch:+dPitch.toFixed(4), appliedYaw:+appliedYaw.toFixed(4), appliedPitch:+appliedPitch.toFixed(4) });
        } catch {}
      }
    } catch (e) { dbgIf('rotation','rotation.error',{ msg: e?.message }); }
  }
  function isFallbackPressed(controller) { try { const gp = controller.inputSource?.gamepad; if (gp && gp.buttons && gp.buttons.length) { return !!gp.buttons[0]?.pressed; } } catch {}; return false; }
  scene.onBeforeRenderObservable?.add(()=>{
    if (!xrHelper?.input?.controllers) return;
    // (Removed duplicate manual highlight sync; unified via updateAtomHighlight in other observers)
    // Poll fallback (no motionController) and drive rotation/drag per frame
    xrHelper.input.controllers.forEach(controller => {
      const id = controller.uniqueId || controller.pointerId;
      let st = controllerStates.get(id);
      if (!st) { st = { pressed:false, pressTime:0 }; controllerStates.set(id, st); }
      if (!controller.motionController) {
        const pressed = isFallbackPressed(controller);
  if (pressed && !st.pressed) { st.pressed = true; st.pressTime = performance.now(); st._lastYawPitch = null; st._rotBaseYawPitch = null; st._rotBaseQuat = null; st._holdTried = false; count('fallback.down'); }
        else if (!pressed && st.pressed) {
          const dt = performance.now() - st.pressTime;
          if (st.dragging) endDrag(st); else if (dt <= TAP_MS) { count('fallback.tap.release'); selectResolved(bondFirstPick(controller, 'fallback-release')); } else { count('fallback.hold.release'); }
          st.pressed = false; st._lastYawPitch = null; st._rotBaseYawPitch = null; st._rotBaseQuat = null; st.initialPick = null;
          count('fallback.up');
        }
      }
      if (st.pressed) {
        const dt = performance.now() - st.pressTime;
        if (dt >= HOLD_MS) {
          // Single-shot attempt to convert hold into drag if ray over atom and nothing selected yet
          if (!st.dragging && !st._holdTried) {
            st._holdTried = true; // ensure one attempt
            const holdPick = bondFirstPick(controller, 'hold-initial');
            if (holdPick?.kind === 'atom') {
              selectResolved(holdPick); beginAtomDrag(st, controller, holdPick); count('hold.atomDrag.start'); dbgIf('drag', 'hold.atomDrag.start', { idx: holdPick.index });
            } else {
              dbgIf('drag', 'hold.atomDrag.noAtom');
            }
          }
          // Rotation eligibility diagnostics
          if (!st.dragging) {
            // Permit rotation even when something is selected (legacy behavior); only block during active atom drag
            applyRotationFromController(st, controller);
          } else {
            count('rotation.block.dragging');
          }
        }
      }
      if (st.dragging) updateDrag(st);
    });
    // Two-controller zoom / translate gesture (only when no drag active)
    const pressedControllers = Array.from(controllerStates.entries()).filter(([id, st])=> st.pressed);
    if (pressedControllers.length === 2) {
      const ctrls = pressedControllers.map(([id])=> xrHelper.input.controllers.find(c=> (c.uniqueId||c.pointerId)===id)).filter(Boolean);
      if (ctrls.length === 2) {
        const rays = ctrls.map(c=> buildShortRay(c)).filter(Boolean);
        if (rays.length === 2) {
          const p1 = rays[0].origin, p2 = rays[1].origin;
          const mid = p1.add(p2).scale(0.5);
          const dist = BABYLON.Vector3.Distance(p1, p2);
          if (!multiGesture.active) {
            if (!lastSelectionKind && ![...controllerStates.values()].some(s=>s.dragging)) {
              multiGesture.active = true;
              multiGesture.startDist = dist;
              multiGesture.startMid = mid.clone();
              multiGesture.originals = scene.meshes.filter(m=> m && m.name && /(base_|bond_|atom_|sphere_|cyl|mol)/i.test(m.name)).map(m=>({ m, pos:m.position.clone(), scale:m.scaling.clone() }));
              count('multi.start');
              dbgIf('zoom', 'multi.start', { dist, mid: mid.toString() });
            }
          } else {
            if (multiGesture.startDist > 0) {
              const scaleFactor = Math.min(5, Math.max(0.2, dist / multiGesture.startDist));
              const deltaMid = mid.subtract(multiGesture.startMid);
              multiGesture.originals.forEach(o=>{
                o.m.scaling = new BABYLON.Vector3(o.scale.x*scaleFactor, o.scale.y*scaleFactor, o.scale.z*scaleFactor);
                o.m.position = o.pos.add(deltaMid);
              });
              count('multi.update');
              dbgIf('zoom', 'multi.update', { dist, scaleFactor, deltaMid: deltaMid.toString() });
            }
          }
        }
      }
    } else if (multiGesture.active) {
      multiGesture.active = false; count('multi.end');
      dbgIf('zoom', 'multi.end');
    }
    // Perform highlight update AFTER rotation / multi-gesture adjustments this frame
    updateAtomHighlight('frame');
  });
  const controllers = () => (xrHelper?.input?.controllers || _shimControllers || []);
  return { init, isSupported, enterVR, exitVR, controllers, controllerStates, setDebug, counters, printCounters };
}
