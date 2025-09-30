// Unified manipulation service for atom dragging & bond rotation.
// Works headless; caller supplies ray -> plane intersection utilities from rendering layer.
// API:
//   beginDrag(pointerRay) if atom selected
//   updateDrag(pointerRay)
//   endDrag()
//   rotateBond(deltaAngleRadians) (uses current bond orientation side from selection)
// Emits position change via molState.markPositionsChanged(); does not recompute bonds.

import { orientationToSide } from '../selection-model.js';

// Optional deps: { bondService } to trigger bond recomputation after an atom drag.
export function createManipulationService(molState, { bondService } = {}) {
  let dragState = null; // { atomIndex, startPos:{x,y,z}, grabOffset:{x,y,z}, planeNormal:{x,y,z}, planePoint:{x,y,z} }
  // Debug instrumentation (re-added to surface desktop drag logs after migration).
  const DRAG_LOG = false; // silenced (set to true locally if detailed drag diagnostics needed)

  function getAtomPosition(i) { return molState.positions[i]; }
  function setAtomPosition(i, p) { molState.positions[i].x=p.x; molState.positions[i].y=p.y; molState.positions[i].z=p.z; }

  function beginDrag(intersector) {
    if (dragState) return false;
    if (!molState.selection || molState.selection.kind !== 'atom') return false;
    const atomIndex = molState.selection.data.index;
    const pos = getAtomPosition(atomIndex);
  // Use Y-up horizontal plane so dragging feels natural relative to ground.
  const planeNormal = { x:0, y:1, z:0 }; // default; caller can override later via setDragPlane
    const planePoint = { ...pos };
    const hit = intersector(planePoint, planeNormal);
    const grabOffset = hit ? { x:pos.x-hit.x, y:pos.y-hit.y, z:pos.z-hit.z } : { x:0,y:0,z:0 };
    dragState = { atomIndex, startPos: { ...pos }, grabOffset, planeNormal, planePoint };
  if (DRAG_LOG) console.log('[drag][start]', { atomIndex, startPos: { ...pos }, planePoint: { ...planePoint }, planeNormal: { ...planeNormal } }); // retained behind flag
    return true;
  }
  function setDragPlane(point, normal) {
    if (!dragState) return;
    dragState.planePoint = { ...point };
    dragState.planeNormal = { ...normal };
  }
  function updateDrag(intersector) {
    if (!dragState) return false;
    const hit = intersector(dragState.planePoint, dragState.planeNormal);
    if (!hit) return false;
    const prev = { ...molState.positions[dragState.atomIndex] };
    const newPos = { x: hit.x + dragState.grabOffset.x, y: hit.y + dragState.grabOffset.y, z: hit.z + dragState.grabOffset.z };
    setAtomPosition(dragState.atomIndex, newPos);
    molState.markPositionsChanged();
    if (DRAG_LOG) {
      const dx = +(newPos.x - prev.x).toFixed(3);
      const dy = +(newPos.y - prev.y).toFixed(3);
      const dz = +(newPos.z - prev.z).toFixed(3);
      console.log('[drag][move]', { atomIndex: dragState.atomIndex, from: prev, to: { ...newPos }, d: { dx, dy, dz } });
    } // suppressed by default
    return true;
  }
  function endDrag() {
    if (!dragState) return;
    // Capture whether the atom actually moved (compare to startPos)
    const { atomIndex, startPos } = dragState;
    const cur = molState.positions[atomIndex];
    const moved = Math.abs(cur.x-startPos.x)>1e-9 || Math.abs(cur.y-startPos.y)>1e-9 || Math.abs(cur.z-startPos.z)>1e-9;
    if (DRAG_LOG) {
      const dx = +(cur.x - startPos.x).toFixed(3);
      const dy = +(cur.y - startPos.y).toFixed(3);
      const dz = +(cur.z - startPos.z).toFixed(3);
      console.log('[drag][drop]', { atomIndex, moved, delta: { dx, dy, dz }, startPos: { ...startPos }, finalPos: { ...cur } });
    } // suppressed by default
    dragState = null;
    if (moved && bondService && typeof bondService.recomputeAndStore === 'function') {
      bondService.recomputeAndStore();
    }
  }

  // Bond rotation: rotate all atoms on one side around axis defined by bond i-j passing through atom i (if side==='j') or j (if side==='i').
  function rotateBond(deltaAngle, sideOverride) {
    if (!molState.selection || molState.selection.kind !== 'bond') return false;
    const { i, j, orientation } = molState.selection.data;
    const side = sideOverride || orientationToSide(orientation); // 'i' or 'j'
    const anchor = side === 'i' ? j : i; // rotate the opposite side around anchor->moving axis
    const movingRoot = side === 'i' ? i : j; // starting atom of moving side
    // Build adjacency
    const adj = Array.from({ length: molState.positions.length }, () => []);
    for (const b of molState.bonds) { adj[b.i].push(b.j); adj[b.j].push(b.i); }
    // Collect side atoms via BFS starting at movingRoot but do not cross the anchor
    const sideAtoms = [];
    const queue=[movingRoot];
    const visited=new Set([anchor]);
    visited.add(movingRoot);
    while(queue.length) {
      const a=queue.shift(); sideAtoms.push(a);
      for (const nb of adj[a]) if (!visited.has(nb)) { visited.add(nb); queue.push(nb); }
    }
    // Axis from anchor to movingRoot
    const pAnchor = molState.positions[anchor];
    const pMove = molState.positions[movingRoot];
    const ax = pMove.x - pAnchor.x; const ay = pMove.y - pAnchor.y; const az = pMove.z - pAnchor.z;
    const len = Math.hypot(ax,ay,az) || 1;
    const ux = ax/len, uy = ay/len, uz = az/len;
    const s = Math.sin(deltaAngle), c = Math.cos(deltaAngle);
    function rotatePoint(p) {
      // Translate
      const dx = p.x - pAnchor.x; const dy = p.y - pAnchor.y; const dz = p.z - pAnchor.z;
      // Rodrigues rotation formula
      const dot = dx*ux + dy*uy + dz*uz;
      const rx = ux*dot + (dx - ux*dot)*c + (-uz*dy + uy*dz)*s;
      const ry = uy*dot + (dy - uy*dot)*c + (uz*dx - ux*dz)*s;
      const rz = uz*dot + (dz - uz*dot)*c + (-uy*dx + ux*dy)*s;
      p.x = pAnchor.x + rx; p.y = pAnchor.y + ry; p.z = pAnchor.z + rz;
    }
    for (const idx of sideAtoms) rotatePoint(molState.positions[idx]);
    molState.markPositionsChanged();
    if (bondService && typeof bondService.recomputeAndStore === 'function') {
      bondService.recomputeAndStore();
    }
    return true;
  }

  return { beginDrag, updateDrag, endDrag, setDragPlane, rotateBond };
}
