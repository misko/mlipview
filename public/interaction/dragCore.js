// mlipviewer2/public/interaction/dragCore.js
// Platform-agnostic atom drag service. Adapters provide world points.
// Emits events through provided bus (expected API: bus.emit(name, payload)).

export function createDragCore({ moleculeState, selectionService, manipulationService, bus }) {
  // Internal tracking of active drag
  let active = null; // { atomIndex, startPos: {x,y,z} }

  function begin(atomIndex) {
    if (active) end(false);
    if (atomIndex == null || atomIndex < 0 || atomIndex >= moleculeState.positions.length) return false;
    active = { atomIndex, startPos: { ...moleculeState.positions[atomIndex] } };
    // Select atom if not already selected
    if (!selectionService.getSelection() || selectionService.getSelection().kind !== 'atom' || selectionService.getSelection().data.index !== atomIndex) {
      selectionService.selectAtom(atomIndex);
    }
    bus.emit('dragStarted', { atomIndex });
    return true;
  }

  function move(worldPoint) {
    if (!active) return;
    const { atomIndex } = active;
    const p = moleculeState.positions[atomIndex];
    if (!p) return;
    p.x = worldPoint.x; p.y = worldPoint.y; p.z = worldPoint.z;
    // Mark positions changed (downstream: moleculeView updates matrices)
    moleculeState.bus.emit('positionsChanged');
    bus.emit('dragMoved', { atomIndex });
  }

  function end(commit = true) {
    if (!active) return;
    const { atomIndex, startPos } = active;
    if (!commit) {
      // Revert
      const p = moleculeState.positions[atomIndex];
      p.x = startPos.x; p.y = startPos.y; p.z = startPos.z;
      moleculeState.bus.emit('positionsChanged');
    } else {
      // Optionally mark changed for physics caches
      if (typeof moleculeState.markPositionsChanged === 'function') moleculeState.markPositionsChanged();
    }
    bus.emit('dragEnded', { atomIndex, commit });
    active = null;
  }

  function isDragging() { return !!active; }

  return { begin, move, end, isDragging };
}
