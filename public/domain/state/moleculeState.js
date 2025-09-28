// MoleculeState and EventBus (Phase 1 scaffold)
// Pure data container + minimal mutation helpers + pub/sub.

export function createEventBus() {
  const listeners = new Map(); // event -> Set(fn)
  return {
    on(event, fn) { if (!listeners.has(event)) listeners.set(event, new Set()); listeners.get(event).add(fn); return () => listeners.get(event)?.delete(fn); },
    emit(event, payload) { const set = listeners.get(event); if (set) for (const fn of Array.from(set)) { try { fn(payload); } catch(e){ /* swallow for stability */ } } },
    clear() { listeners.clear(); }
  };
}

export function createMoleculeState(initial = {}) {
  const state = {
    atoms: initial.atoms || [], // [{ element, pos: BABYLON.Vector3, id }]
    bonds: initial.bonds || [], // [{ i,j,weight,opacity,flags }]
    cell: initial.cell || { a:null, b:null, c:null, originOffset:null, visible:false },
    selection: initial.selection || { kind:null, data:{} },
    dynamics: initial.dynamics || { mode:'idle', step:0, temperature:0 },
    versions: { atoms:0, bonds:0, cell:0, selection:0, dynamics:0 },
    meta: initial.meta || {},
    bus: createEventBus()
  };
  return state;
}

// Mutation helpers; each increments version + emits event for render/services.
export const StateMutators = {
  setAtoms(state, atoms) {
    state.atoms = atoms; state.versions.atoms++; state.bus.emit('atomsChanged', { atoms });
  },
  updateAtomPositions(state, positions) { // positions: [{id,pos}]
    const byId = new Map(state.atoms.map(a=>[a.id,a]));
    for (const { id, pos } of positions) { const a = byId.get(id); if (a) a.pos = pos; }
    state.versions.atoms++; state.bus.emit('atomsChanged', { atoms: state.atoms });
  },
  setBonds(state, bonds) { state.bonds = bonds; state.versions.bonds++; state.bus.emit('bondsChanged', { bonds }); },
  setCell(state, cellPatch) { state.cell = { ...state.cell, ...cellPatch }; state.versions.cell++; state.bus.emit('cellChanged', { cell: state.cell }); },
  setSelection(state, sel) { state.selection = sel; state.versions.selection++; state.bus.emit('selectionChanged', { selection: sel }); },
  setDynamics(state, patch) { state.dynamics = { ...state.dynamics, ...patch }; state.versions.dynamics++; state.bus.emit('dynamicsChanged', { dynamics: state.dynamics }); }
};

export function getDirtyFlags(prevVersions, current) {
  const dirty = {};
  for (const k of Object.keys(prevVersions)) { if (prevVersions[k] !== current.versions[k]) dirty[k] = true; }
  return dirty; // e.g., { atoms:true, bonds:true }
}
