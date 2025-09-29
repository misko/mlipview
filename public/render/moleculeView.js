import { elInfo } from '../elements.js';
import { computeBondsNoState } from '../bond_render.js';

export function createMoleculeView(scene, molState) {
  const atomGroups = new Map();
  const bondGroups = new Map();
  const ghostAtomGroups = new Map(); // separate to keep pickable flag false
  const ghostBondGroups = new Map();
  let cellLines = null; // line system for cell outline
  // Highlight meshes (single reusable shell each for selected atom or bond)
  const highlight = { atom:null, bond:null };
  function keyOf(i,j){ const a=molState.elements[i], b=molState.elements[j]; return [a,b].sort().join('-'); }
  function ensureAtomGroup(el) {
    if (atomGroups.has(el)) return atomGroups.get(el);
    const info = elInfo(el);
    const mat = new BABYLON.StandardMaterial('mat_'+el, scene);
    mat.diffuseColor = info.color.clone();
    mat.emissiveColor = info.color.scale(0.06);
    const master = BABYLON.MeshBuilder.CreateSphere('atom_'+el,{diameter:1,segments:24},scene);
    master.material = mat; master.isPickable = true; master.thinInstanceEnablePicking = true;
    // Keep master enabled; thin instances render on top but some pipelines rely on master visibility for picking bounds.
    if (typeof master.setEnabled === 'function') master.setEnabled(true); else master.isVisible = true;
    const g = { master, mats:[], indices:[] }; atomGroups.set(el,g); return g;
  }
  function ensureGhostAtomGroup(el) {
    if (ghostAtomGroups.has(el)) return ghostAtomGroups.get(el);
    const info = elInfo(el);
    const mat = new BABYLON.StandardMaterial('ghost_atom_'+el, scene);
    mat.diffuseColor = info.color.clone();
    mat.emissiveColor = info.color.scale(0.02);
    mat.alpha = 0.35;
    const master = BABYLON.MeshBuilder.CreateSphere('ghost_atom_'+el,{diameter:1,segments:16},scene);
    master.material = mat; master.isPickable = false; master.thinInstanceEnablePicking = false;
    const g = { master, mats:[], indices:[] }; ghostAtomGroups.set(el,g); return g;
  }
  function ensureGhostBondGroup(key) {
    if (ghostBondGroups.has(key)) return ghostBondGroups.get(key);
    const mat = new BABYLON.StandardMaterial('ghost_bond_'+key, scene);
    mat.diffuseColor = new BABYLON.Color3(0.55,0.58,0.60);
    mat.emissiveColor = mat.diffuseColor.scale(0.02);
    mat.alpha = 0.25;
    const master = BABYLON.MeshBuilder.CreateCylinder('ghost_bond_'+key,{height:1,diameter:1,tessellation:12},scene);
    master.material=mat; master.isPickable=false; master.thinInstanceEnablePicking=false;
    const g = { master, mats:[], indices:[] }; ghostBondGroups.set(key,g); return g;
  }
  function ensureBondGroup(key) {
    if (bondGroups.has(key)) return bondGroups.get(key);
    const mat = new BABYLON.StandardMaterial('bond_'+key, scene);
    mat.diffuseColor = new BABYLON.Color3(0.75,0.78,0.80);
    mat.emissiveColor = mat.diffuseColor.scale(0.05);
  const master = BABYLON.MeshBuilder.CreateCylinder('bond_'+key,{height:1,diameter:1,tessellation:22},scene);
  // IMPORTANT: bonds must be pickable for bond selection highlight (was false previously, preventing bond highlight in runtime)
  master.material=mat; master.isPickable=true; master.thinInstanceEnablePicking=true;
    const g = { master, mats:[], indices:[] }; bondGroups.set(key,g); return g;
  }
  function buildInitial() {
    for (let i=0;i<molState.positions.length;i++) {
      const el = molState.elements[i];
      const g = ensureAtomGroup(el);
      const info = elInfo(el);
      const d = info.scale;
      const p = molState.positions[i];
      const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(d,d,d), BABYLON.Quaternion.Identity(), new BABYLON.Vector3(p.x,p.y,p.z));
      g.mats.push(mat); g.indices.push(i);
    }
    for (const g of atomGroups.values()) g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    rebuildBonds();
  }
  function flattenMatrices(mats){ const arr=new Float32Array(mats.length*16); mats.forEach((m,i)=>{arr.set(m.m,i*16);}); return arr; }
  function bondMatrix(pA,pB,radius){
    const mid = new BABYLON.Vector3((pA.x+pB.x)/2,(pA.y+pB.y)/2,(pA.z+pB.z)/2);
    const v = new BABYLON.Vector3(pB.x-pA.x,pB.y-pA.y,pB.z-pA.z);
    const len = v.length(); const up = new BABYLON.Vector3(0,1,0);
    let rot; const d=v.normalizeToNew(); const dot=BABYLON.Vector3.Dot(up,d);
    if (dot>0.9999) rot=BABYLON.Quaternion.Identity(); else if (dot<-0.9999) rot=BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1,0,0),Math.PI); else { const axis=BABYLON.Vector3.Cross(up,d).normalize(); rot=BABYLON.Quaternion.RotationAxis(axis,Math.acos(dot)); }
    return BABYLON.Matrix.Compose(new BABYLON.Vector3(radius*2,len,radius*2), rot, mid);
  }
  function rebuildBonds(bondData) {
    for (const g of bondGroups.values()) { g.mats=[]; g.indices=[]; }
    const source = bondData || molState.bonds;
    const debug = (typeof window !== 'undefined' && window.location && window.location.search.includes('debug=1'));
    for (const b of source) {
      const g = ensureBondGroup(keyOf(b.i,b.j));
      const pA = molState.positions[b.i]; const pB = molState.positions[b.j];
      const mat = bondMatrix(pA,pB,0.1);
      g.mats.push(mat); g.indices.push(b);
    }
    for (const [key,g] of bondGroups) {
      g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
      // Build or reuse color buffer
      const cols = new Float32Array(g.mats.length * 4);
      for (let i=0;i<g.indices.length;i++) {
        const b = g.indices[i];
        const alpha = (b.opacity != null) ? b.opacity : 1.0;
        // Base diffuse color pulled from material diffuseColor
        const dc = g.master.material.diffuseColor || { r:0.75,g:0.78,b:0.80 };
        cols[i*4+0] = dc.r;
        cols[i*4+1] = dc.g;
        cols[i*4+2] = dc.b;
        cols[i*4+3] = alpha;
      }
      // Assign color buffer (Babylon expects 'color' for per-instance colors)
      try { g.master.thinInstanceSetBuffer('color', cols, 4); } catch(e) { /* ignore in stub */ }
      // Ensure material blending if any alpha < 1
      if (g.indices.some(b=> b.opacity != null && b.opacity < 0.999)) {
        g.master.material.alpha = 1.0; // keep host opaque; per-instance alpha used
        if (BABYLON.Material && typeof BABYLON.Material.MATERIAL_ALPHABLEND !== 'undefined') {
          g.master.material.transparencyMode = BABYLON.Material.MATERIAL_ALPHABLEND;
        }
      } else {
        g.master.material.alpha = 1.0;
      }
      // Hide master mesh when it has zero instances to avoid stray default cylinder (legacy behavior parity)
      if (g.mats.length === 0) {
        if (typeof g.master.setEnabled === 'function') g.master.setEnabled(false); else g.master.isVisible=false;
      } else {
        if (typeof g.master.setEnabled === 'function') g.master.setEnabled(true); else g.master.isVisible=true;
      }
      if (debug) console.debug('[debug][rebuildBonds] group', key, 'instances', g.mats.length);
    }
  }
  function rebuildGhosts() {
    // Clear previous ghost buffers
    for (const g of ghostAtomGroups.values()) { g.mats.length = 0; g.indices.length = 0; }
    for (const g of ghostBondGroups.values()) { g.mats.length = 0; g.indices.length = 0; }
    if (!molState.showGhostCells || !molState.showCell || !molState.cell?.enabled) {
      for (const g of ghostAtomGroups.values()) g.master.thinInstanceSetBuffer('matrix', new Float32Array());
      for (const g of ghostBondGroups.values()) g.master.thinInstanceSetBuffer('matrix', new Float32Array());
      return;
    }
    const { a, b, c } = molState.cell;
    const shifts = [ {x:0,y:0,z:0}, {x:1,y:0,z:0}, {x:-1,y:0,z:0}, {x:0,y:1,z:0}, {x:0,y:-1,z:0}, {x:0,y:0,z:1}, {x:0,y:0,z:-1} ];
    // Build augmented atom list (first block is original atoms at shift 0)
    const augAtoms = [];
    for (const S of shifts) {
      for (let i=0;i<molState.positions.length;i++) {
        const base = molState.positions[i];
        const pos = { x: base.x + S.x*a.x + S.y*b.x + S.z*c.x, y: base.y + S.x*a.y + S.y*b.y + S.z*c.y, z: base.z + S.x*a.z + S.y*b.z + S.z*c.z };
        augAtoms.push({ element: molState.elements[i], baseIndex:i, shift:S, pos });
        if (S.x!==0 || S.y!==0 || S.z!==0) {
          // Only create thin instances for non-primary images (primary already rendered in main atom groups)
          const g = ensureGhostAtomGroup(molState.elements[i]);
          const info = elInfo(molState.elements[i]); const d = info.scale;
          const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(d,d,d), BABYLON.Quaternion.Identity(), new BABYLON.Vector3(pos.x,pos.y,pos.z));
          g.mats.push(mat); g.indices.push({ base:i, shift:[S.x,S.y,S.z] });
        }
      }
    }
    for (const g of ghostAtomGroups.values()) g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    // Run bond calculator over augmented list (stateless) and then extract bonds that involve at least one ghost image OR cross-image pair.
    // We import lazily to avoid circular dependency; computeBondsNoState is already globally available in legacy path.
    const augSimple = augAtoms.map(a=>({ element:a.element, pos:[a.pos.x,a.pos.y,a.pos.z] }));
    const augBonds = computeBondsNoState(augSimple);
    // Build mapping from augmented atom index back to (baseIndex, shift)
    // augmented index layout: for shift index si and local atom i => idx = si * nAtoms + i
    const nAtoms = molState.positions.length;
    function decode(idx){ const si = Math.floor(idx / nAtoms); const local = idx % nAtoms; return { si, local, shift: shifts[si] }; }
    for (const eb of augBonds) {
      const A = decode(eb.i); const B = decode(eb.j);
      // Skip primary-primary (both shift 0) because those are already rendered as real bonds
      if (A.si===0 && B.si===0) continue;
      const baseI = A.local; const baseJ = B.local;
      const key = keyOf(baseI, baseJ);
      const g = ensureGhostBondGroup(key);
      const pA = augAtoms[eb.i].pos; const pB = augAtoms[eb.j].pos;
      const mat = bondMatrix(pA,pB,0.1);
      g.mats.push(mat); g.indices.push({ i:baseI, j:baseJ, shiftA:[A.shift.x,A.shift.y,A.shift.z], shiftB:[B.shift.x,B.shift.y,B.shift.z] });
    }
    for (const g of ghostBondGroups.values()) g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
  }
  function rebuildCellLines() {
    if (!molState.showCell || !molState.cell?.enabled) { if (cellLines) { cellLines.dispose(); cellLines=null; } return; }
    const { a,b,c, originOffset } = molState.cell;
    const O = originOffset || {x:0,y:0,z:0};
    const pts = [
      O,
      {x:O.x+a.x,y:O.y+a.y,z:O.z+a.z},
      {x:O.x+a.x+b.x,y:O.y+a.y+b.y,z:O.z+a.z+b.z},
      {x:O.x+b.x,y:O.y+b.y,z:O.z+b.z},
      O,
      {x:O.x+c.x,y:O.y+c.y,z:O.z+c.z},
      {x:O.x+a.x+c.x,y:O.y+a.y+c.y,z:O.z+a.z+c.z},
      {x:O.x+a.x+b.x+c.x,y:O.y+a.y+b.y+c.y,z:O.z+a.z+b.z+c.z},
      {x:O.x+b.x+c.x,y:O.y+b.y+c.y,z:O.z+b.z+c.z},
      {x:O.x+c.x,y:O.y+c.y,z:O.z+c.z},
      {x:O.x+a.x+c.x,y:O.y+a.y+c.y,z:O.z+a.z+c.z},
      {x:O.x+a.x+b.x+c.x,y:O.y+a.y+b.y+c.y,z:O.z+a.z+b.z+c.z},
      {x:O.x+a.x+b.x,y:O.y+a.y+b.y,z:O.z+a.z+b.z},
      {x:O.x+b.x+c.x,y:O.y+b.y+c.y,z:O.z+b.z+c.z},
      {x:O.x+b.x,y:O.y+b.y,z:O.z+b.z}
    ].map(p=> new BABYLON.Vector3(p.x,p.y,p.z));
    if (cellLines) cellLines.dispose();
    cellLines = BABYLON.MeshBuilder.CreateLines('cell_lines',{ points: pts }, scene);
    cellLines.color = new BABYLON.Color3(0.9,0.8,0.25);
    cellLines.isPickable = false;
  }
  function rebuildAtoms() {
    // Clear current groups and rebuild from molState.elements/positions
    for (const g of atomGroups.values()) { g.mats.length = 0; g.indices.length = 0; }
    for (let i=0;i<molState.positions.length;i++) {
      const el = molState.elements[i];
      const g = ensureAtomGroup(el);
      const info = elInfo(el); const d = info.scale; const p = molState.positions[i];
      const mat = BABYLON.Matrix.Compose(new BABYLON.Vector3(d,d,d), BABYLON.Quaternion.Identity(), new BABYLON.Vector3(p.x,p.y,p.z));
      g.mats.push(mat); g.indices.push(i);
    }
    for (const g of atomGroups.values()) g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    const debug = (typeof window !== 'undefined' && window.location && window.location.search.includes('debug=1'));
    if (debug) console.debug('[debug][rebuildAtoms] groups', Array.from(atomGroups.entries()).map(([el,g])=>({el, count:g.mats.length, visible:g.master.isVisible})));
  }
  function updatePositions() {
    for (const [el,g] of atomGroups) {
      for (let k=0;k<g.indices.length;k++) {
        const idx = g.indices[k];
        const p = molState.positions[idx];
        if (!p) continue; // guard during transient molecule swap
        const info = elInfo(el); const d = info.scale;
        g.mats[k] = BABYLON.Matrix.Compose(new BABYLON.Vector3(d,d,d), BABYLON.Quaternion.Identity(), new BABYLON.Vector3(p.x,p.y,p.z));
      }
      g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    }
  }
  molState.bus.on('positionsChanged', updatePositions);
  molState.bus.on('bondsChanged', () => {
    // If atom groups length total != positions length, we loaded a new molecule; rebuild atoms first.
    const currentCount = Array.from(atomGroups.values()).reduce((s,g)=>s+g.indices.length,0);
    const topologyChanged = currentCount !== molState.positions.length;
    if (topologyChanged) rebuildAtoms();
    rebuildBonds();
    rebuildGhosts();
    // Decide if current selection is still valid
    let invalidate = false;
    if (molState.selection && molState.selection.kind) {
      if (molState.selection.kind === 'atom') {
        const idx = molState.selection.data?.index;
        if (idx == null || idx < 0 || idx >= molState.positions.length) invalidate = true;
      } else if (molState.selection.kind === 'bond') {
        const { i, j } = molState.selection.data || {};
        if (i == null || j == null || i < 0 || j < 0 || i >= molState.positions.length || j >= molState.positions.length) {
          invalidate = true;
        } else {
          // Bond must still exist (unordered match)
          const exists = molState.bonds.some(b => (b.i===i && b.j===j) || (b.i===j && b.j===i));
          if (!exists) invalidate = true;
        }
      }
    }
    if (invalidate || topologyChanged) {
      molState.selection = { kind:null, data:null };
      if (highlight.bond) highlight.bond.isVisible=false;
      if (highlight.atom) highlight.atom.isVisible=false;
    } else {
      // Selection remains; refresh highlight transform/visibility
      updateSelectionHighlight();
    }
    const debug = (typeof window !== 'undefined' && window.location && window.location.search.includes('debug=1'));
    if (debug) {
      console.debug('[debug][bondsChanged] atomCount', molState.positions.length, 'groups', Array.from(atomGroups.keys()));
      for (const [el,g] of atomGroups) {
        console.debug('[debug][atomGroup]', el, { instances:g.mats.length, masterVisible:g.master.isVisible, hasColorBuffer: !!g.master?._thinInstanceBufferUpdated });
      }
      console.debug('[debug][highlight]', { atomVisible: highlight.atom?.isVisible, bondVisible: highlight.bond?.isVisible });
    }
  });
  molState.bus.on('cellChanged', () => { rebuildCellLines(); rebuildGhosts(); });
  // Selection highlight management
  function ensureHighlightMeshes() {
    if (!highlight.atom) {
      const mat = new BABYLON.StandardMaterial('highlight_atom', scene);
      mat.diffuseColor = new BABYLON.Color3(0.0, 0.9, 0.95); // cyan tint
      mat.emissiveColor = new BABYLON.Color3(0.0, 0.45, 0.5);
      mat.alpha = 0.45; // within spec test range (0.3 - 0.6)
      mat.disableLighting = true;
      mat.backFaceCulling = false;
      const sphere = BABYLON.MeshBuilder.CreateSphere('highlight_atom_mesh', { diameter:1, segments:24 }, scene);
      sphere.material = mat;
      sphere.isPickable = false;
      sphere.isVisible = false;
      sphere.alwaysSelectAsActiveMesh = true;
      sphere.renderingGroupId = 2; // render after atoms
      highlight.atom = sphere;
    }
    if (!highlight.bond) {
      const mat = new BABYLON.StandardMaterial('highlight_bond', scene);
      mat.diffuseColor = new BABYLON.Color3(0.0, 0.9, 0.95);
      mat.emissiveColor = new BABYLON.Color3(0.0, 0.45, 0.5);
      mat.alpha = 0.40; // similar look to atom highlight, slightly less opaque
      mat.disableLighting = true;
      mat.backFaceCulling = false;
      const cyl = BABYLON.MeshBuilder.CreateCylinder('highlight_bond_mesh', { height:1, diameter:1, tessellation:22 }, scene);
      cyl.material = mat;
      cyl.isPickable = false;
      // Keep fully disabled until an actual bond is selected to avoid stray cylinder at origin
      cyl.isVisible = false;
      cyl.alwaysSelectAsActiveMesh = true;
      cyl.renderingGroupId = 2;
      highlight.bond = cyl;
    }
  }
  function updateSelectionHighlight() {
    ensureHighlightMeshes();
    const sel = molState.selection;
    // Hide both first
    highlight.atom.isVisible = false;
    highlight.bond.isVisible = false;
    // Also disable bond mesh so it cannot render with default scale at origin
  // Avoid relying on setEnabled (in some stubs it may be missing); visibility alone controls rendering here.
    if (!sel || !sel.kind) return;
    if (sel.kind === 'atom') {
      const idx = sel.data.index;
      const p = molState.positions[idx];
      const el = molState.elements[idx];
      const info = elInfo(el);
      // Slight enlargement to form a subtle shell (previously 2.2x which was visually bulky)
      const scale = info.scale * 1.35;
      // IMPORTANT (VR rotation fix): parent the highlight sphere to the atom master mesh so
      // it inherits any scene-level rotations applied to masters (e.g. via WebXR controller).
      // Previously the highlight sphere lived at world coordinates and did not rotate with
      // the molecule, causing misalignment after user rotations in VR.
      const group = atomGroups.get(el);
      if (group && group.master) {
        // Only (re)parent if different to avoid triggering Babylon matrix recomputation unnecessarily.
        if (highlight.atom.parent !== group.master) {
          highlight.atom.parent = group.master;
        }
        // Because we parent to the master (whose thin instances represent atoms), the highlight
        // position must be expressed in the master's local coordinate system (identical to the
        // atom's original position vector). Rotation/scaling of the master will now also affect
        // the highlight automatically.
        highlight.atom.position = new BABYLON.Vector3(p.x, p.y, p.z);
      } else {
        // Fallback: no master found (should not happen) – keep unparented world position path.
        if (highlight.atom.parent) highlight.atom.parent = null;
        highlight.atom.position = new BABYLON.Vector3(p.x, p.y, p.z);
      }
      highlight.atom.scaling = new BABYLON.Vector3(scale, scale, scale);
      highlight.atom.isVisible = true;
    } else if (sel.kind === 'bond') {
      const { i, j } = sel.data;
      const pA = molState.positions[i];
      const pB = molState.positions[j];
      if (!pA || !pB) return;
      const mid = new BABYLON.Vector3((pA.x+pB.x)/2,(pA.y+pB.y)/2,(pA.z+pB.z)/2);
      const v = new BABYLON.Vector3(pB.x-pA.x,pB.y-pA.y,pB.z-pA.z);
      const len = v.length();
      const up = new BABYLON.Vector3(0,1,0);
      let rotQ; const d=v.normalizeToNew(); const dot=BABYLON.Vector3.Dot(up,d);
      if (dot>0.9999) rotQ=BABYLON.Quaternion.Identity(); else if (dot<-0.9999) rotQ=BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1,0,0),Math.PI); else { const axis=BABYLON.Vector3.Cross(up,d).normalize(); rotQ=BABYLON.Quaternion.RotationAxis(axis,Math.acos(dot)); }
      highlight.bond.position = mid;
      highlight.bond.rotationQuaternion = rotQ;
      const radius = 0.16; // proportionate shell around default 0.1 bond radius
      highlight.bond.scaling = new BABYLON.Vector3(radius*2, len, radius*2);
      highlight.bond.isVisible = true;
    }
  }
  molState.bus.on('selectionChanged', updateSelectionHighlight);
  // In some minimal test environments, event loop timing can differ; proactively call once if a selection already exists.
  if (molState.selection && molState.selection.kind) { try { updateSelectionHighlight(); } catch(e){} }
  // If positions of a selected item change, update highlight transform.
  molState.bus.on('positionsChanged', () => {
    if (molState.selection && molState.selection.kind) updateSelectionHighlight();
  });
  // Frame-level safety: ensure we never show a highlight when selection is null/invalid.
  // This guards against race conditions during molecule switches leaving a stray mesh enabled.
  // Some test environments pass a minimal scene stub without onBeforeRenderObservable.
  // Provide a lightweight shim if absent so highlight safety still functions.
  if (!scene.onBeforeRenderObservable) {
    scene.onBeforeRenderObservable = { _c:[], add(fn){ this._c.push(fn); }, run(){ for(const f of this._c) try{f();}catch(_){} } };
  }
  scene.onBeforeRenderObservable.add(() => {
    if (!highlight.atom || !highlight.bond) return; // not yet created
    const sel = molState.selection;
    if (!sel || !sel.kind) {
      if (highlight.atom.isVisible) highlight.atom.isVisible = false;
      if (highlight.bond.isVisible) highlight.bond.isVisible = false;
      return;
    }
    // Additional integrity checks: indices must be in range, bond must still exist
    if (sel.kind === 'atom') {
      const idx = sel.data?.index;
      if (!(Number.isInteger(idx) && idx >= 0 && idx < molState.positions.length)) {
        highlight.atom.isVisible = false;
        highlight.bond.isVisible = false;
        molState.selection = { kind:null, data:null };
      }
    } else if (sel.kind === 'bond') {
      const { i, j } = sel.data || {};
      const valid = Number.isInteger(i) && Number.isInteger(j) && i >=0 && j >=0 && i < molState.positions.length && j < molState.positions.length && molState.bonds.some(b => (b.i===i && b.j===j) || (b.i===j && b.j===i));
      if (!valid) {
        highlight.atom.isVisible = false;
        highlight.bond.isVisible = false;
        molState.selection = { kind:null, data:null };
      }
    }
  });
  buildInitial();
  // Initialize highlight once after initial build
  ensureHighlightMeshes();
  function resolveAtomPick(pick) {
    if (!pick?.hit || !pick.pickedMesh) return null;
    for (const [el,g] of atomGroups) {
      if (g.master === pick.pickedMesh) {
        const idx = pick.thinInstanceIndex;
        if (idx == null || idx < 0 || idx >= g.indices.length) return null;
        const atomIndex = g.indices[idx];
        return { kind:'atom', index: atomIndex, element: el };
      }
    }
    return null;
  }
  function resolveBondPick(pick) {
    if (!pick?.hit || !pick.pickedMesh) return null;
    for (const [key,g] of bondGroups) {
      if (g.master === pick.pickedMesh) {
        const idx = pick.thinInstanceIndex;
        if (idx == null || idx < 0 || idx >= g.indices.length) return null;
        const b = g.indices[idx];
        return { kind:'bond', i:b.i, j:b.j, key, index: idx };
      }
    }
    return null;
  }
  return { rebuildBonds, rebuildGhosts, _internals:{ atomGroups, bondGroups, ghostAtomGroups, ghostBondGroups, highlight }, resolveAtomPick, resolveBondPick };
}
