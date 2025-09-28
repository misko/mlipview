import { elInfo } from '../elements.js';

export function createMoleculeView(scene, molState) {
  const atomGroups = new Map();
  const bondGroups = new Map();
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
    const g = { master, mats:[], indices:[] }; atomGroups.set(el,g); return g;
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
    // Collect per-group color arrays (RGBA) aligned with g.mats
    const colorsPerGroup = new Map();
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
    }
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
  }
  function updatePositions() {
    for (const [el,g] of atomGroups) {
      for (let k=0;k<g.indices.length;k++) {
        const idx = g.indices[k];
        const p = molState.positions[idx];
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
    if (currentCount !== molState.positions.length) rebuildAtoms();
    rebuildBonds();
  });
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
    if (!sel || !sel.kind) return;
    if (sel.kind === 'atom') {
      const idx = sel.data.index;
      const p = molState.positions[idx];
      const el = molState.elements[idx];
      const info = elInfo(el);
      // Slight enlargement to form a subtle shell (previously 2.2x which was visually bulky)
      const scale = info.scale * 1.35;
      highlight.atom.position = new BABYLON.Vector3(p.x, p.y, p.z);
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
  // If positions of a selected item change, update highlight transform.
  molState.bus.on('positionsChanged', () => {
    if (molState.selection && molState.selection.kind) updateSelectionHighlight();
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
  return { rebuildBonds, _internals:{ atomGroups, bondGroups, highlight }, resolveAtomPick, resolveBondPick };
}
