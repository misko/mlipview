import { elInfo } from '../elements.js';
import { computeBondsNoState } from '../bond_render.js';
import { __count } from '../util/funcCount.js';

export function createMoleculeView(scene, molState) {
  __count('moleculeView#createMoleculeView');
  // === Molecule Master Root & Registry =====================================
  // We introduce a single transform node (moleculeRoot) that becomes the parent
  // of every master thin‑instance host (atoms, bonds, forces, ghosts, highlights).
  // VR / AR systems should manipulate ONLY this root when applying global
  // rotations / translations / scaling. This eliminates regex heuristics over
  // mesh names and provides a stable anchor even if internal naming changes.
  //
  // A lightweight registry is also exposed on molState.__masters for debugging
  // and for fallback consumers. Each entry is { kind, mesh }. The registry is
  // append‑only during a molecule’s lifetime; on a full molecule swap a new
  // view (with a new root & registry) is constructed.
  // ===========================================================================
  // Some test environments provide a very small Babylon stub without TransformNode.
  // Fallback: use a MeshBuilder sphere (disabled) as the root parent container.
  let moleculeRoot;
  try {
    if (BABYLON.TransformNode) {
      moleculeRoot = new BABYLON.TransformNode('molecule_root', scene);
    } else {
      // minimal fallback; flagged via metadata for VR utilities
      moleculeRoot = BABYLON.MeshBuilder?.CreateSphere ? BABYLON.MeshBuilder.CreateSphere('molecule_root',{diameter:0.0001},scene) : { name:'molecule_root', position:new BABYLON.Vector3(0,0,0), scaling:new BABYLON.Vector3(1,1,1) };
      if (moleculeRoot.setEnabled) moleculeRoot.setEnabled(false); else moleculeRoot.isVisible=false;
    }
  } catch {
    moleculeRoot = { name:'molecule_root_stub', position:{x:0,y:0,z:0}, scaling:{x:1,y:1,z:1} };
  }
  moleculeRoot.metadata = { role: 'moleculeRoot' };
  const mastersRegistry = []; // array of { kind: 'atomMaster'|'bondMaster'|'forceMaster'|'ghostAtomMaster'|'ghostBondMaster', mesh }
  function registerMaster(kind, mesh) {
    try { if (mesh && mesh.parent !== moleculeRoot) mesh.parent = moleculeRoot; } catch {}
    mastersRegistry.push({ kind, mesh });
  }
  // Expose for VR/services (read‑only intended):
  molState.moleculeRoot = moleculeRoot;
  molState.__masters = mastersRegistry;
  const atomGroups = new Map();
  const bondGroups = new Map();
  // Force vectors rendered as thin instances (similar to bonds) for VR/AR consistency
  const forceGroups = new Map(); // single group keyed by 'force' currently; extensible for per-type later
  const ghostAtomGroups = new Map(); // separate to keep pickable flag false
  const ghostBondGroups = new Map();
  // Geometry version stamping: increment outside (molState.geometryVersion++) when topology changes
  if (typeof molState.geometryVersion !== 'number') molState.geometryVersion = 0;
  let cachedAtomVersion = -1;
  let cachedBondVersion = -1;
  let cellLines = null; // line system for cell outline
  // Highlight meshes (single reusable shell each for selected atom or bond)
  const highlight = { atom:null, bond:null };
  function keyOf(i,j){ const a=molState.elements[i], b=molState.elements[j]; return [a,b].sort().join('-'); }
  function ensureGroup(map, key, kind){
    // kind: 'atom' | 'ghostAtom' | 'bond' | 'ghostBond'
    const store = map;
    if (store.has(key)) return store.get(key);
    let mat, master;
    if (kind === 'atom') {
      const info = elInfo(key);
      mat = new BABYLON.StandardMaterial('mat_'+key, scene);
      mat.diffuseColor = info.color.clone();
      mat.emissiveColor = info.color.scale(0.06);
      master = BABYLON.MeshBuilder.CreateSphere('atom_'+key,{diameter:1,segments:24},scene);
      master.isPickable = true; master.thinInstanceEnablePicking = true;
      if (typeof master.setEnabled === 'function') master.setEnabled(true); else master.isVisible = true;
      registerMaster('atomMaster', master);
    } else if (kind === 'ghostAtom') {
      const info = elInfo(key);
      mat = new BABYLON.StandardMaterial('ghost_atom_'+key, scene);
      mat.diffuseColor = info.color.clone();
      mat.emissiveColor = info.color.scale(0.02);
      mat.alpha = 0.35;
      master = BABYLON.MeshBuilder.CreateSphere('ghost_atom_'+key,{diameter:1,segments:16},scene);
      master.isPickable=false; master.thinInstanceEnablePicking=false;
      registerMaster('ghostAtomMaster', master);
    } else if (kind === 'ghostBond') {
      mat = new BABYLON.StandardMaterial('ghost_bond_'+key, scene);
      mat.diffuseColor = new BABYLON.Color3(0.55,0.58,0.60);
      mat.emissiveColor = mat.diffuseColor.scale(0.02); mat.alpha = 0.25;
      master = BABYLON.MeshBuilder.CreateCylinder('ghost_bond_'+key,{height:1,diameter:1,tessellation:12},scene);
      master.isPickable=false; master.thinInstanceEnablePicking=false;
      registerMaster('ghostBondMaster', master);
    } else if (kind === 'bond') {
      mat = new BABYLON.StandardMaterial('bond_'+key, scene);
      mat.diffuseColor = new BABYLON.Color3(0.75,0.78,0.80);
      mat.emissiveColor = mat.diffuseColor.scale(0.05);
      master = BABYLON.MeshBuilder.CreateCylinder('bond_'+key,{height:1,diameter:1,tessellation:22},scene);
      master.isPickable=true; master.thinInstanceEnablePicking=true;
      registerMaster('bondMaster', master);
    }
    if (master && mat) master.material = mat;
    const g = { master, mats:[], indices:[] };
    store.set(key,g);
    // (debug/backfill removed)
    return g;
  }
  const ensureAtomGroup = (el)=>{ __count('moleculeView#ensureAtomGroup'); return ensureGroup(atomGroups, el, 'atom'); };
  const ensureGhostAtomGroup = (el)=>{ __count('moleculeView#ensureGhostAtomGroup'); return ensureGroup(ghostAtomGroups, el, 'ghostAtom'); };
  const ensureGhostBondGroup = (key)=>{ __count('moleculeView#ensureGhostBondGroup'); return ensureGroup(ghostBondGroups, key, 'ghostBond'); };
  const ensureBondGroup = (key)=>{ __count('moleculeView#ensureBondGroup'); return ensureGroup(bondGroups, key, 'bond'); };
  function ensureForceGroup(){
    __count('moleculeView#ensureForceGroup');
    if (forceGroups.has('force')) return forceGroups.get('force');
    const mat = new BABYLON.StandardMaterial('force_vec_mat', scene);
    mat.diffuseColor = new BABYLON.Color3(0.95,0.20,0.20);
    mat.emissiveColor = mat.diffuseColor.scale(0.55);
    mat.specularColor = new BABYLON.Color3(0.1,0.1,0.1);
    mat.disableLighting = true; // keep vivid in dim VR/AR lighting
    // Use unit diameter & scale per-instance so radius logic matches bonds; previous tiny base made arrows microscopic
    const master = BABYLON.MeshBuilder.CreateCylinder('force_vector_master',{height:1,diameter:1,tessellation:14},scene);
    master.isPickable = false; master.thinInstanceEnablePicking = false;
    master.material = mat;
    // Hide until we have instances
    if (typeof master.setEnabled === 'function') master.setEnabled(false); else master.isVisible=false;
    const g = { master, mats:[], indices:[], colors:[] };
    forceGroups.set('force', g);
    registerMaster('forceMaster', master);
    return g;
  }
  function buildInitial() {
    __count('moleculeView#buildInitial');
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
    __count('moleculeView#rebuildBonds');
    for (const g of bondGroups.values()) { g.mats=[]; g.indices=[]; }
    const source = bondData || molState.bonds;
    // Debug logging removed (previously controlled by debug=1 query param) to reduce console noise.
  const ODBG_ENABLED = (typeof window === 'undefined') ? false : (window.O_BOND_DEBUG === true); // default OFF
    let oBondCount = 0;
    for (const b of source) {
      const g = ensureBondGroup(keyOf(b.i,b.j));
      const pA = molState.positions[b.i]; const pB = molState.positions[b.j];
      const mat = bondMatrix(pA,pB,0.1);
      g.mats.push(mat); g.indices.push(b);
      if (ODBG_ENABLED) {
        try {
          const elI = molState.elements[b.i]; const elJ = molState.elements[b.j];
          const involvesO = (elI === 'O' || elJ === 'O');
          if (involvesO) {
            oBondCount++;
            // Extract rotation & mid info from matrix for quick diagnostics
            const midX = (pA.x + pB.x)/2, midY = (pA.y + pB.y)/2, midZ = (pA.z + pB.z)/2;
            console.log(`[O-BOND-DBG][rebuildBonds] bond ${elI}-${elJ} i=${b.i} j=${b.j} key=${keyOf(b.i,b.j)} mid=(${midX.toFixed(3)},${midY.toFixed(3)},${midZ.toFixed(3)}) len=${Math.hypot(pB.x-pA.x,pB.y-pA.y,pB.z-pA.z).toFixed(3)} opacity=${b.opacity!=null?b.opacity:1}`);
            // Additional request: for O-C bonds specifically, log local & world positions of both atoms.
            if ((elI === 'O' && elJ === 'C') || (elI === 'C' && elJ === 'O')) {
              try {
                // Local positions are just pA / pB (molecule space)
                const lpA = pA, lpB = pB;
                let wpA = lpA, wpB = lpB;
                // If atom masters exist & have world matrices (after scene graph build), transform local to world.
                // We find the atom master for each element; thin instance positions are encoded in the matrix buffer;
                // world position of a given thin instance = Transform(localPosition, masterWorldMatrix), but because we
                // parent masters to moleculeRoot (and don't apply per-instance extra transforms beyond translation), we can approximate
                // by applying moleculeRoot's world matrix if available.
                const root = molState.moleculeRoot;
                if (root && root.getWorldMatrix) {
                  const wm = root.getWorldMatrix();
                  const vA = BABYLON.Vector3.TransformCoordinates(new BABYLON.Vector3(lpA.x,lpA.y,lpA.z), wm);
                  const vB = BABYLON.Vector3.TransformCoordinates(new BABYLON.Vector3(lpB.x,lpB.y,lpB.z), wm);
                  wpA = { x: vA.x, y: vA.y, z: vA.z };
                  wpB = { x: vB.x, y: vB.y, z: vB.z };
                } else if (root && root._worldMatrix) { // fallback if using stub storing _worldMatrix
                  try {
                    const wm = root._worldMatrix;
                    const vA = BABYLON.Vector3.TransformCoordinates(new BABYLON.Vector3(lpA.x,lpA.y,lpA.z), wm);
                    const vB = BABYLON.Vector3.TransformCoordinates(new BABYLON.Vector3(lpB.x,lpB.y,lpB.z), wm);
                    wpA = { x: vA.x, y: vA.y, z: vA.z };
                    wpB = { x: vB.x, y: vB.y, z: vB.z };
                  } catch {}
                }
                console.log(
                  `[O-BOND-DBG][OC-pos] ${elI}-${elJ} i=${b.i} j=${b.j} ` +
                  `localA=(${lpA.x.toFixed(3)},${lpA.y.toFixed(3)},${lpA.z.toFixed(3)}) ` +
                  `localB=(${lpB.x.toFixed(3)},${lpB.y.toFixed(3)},${lpB.z.toFixed(3)}) ` +
                  `worldA=(${wpA.x.toFixed(3)},${wpA.y.toFixed(3)},${wpA.z.toFixed(3)}) ` +
                  `worldB=(${wpB.x.toFixed(3)},${wpB.y.toFixed(3)},${wpB.z.toFixed(3)})`
                );
              } catch (e) {
                console.log('[O-BOND-DBG][OC-pos][error]', e);
              }
            }
          }
        } catch {}
      }
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
      // (debug log suppressed)
    }
    if (ODBG_ENABLED && oBondCount>0) {
      try {
        // After rebuild, verify parenting & log rotationQuaternion of moleculeRoot & each O-* bond master
        const rootQ = molState.moleculeRoot && molState.moleculeRoot.rotationQuaternion;
        const rootRotStr = rootQ ? `rootQ=(${rootQ.x.toFixed(4)},${rootQ.y.toFixed(4)},${rootQ.z.toFixed(4)},${rootQ.w.toFixed(4)})` : 'rootQ=none';
        for (const [key,g] of bondGroups) {
          if (/O-/.test(key) || /-O/.test(key)) {
            const parentOk = g.master.parent === molState.moleculeRoot;
            const q = g.master.rotationQuaternion;
            const qStr = q ? `q=(${q.x.toFixed(4)},${q.y.toFixed(4)},${q.z.toFixed(4)},${q.w.toFixed(4)})` : 'q=none';
            // Also log Euler rotation if quaternion absent, and world matrix first row as a quick change signature
            let eulerStr = '';
            try {
              if (!q && g.master.rotation) {
                const r = g.master.rotation; eulerStr = ` rotEuler=(${(r.x||0).toFixed(3)},${(r.y||0).toFixed(3)},${(r.z||0).toFixed(3)})`;
              }
            } catch {}
            let wmSig='';
            try {
              const wm = g.master.getWorldMatrix && g.master.getWorldMatrix();
              if (wm && wm.m) wmSig = ` wm0=(${wm.m[0].toFixed(3)},${wm.m[1].toFixed(3)},${wm.m[2].toFixed(3)})`;
            } catch {}
            console.log(`[O-BOND-DBG][postRebuild] key=${key} instances=${g.mats.length} parentOk=${parentOk} ${qStr} ${rootRotStr}${eulerStr}${wmSig}`);
          }
        }
      } catch {}
    }
    // After bonds are rebuilt, rebuild forces to stay in sync with transforms (shared masters for VR)
    try { rebuildForces(); } catch {}
  }
  // Force vector thin instance matrix builder (arrow style: cylinder scaled & rotated)
  function forceMatrix(p, f, length, radius){
    // p: {x,y,z}; f: [fx,fy,fz]; length: fixed length for debug visualization
    const fx=f[0], fy=f[1], fz=f[2];
    const mag = Math.hypot(fx,fy,fz) || 1e-9;
    const dir = new BABYLON.Vector3(fx/mag, fy/mag, fz/mag);
    const up = new BABYLON.Vector3(0,1,0);
    let rot; const dot = BABYLON.Vector3.Dot(up, dir);
    if (dot>0.9999) rot = BABYLON.Quaternion.Identity();
    else if (dot<-0.9999) rot = BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1,0,0), Math.PI);
    else { const axis = BABYLON.Vector3.Cross(up, dir).normalize(); rot = BABYLON.Quaternion.RotationAxis(axis, Math.acos(Math.min(1,Math.max(-1,dot)))); }
    // Position: offset so base sits at atom, cylinder centered => shift half-length along dir
    const mid = new BABYLON.Vector3(p.x + dir.x*length/2, p.y + dir.y*length/2, p.z + dir.z*length/2);
    return BABYLON.Matrix.Compose(new BABYLON.Vector3(radius*2, length, radius*2), rot, mid);
  }
  function rebuildForces(){
    __count('moleculeView#rebuildForces');
    const DBG = (typeof window !== 'undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search||''));
    if (DBG) console.log('[Forces][rebuild] start');
    const g = ensureForceGroup();
    g.mats.length = 0; g.indices.length = 0;
    // Accept forces from molState.forces OR molState.dynamics?.forces OR global window.__RELAX_FORCES for backward compat
    let forces = molState.forces || (molState.dynamics && molState.dynamics.forces && molState.dynamics.forces.map(f=>[f.x,f.y,f.z])) || (typeof window!=='undefined' && window.__RELAX_FORCES) || [];
    if (!Array.isArray(forces)) forces = [];
    const n = Math.min(forces.length, molState.positions.length);
    if (!n) {
      if (DBG) console.log('[Forces][rebuild] no forces or positions (n=0)');
      g.master.thinInstanceSetBuffer('matrix', new Float32Array());
      if (typeof g.master.setEnabled === 'function') g.master.setEnabled(false); else g.master.isVisible=false;
      return;
    }
    // Scaling controls
  const hasWin = (typeof window !== 'undefined');
  const search = hasWin ? (window.location?.search||'') : '';
  const qs = Object.fromEntries(Array.from(new URLSearchParams(search)).map(([k,v])=>[k,v]));
  const forcedFixed = ('forceFixed' in qs) || (hasWin && window.FORCE_FIXED); // keep constant length if set
  const forceScale = parseFloat(qs.forceScale || (hasWin && window.FORCE_SCALE) || '0.4'); // length per |f| unit
  const maxLen = parseFloat(qs.forceMax || (hasWin && window.FORCE_MAX) || '1.2');
  const minLen = parseFloat(qs.forceMin || (hasWin && window.FORCE_MIN) || '0.12');
  const radius = parseFloat(qs.forceRadius || (hasWin && window.FORCE_RADIUS) || '0.05');
    let maxMag = 0; if (!forcedFixed) { for (let i=0;i<n;i++){ const f=forces[i]; if(!f) continue; const m=Math.hypot(f[0],f[1],f[2]); if(m>maxMag) maxMag=m; } }
    if (DBG) console.log('[Forces][params]', { forcedFixed, forceScale, maxLen, minLen, radius, maxMag: maxMag.toFixed ? maxMag.toFixed(4):maxMag });
    const fixedLen = forcedFixed ? (parseFloat(qs.forceFixed) || 0.9) : null;
    for (let i=0;i<n;i++){
      const p = molState.positions[i]; const f = forces[i];
      if(!p || !f) continue;
      const mag = Math.hypot(f[0],f[1],f[2]);
      if(mag < 1e-10) continue;
      let drawLen;
      if (fixedLen != null) drawLen = fixedLen; else drawLen = Math.min(maxLen, Math.max(minLen, mag * forceScale));
      const mat = forceMatrix(p, f, drawLen, radius);
      g.mats.push(mat); g.indices.push({ atom:i, mag });
      if (DBG && i < 8) console.log('[Forces][rebuild] atom', i, 'f=', f, 'mag=', mag.toFixed(4));
    }
    g.master.thinInstanceSetBuffer('matrix', flattenMatrices(g.mats));
    try { g.master.thinInstanceRefreshBoundingInfo && g.master.thinInstanceRefreshBoundingInfo(); } catch {}
    // Per-instance color (solid red, alpha 1)
    if (g.mats.length) {
      const cols = new Float32Array(g.mats.length * 4);
      for (let i=0;i<g.mats.length;i++) { cols[i*4+0]=0.95; cols[i*4+1]=0.05; cols[i*4+2]=0.05; cols[i*4+3]=1.0; }
      try { g.master.thinInstanceSetBuffer('color', cols, 4); } catch {}
    } else {
      try { g.master.thinInstanceSetBuffer('color', new Float32Array(), 4); } catch {}
    }
    if (DBG) console.log('[Forces][rebuild] instances=', g.mats.length);
    if (g.mats.length) {
      if (typeof g.master.setEnabled === 'function') g.master.setEnabled(true); else g.master.isVisible=true;
    } else {
      if (typeof g.master.setEnabled === 'function') g.master.setEnabled(false); else g.master.isVisible=false;
    }
    if (DBG) console.log('[Forces][rebuild] complete visible=', g.mats.length>0);
  }
  function rebuildGhosts() {
    __count('moleculeView#rebuildGhosts');
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
    // (ghost bond debug removed)
  }
  function rebuildCellLines() {
    __count('moleculeView#rebuildCellLines');
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
    __count('moleculeView#rebuildAtoms');
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
    // (debug log suppressed)
  }
  function updatePositions() {
    __count('moleculeView#updatePositions');
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
  try { rebuildForces(); } catch {}
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
    // (debug logs suppressed)
  });
  molState.bus.on('cellChanged', () => { rebuildCellLines(); rebuildGhosts(); });
  // Forces update event (consumer code can emit this when new forces computed)
  if (molState.bus && typeof molState.bus.on === 'function') {
    molState.bus.on('forcesChanged', () => {
      const DBG = (typeof window !== 'undefined') && (window.FORCE_DEBUG || /[?&]forceDebug=1/.test(window.location?.search||''));
      if (DBG) console.log('[Forces][event] forcesChanged received');
      try { rebuildForces(); } catch(e){ if (DBG) console.warn('[Forces][event] rebuild error', e); }
    });
  }
  // Selection highlight management
  function ensureHighlightMeshes() {
    __count('moleculeView#ensureHighlightMeshes');
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
    __count('moleculeView#updateSelectionHighlight');
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
      // VR Regression Fix: Previously bond highlight used world coordinates once at selection time
      // and was never updated when masters rotated (after removing VR's own per-frame marker code).
      // We now parent the bond highlight to the same master mesh that hosts the bond thin instances.
      // This lets it inherit rotation/scale automatically just like the atoms/bonds.
      let bondMaster = null;
      for (const [_, g] of bondGroups) {
        if (g && g.indices && g.indices.some(b => (b.i===i && b.j===j) || (b.i===j && b.j===i))) { bondMaster = g.master; break; }
      }
      if (bondMaster && highlight.bond.parent !== bondMaster) {
        highlight.bond.parent = bondMaster;
      }
      if (!bondMaster && highlight.bond.parent) {
        // Fallback: detach if master vanished (e.g., molecule switch)
        highlight.bond.parent = null;
      }
      // Compute local midpoint and vector using molecule-space coordinates (positions[] already local).
      const midLocal = new BABYLON.Vector3((pA.x+pB.x)/2,(pA.y+pB.y)/2,(pA.z+pB.z)/2);
      const vLocal = new BABYLON.Vector3(pB.x-pA.x,pB.y-pA.y,pB.z-pA.z);
      const lenLocal = vLocal.length();
      const up = new BABYLON.Vector3(0,1,0);
      let rotQ; const d=vLocal.normalizeToNew(); const dot=BABYLON.Vector3.Dot(up,d);
      if (dot>0.9999) rotQ=BABYLON.Quaternion.Identity();
      else if (dot<-0.9999) rotQ=BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1,0,0),Math.PI);
      else { const axis=BABYLON.Vector3.Cross(up,d).normalize(); rotQ=BABYLON.Quaternion.RotationAxis(axis,Math.acos(dot)); }
      // Apply local transform; master rotation/scale (if any) will be inherited.
      highlight.bond.position = midLocal;
      highlight.bond.rotationQuaternion = rotQ;
      const radius = 0.16; // local radius shell; master scaling will modify in world space
      highlight.bond.scaling = new BABYLON.Vector3(radius*2, lenLocal, radius*2);
      highlight.bond.isVisible = true;
      // O-* bond highlight debug
      try {
  const ODBG = (typeof window === 'undefined') ? false : (window.O_BOND_DEBUG === true);
        if (ODBG) {
          const elI = molState.elements[i]; const elJ = molState.elements[j];
            if (elI==='O' || elJ==='O') {
              const bm = bondMaster && bondMaster.name;
              console.log(`[O-BOND-DBG][highlight] sel ${elI}-${elJ} i=${i} j=${j} master=${bm} localLen=${lenLocal.toFixed(3)}`);
            }
        }
      } catch {}
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
    // If a bond remains selected, recompute its local transform each frame so any
    // underlying atom position mutations (e.g., dragging or dynamics) reflect immediately.
    if (sel.kind === 'bond' && highlight.bond.isVisible) {
      try {
        const { i, j } = sel.data;
        const pA = molState.positions[i];
        const pB = molState.positions[j];
        if (pA && pB) {
          const midLocal = new BABYLON.Vector3((pA.x+pB.x)/2,(pA.y+pB.y)/2,(pA.z+pB.z)/2);
          const vLocal = new BABYLON.Vector3(pB.x-pA.x,pB.y-pA.y,pB.z-pA.z);
          const lenLocal = vLocal.length();
          if (lenLocal > 1e-6) {
            const up = new BABYLON.Vector3(0,1,0);
            let rotQ; const d=vLocal.normalizeToNew(); const dot=BABYLON.Vector3.Dot(up,d);
            if (dot>0.9999) rotQ=BABYLON.Quaternion.Identity();
            else if (dot<-0.9999) rotQ=BABYLON.Quaternion.RotationAxis(new BABYLON.Vector3(1,0,0),Math.PI);
            else { const axis=BABYLON.Vector3.Cross(up,d).normalize(); rotQ=BABYLON.Quaternion.RotationAxis(axis,Math.acos(dot)); }
            highlight.bond.position = midLocal;
            highlight.bond.rotationQuaternion = rotQ;
            // Keep scaling local; world scaling inherited
            highlight.bond.scaling.y = lenLocal; // preserve radius components
            // (per-frame O-bond debug removed)
          }
        }
      } catch {}
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

    // (advanced rotation logging removed)
  });
  buildInitial();
  // Initialize highlight once after initial build
  ensureHighlightMeshes();
  function resolveAtomPick(pick) {
    __count('moleculeView#resolveAtomPick');
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
    __count('moleculeView#resolveBondPick');
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
  return { rebuildBonds, rebuildGhosts, rebuildForces, _internals:{ atomGroups, bondGroups, forceGroups, ghostAtomGroups, ghostBondGroups, highlight }, resolveAtomPick, resolveBondPick };
}
