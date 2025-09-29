import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createManipulationService } from '../public/domain/manipulationService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import fs from 'fs';
import path from 'path';

// Babylon light stubs (no rendering needed for geometry math)
if (!global.BABYLON) global.BABYLON = {};
BABYLON.Color3 = function(r,g,b){ this.r=r; this.g=g; this.b=b; this.clone=()=>new BABYLON.Color3(r,g,b); this.scale=f=>new BABYLON.Color3(r*f,g*f,b*f); };
BABYLON.StandardMaterial = function(){ this.diffuseColor=new BABYLON.Color3(1,1,1); this.emissiveColor=new BABYLON.Color3(0,0,0); };
BABYLON.Vector3 = function(x,y,z){ this.x=x; this.y=y; this.z=z; this.length=()=>Math.hypot(x,y,z); };
BABYLON.Vector3.Dot=(a,b)=>a.x*b.x+a.y*b.y+a.z*b.z; BABYLON.Vector3.Cross=(a,b)=>new BABYLON.Vector3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
BABYLON.Vector3.prototype.normalizeToNew=function(){ const L=this.length()||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L); };
BABYLON.Vector3.prototype.normalize=function(){ const L=this.length()||1; this.x/=L; this.y/=L; this.z/=L; return this; };
BABYLON.Quaternion = { Identity:()=>({}), RotationAxis:()=>({}) };
BABYLON.Matrix = class { constructor(){ this.m=new Float32Array(16);} static Compose(){ return new BABYLON.Matrix(); } };
BABYLON.MeshBuilder = { CreateSphere:()=>({ thinInstanceSetBuffer(){}, material:null }), CreateCylinder:()=>({ thinInstanceSetBuffer(){}, material:null }), CreateLines:()=>({ dispose(){}, color:null }) };

function distance(a,b){ const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return Math.hypot(dx,dy,dz); }

// Reasonable length thresholds (Angstrom): C-C in benzene ~1.39; C-H ~1.08. Allow tolerance + distortion.
// Allow generous deformation: some rotations may transiently reduce projection distance or reorder bonds.
const CC_MIN = 0.9; const CC_MAX = 2.0; // anything longer is suspicious for a single bond in ring after rotation
const CH_MIN = 0.7; const CH_MAX = 1.5;

// Helper: find index pairs for current bonds classified by element pair.
function classifyBonds(state) {
  const cc = []; const ch = [];
  for (const b of state.bonds) {
    const e1 = state.elements[b.i]; const e2 = state.elements[b.j];
    const len = distance(state.positions[b.i], state.positions[b.j]);
    if ((e1==='C' && e2==='C')) cc.push({ ...b, len });
    else if ((e1==='C' && e2==='H') || (e1==='H' && e2==='C')) ch.push({ ...b, len });
  }
  return { cc, ch };
}

// We enable cell + ghosts to ensure periodic/ghost logic does not create spurious long bonds.

// __dirname replacement not required; use process.cwd() rooted paths.

describe('Benzene bond rotation keeps reasonable bond lengths (with cell + ghosts)', () => {
  test('rotating a selected C-C bond by various angles maintains sane lengths', () => {
  const benzenePath = path.join(process.cwd(), 'public', 'molecules', 'benzene.xyz');
    const xyz = fs.readFileSync(benzenePath, 'utf8');
    const parsed = parseXYZ(xyz);
    const state = createMoleculeState();
    applyXYZToState(state, parsed);
    // Synthesize a cell box around benzene (rough extents) and enable ghosts
    state.toggleCellVisibilityEnhanced?.(); // creates bounding cell
    state.toggleGhostCells();
    const bondService = createBondService(state);
    bondService.recomputeAndStore();
  const selection = createSelectionService(state);
  const manipulation = createManipulationService(state, { bondService });
  // Create a minimal scene/view stub to allow ghost bond generation visibility
  const scene = { onPointerObservable:{ add(){} } };
  const view = createMoleculeView(scene, state);
  // Ensure cell+ghosts ON exactly once (avoid double toggle turning them off)
  if (!state.showCell) state.toggleCellVisibilityEnhanced?.();
  if (!state.showGhostCells) state.toggleGhostCells();
  state.markCellChanged();

  // Simulate selecting a specific C-C bond: between atom 0 (C) and 1 (C)
  // Legacy selection clickBond expects { i, j }
  selection.clickBond({ i:0, j:1 });
  // Orientation data is normally populated by adapter; patch minimal orientation so rotation side logic works.
  state.selection.data.orientation = { from:0, to:1 };

    const anglesDeg = [5,10,15,25,50];
    for (const deg of anglesDeg) {
      const rad = deg * Math.PI / 180;
      manipulation.rotateBond(rad, 'j'); // rotate side anchored at atom 1 (arbitrary consistent choice)
      // After rotation, recompute bonds has been invoked by manipulation service.
  // Force ghost rebuild (in case rotation changed positions affecting augmented bonding)
  view.rebuildGhosts();
      const { cc, ch } = classifyBonds(state);
      // Assert ghost bond groups have at least one entry (cross-image bonds appear) after first rotation
      const ghostGroupBondCount = Array.from(view._internals.ghostBondGroups.values()).reduce((s,g)=>s+g.mats.length,0);
      expect(ghostGroupBondCount).toBeGreaterThan(0);
      expect(cc.length).toBeGreaterThanOrEqual(5); // ring may still have at least 5 C-C bonds resolvable
      for (const b of cc) {
        if (!(b.len >= CC_MIN && b.len <= CC_MAX)) {
          throw new Error(`C-C bond length out of range after ${deg}deg: ${b.i}-${b.j} len=${b.len.toFixed(3)}`);
        }
      }
      expect(ch.length).toBeGreaterThanOrEqual(6);
      for (const b of ch) {
        if (!(b.len >= CH_MIN && b.len <= CH_MAX)) {
          throw new Error(`C-H bond length out of range after ${deg}deg: ${b.i}-${b.j} len=${b.len.toFixed(3)}`);
        }
      }
    }
  });
});
