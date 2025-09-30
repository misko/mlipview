// Core spherical drag radial delta computations (framework agnostic for testing)
// Exports computeRadialDelta(state, inputs, config)
// state: { initialRadius, controllerDist0, ctrlVec0:[x,y,z], initialCtrlForward:[x,y,z] }
// inputs: { controllerPos:[x,y,z], controllerForward:[x,y,z], cameraPos:[x,y,z] }
// config: { mode, forwardGain, adaptiveGain, adaptiveMaxFrac, expK, expGain }
// Returns { deltaR, components:{ distance, projection, forward, adaptive, exp }, ratio }

function vecSub(a,b){ return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
function vecDot(a,b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
function vecLen(a){ return Math.hypot(a[0],a[1],a[2]); }
function vecNorm(a){ const l=vecLen(a)||1; return [a[0]/l,a[1]/l,a[2]/l]; }

export function computeRadialDelta(state, inputs, config={}){
  const mode = config.mode || 'distance';
  const fGain = config.forwardGain!=null?config.forwardGain:2.5;
  const adaptGain = config.adaptiveGain!=null?config.adaptiveGain:0.15;
  const adaptMax = config.adaptiveMaxFrac!=null?config.adaptiveMaxFrac:0.75;
  const expK = config.expK!=null?config.expK:1.2;
  const expGain = config.expGain!=null?config.expGain:0.4;
  const center = inputs.cameraPos;
  const ctrlPos = inputs.controllerPos;
  const ctrlForward = inputs.controllerForward;
  const ctrlVec = vecSub(ctrlPos, center);
  const distCtrl = vecLen(ctrlVec);
  const distance = distCtrl - state.controllerDist0;
  let projection = 0; try { const baseDir = vecNorm(state.ctrlVec0); projection = vecDot(ctrlVec, baseDir) - vecDot(state.ctrlVec0, baseDir); } catch {}
  let forward = 0; try { const camForward = vecNorm(ctrlForward); forward = vecDot(ctrlVec, camForward) - vecDot(state.ctrlVec0, camForward); } catch {}
  let adaptive = 0; let ratio = 1;
  if(mode==='adaptive'){ const base=state.controllerDist0||1e-4; ratio = distCtrl/base; let frac=(ratio-1)*adaptGain; if(frac>adaptMax) frac=adaptMax; else if(frac<-adaptMax) frac=-adaptMax; adaptive = frac*state.initialRadius; }
  let exp=0; if(mode==='exp'){ const base=state.controllerDist0||1e-4; ratio = distCtrl/base; const f = Math.pow(ratio, expK)-1; exp = f*expGain*state.initialRadius; }
  let deltaR;
  if(mode==='projection') deltaR=projection;
  else if(mode==='hybrid') deltaR=(distance+projection)*0.5;
  else if(mode==='adaptive') deltaR=adaptive;
  else if(mode==='forward') deltaR=forward*fGain;
  else if(mode==='exp') deltaR=exp;
  else deltaR=distance;
  return { deltaR, components:{ distance, projection, forward, adaptive, exp }, ratio };
}

export default { computeRadialDelta };
