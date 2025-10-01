/** @jest-environment jsdom */
// Continuous relax run test: uses startRelaxContinuous pacing/backoff and compares
// final energy to reference aggregate /relax multi-step result (ASE backend parity).
import http from 'http';
import https from 'https';

if (typeof fetch === 'undefined') {
  global.fetch = function(nodeUrl, opts={}){
    return new Promise((resolve, reject) => {
      try {
        const url = new URL(nodeUrl);
        const lib = url.protocol === 'https:' ? https : http;
        const req = lib.request(url, { method: opts.method||'GET', headers: opts.headers||{} }, res => {
          const chunks = [];
            res.on('data', d=> chunks.push(d));
            res.on('end', ()=> {
              const body = Buffer.concat(chunks).toString('utf8');
              resolve({ ok: res.statusCode>=200 && res.statusCode<300, status: res.statusCode, json: async ()=> JSON.parse(body), text: async ()=> body });
            });
        });
        req.on('error', reject);
        if (opts.body) req.write(opts.body);
        req.end();
      } catch(e) { reject(e); }
    });
  };
}

function haveServer(url='http://127.0.0.1:8000'){ return new Promise(resolve=>{ const body=JSON.stringify({atomic_numbers:[1],coordinates:[[0,0,0]],properties:['energy'],calculator:'lj'}); const req=http.request(url+'/simple_calculate',{method:'POST',headers:{'Content-Type':'application/json'}},res=>resolve(res.statusCode>=200&&res.statusCode<300)); req.on('error',()=>resolve(false)); req.setTimeout(500,()=>{ try{req.destroy();}catch{} resolve(false);}); req.end(body); }); }

async function callRelax(calculator, steps){ const atomic_numbers=[8,1,1]; const coordinates=[[0,0,0],[0.9575,0,0],[-0.2399872,0.92662721,0]]; const body={ atomic_numbers, coordinates, steps, calculator }; const resp = await fetch('http://127.0.0.1:8000/relax',{ method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body) }); if(!resp.ok) throw new Error('relax failed '+resp.status); return await resp.json(); }

beforeAll(()=>{
  if(!global.BABYLON){
    global.BABYLON = { TransformNode: class {}, MeshBuilder: { CreateCylinder: ()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial: class { constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }}, Color3: class {}, Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new global.BABYLON.Vector3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);} normalize(){ const m=Math.hypot(this.x,this.y,this.z)||1; this.x/=m; this.y/=m; this.z/=m; return this;} rotateByQuaternionToRef(_, out){ out.x=this.x; out.y=this.y; out.z=this.z; return out; } }, Quaternion: class { static Identity(){ return {}; } static RotationAxis(){ return {}; } }, Scene: class {} };
  }
});

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn)=>{} },
    scene: { meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } },
    camera: { attachControl: ()=>{} }
  })
}));

async function setupViewer(){
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  window.__MLIP_FEATURES = { RELAX_LOOP:true, MD_LOOP:false, ENERGY_TRACE:true, FORCE_VECTORS:false };
  const canvas=document.createElement('canvas'); canvas.id='viewer'; document.body.appendChild(canvas);
  const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
  const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=300; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
  const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);
  const mod = await import('../public/index.js');
  return await mod.initNewViewer(canvas, { elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.9575,y:0,z:0},{x:-0.2399872,y:0.92662721,z:0}], bonds:[] });
}

describe('continuous relax run parity', () => {
  test('final energy within tolerance of aggregate relax', async () => {
    const up = await haveServer(); if(!up){ console.warn('[relax run parity] server not reachable; skipping'); return; }
    const steps=50; // use fewer steps for test time
    const ref = await callRelax('uma', steps);
    const api = await setupViewer();
    const res = await api.startRelaxContinuous({ maxSteps: steps });
    expect(res.disabled).not.toBe(true);
    // Expect at least steps executed
    expect(res.steps).toBeGreaterThanOrEqual(steps*0.9); // allow early abort if converged
    const finalBrowser = (window.__RELAX_TRACE||[]).slice(-1)[0];
    const tol = 5e-3;
    expect(Math.abs(ref.final_energy - finalBrowser)).toBeLessThan(tol);
  }, 60000);
});
