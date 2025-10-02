/** @jest-environment jsdom */
// Playwright-style parity test simulated under jsdom.
// Compares incremental 1-step relaxStep calls vs aggregate /relax 20-step result.
import http from 'http';
import { haveServer } from './helpers/server.js';
import https from 'https';

// Minimal fetch polyfill for Node <18 or test env without global fetch
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


async function callRelax(calculator){ const atomic_numbers=[8,1,1]; const coordinates=[[0,0,0],[0.9575,0,0],[-0.2399872,0.92662721,0]]; const body={ atomic_numbers, coordinates, steps:20, calculator }; const resp = await fetch('http://127.0.0.1:8000/serve/relax',{ method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body) }); if(!resp.ok) throw new Error('relax failed '+resp.status); return await resp.json(); }

// Minimal BABYLON mocks required by index.js
beforeAll(()=>{
  if(!global.BABYLON){
    global.BABYLON = {
      TransformNode: class {},
      MeshBuilder: { CreateCylinder: ()=>({ parent:null, dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) },
      StandardMaterial: class { constructor(){ this.diffuseColor={}; this.emissiveColor={}; this.specularColor={}; }},
      Color3: class { constructor(){ } },
      Vector3: class { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } static Up(){ return new global.BABYLON.Vector3(0,1,0);} static Dot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z;} static Cross(a,b){ return new global.BABYLON.Vector3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);} normalize(){ const m=Math.hypot(this.x,this.y,this.z)||1; this.x/=m; this.y/=m; this.z/=m; return this;} rotateByQuaternionToRef(_, out){ out.x=this.x; out.y=this.y; out.z=this.z; return out; } },
      Quaternion: class { static Identity(){ return {}; } static RotationAxis(){ return {}; } },
      Scene: class {}
    };
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
  // Ensure base URL configured before viewer init
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  const canvas=document.createElement('canvas'); canvas.id='viewer'; document.body.appendChild(canvas);
  const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
  const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=300; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
  const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);
  const mod = await import('../public/index.js');
  return await mod.initNewViewer(canvas, { elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.9575,y:0,z:0},{x:-0.2399872,y:0.92662721,z:0}], bonds:[] });
}

async function performBrowserSteps(api, n=20){ for(let i=0;i<n;i++){ await api.relaxStep(); } return window.__RELAX_TRACE || []; }

['uma'].forEach(calc => { // relaxStep currently hardcodes UMA backend calculator
  describe(`browser 20-step parity ${calc}`, () => {
    test('first/last energy matches server /relax aggregate call within tolerance', async () => {
      const up = await haveServer();
      if(!up){ console.warn('[browser parity] server not reachable; skipping'); return; }
      // Reference energies from single aggregate /relax 20-step call
      const ref = await callRelax(calc);
      const api = await setupViewer();
      const trace = await performBrowserSteps(api, 20);
      expect(trace.length).toBeGreaterThanOrEqual(1);
      // Fetch initial single-point energy via /simple_calculate for the calculator in use
      const initBody = { atomic_numbers:[8,1,1], coordinates:[[0,0,0],[0.9575,0,0],[-0.2399872,0.92662721,0]], properties:['energy'], calculator:calc };
  const initResp = await fetch('http://127.0.0.1:8000/serve/simple', { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(initBody) });
      const initJson = await initResp.json();
      const initEnergy = initJson.results.energy;
      const finalBrowser = trace[trace.length-1];
      const tolInit = calc==='lj'? 2e-2 : 5e-3;
      const tolFinal = calc==='lj'? 1e-2 : 5e-3;
      expect(Math.abs(ref.initial_energy - initEnergy)).toBeLessThan(tolInit);
      expect(Math.abs(ref.final_energy - finalBrowser)).toBeLessThan(tolFinal);
    }, 30000);
  });
});
