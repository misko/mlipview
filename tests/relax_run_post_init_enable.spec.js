/** @jest-environment jsdom */
// Verifies enabling RELAX_LOOP after viewer init now works (dynamic flags).
import https from 'https';
import { haveServer } from './helpers/server.js';
if (typeof fetch === 'undefined') {
  global.fetch = function(nodeUrl, opts={}){ return new Promise((resolve,reject)=>{ try{ const url=new URL(nodeUrl); const lib=url.protocol==='https:'?https:http; const req=lib.request(url,{method:opts.method||'GET',headers:opts.headers||{}},res=>{ const chunks=[];res.on('data',d=>chunks.push(d));res.on('end',()=>{ const body=Buffer.concat(chunks).toString('utf8'); resolve({ ok:res.statusCode>=200&&res.statusCode<300, status:res.statusCode, json:async()=>JSON.parse(body), text:async()=>body }); });}); req.on('error',reject); if(opts.body) req.write(opts.body); req.end(); } catch(e){ reject(e); } });};
}
// haveServer imported from helpers/server.js (uses /serve/health)

beforeAll(()=>{ if(!global.BABYLON){ global.BABYLON={ TransformNode:class{}, MeshBuilder:{ CreateCylinder:()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial:class{constructor(){this.diffuseColor={};this.emissiveColor={};this.specularColor={};}}, Color3:class{}, Vector3:class{ constructor(x=0,y=0,z=0){this.x=x;this.y=y;this.z=z;} static Up(){return new global.BABYLON.Vector3(0,1,0);} }, Quaternion:class{ static Identity(){return{};}}, Scene:class{} }; }});

jest.mock('../public/render/scene.js', ()=>({ createScene: async ()=>({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(f){this._l.push(f);} } }, camera:{ attachControl:()=>{} } }) }));

async function initViewer(){ window.__MLIPVIEW_SERVER='http://127.0.0.1:8000'; const canvas=document.createElement('canvas'); canvas.id='viewer'; document.body.appendChild(canvas); const energyDiv=document.createElement('div'); energyDiv.id='energyPlot'; document.body.appendChild(energyDiv); const c=document.createElement('canvas'); c.id='energyCanvas'; c.width=200; c.height=60; c.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}}); energyDiv.appendChild(c); const label=document.createElement('div'); label.id='energyLabel'; energyDiv.appendChild(label); const mod = await import('../public/index.js'); return await mod.initNewViewer(canvas,{ elements:[{Z:8},{Z:1},{Z:1}], positions:[{x:0,y:0,z:0},{x:0.95,y:0,z:0},{x:-0.24,y:0.93,z:0}], bonds:[] }); }

describe('post-init feature flag enable', ()=>{
  test('relax run executes after enabling flag dynamically', async ()=>{
    if(!(await haveServer())){ console.warn('[skip] server not up'); return; }
    const api = await initViewer();
  // With auto-enabled defaults, first run should already work
  const first = await api.startRelaxContinuous({ maxSteps:5 });
  expect(first.disabled).not.toBe(true);
  expect(first.steps).toBeGreaterThan(0);
  // Disable then re-enable to prove dynamic control still works
  api.enableFeatureFlag('RELAX_LOOP', false);
  const disabledAttempt = await api.startRelaxContinuous({ maxSteps:3 });
  expect(disabledAttempt.disabled).toBe(true);
  api.enableFeatureFlag('RELAX_LOOP', true);
  const second = await api.startRelaxContinuous({ maxSteps:3 });
  expect(second.disabled).not.toBe(true);
  expect(second.steps).toBeGreaterThan(0);
  }, 20000);
});
