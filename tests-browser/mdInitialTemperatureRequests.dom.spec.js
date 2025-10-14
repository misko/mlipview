/** @jest-environment jsdom */

// Validate that the first 10 MD API requests on page load (continuous run)
// use the temperature from the XYZ comment if present.

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';

// Minimal BABYLON shims to satisfy viewer init in jsdom tests that exercise MD logic indirectly.
beforeAll(()=>{ if(!global.BABYLON){ global.BABYLON={ TransformNode:class{}, MeshBuilder:{ CreateCylinder:()=>({ dispose(){}, setEnabled(){}, position:{ set(){} }, rotationQuaternion:null, scaling:{}, isPickable:false, visibility:1 }) }, StandardMaterial:class{constructor(){this.diffuseColor={};this.emissiveColor={};this.specularColor={};}}, Color3:class{}, Vector3:class{ constructor(x=0,y=0,z=0){this.x=x;this.y=y;this.z=z;} static Up(){return new global.BABYLON.Vector3(0,1,0);} }, Quaternion:class{ static Identity(){return{};}}, Scene:class{} }; } });

// Mock the scene builder so index.js can initialize without a full WebGL context
jest.mock('../public/render/scene.js', ()=>({ createScene: async ()=>({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(f){this._l.push(f);} } }, camera:{ attachControl:()=>{} } }) }));

function buildViewerCanvas(){
  const canvas=document.createElement('canvas'); canvas.id='viewer';
  const energyDiv=document.createElement('div'); energyDiv.id='energyPlot';
  const c=document.createElement('canvas'); c.id='energyCanvas'; c.width=200; c.height=60; c.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}});
  energyDiv.appendChild(c);
  document.body.appendChild(canvas); document.body.appendChild(energyDiv);
}

function fixtureXYZ(){
  return `4\ntemp=500K\nH 0 0 0\nH 0 0 1\nH 0 1 0\nH 1 0 0\n`;
}

describe('MD initial requests honor XYZ temperature', ()=>{
  test('first 10 MD requests use 500K from XYZ', async ()=>{
    // Arrange test mode and canvas
    window.__MLIPVIEW_TEST_MODE = true; // bypass focus gating
    window.__MLIPVIEW_NO_AUTO_MD = true; // prevent auto-start path from index.html
    buildViewerCanvas();

    // Prepare a fake fetch to capture MD requests
    const requests = [];
    const origFetch = global.fetch;
    global.fetch = async (url, opts={})=>{
      try {
        const body = opts && opts.body ? JSON.parse(opts.body) : {};
        if (/\/md$/.test(String(url)) && body && typeof body.temperature !== 'undefined') {
          requests.push(body.temperature);
        }
        // Return a minimal successful MD response
        if (/\/md$/.test(String(url))) {
          return { ok:true, status:200, json: async()=>({ positions:[], forces:[], final_energy:0, velocities:[], temperature: body.temperature||0 }) };
        }
      } catch {}
      return { ok:true, status:200, json: async()=>({ ok:true }), text: async()=>('ok') };
    };

    // Initialize viewer API directly (index.js exports initNewViewer)
    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(document.getElementById('viewer'), { elements:[{Z:1}], positions:[{x:0,y:0,z:0}], bonds:[] });

    // Apply parsed XYZ with temp=500K to seed defaults
    const parsed = parseXYZ(fixtureXYZ());
    applyParsedToViewer(api, parsed);
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(500);

    // Run MD continuously for 10 steps
    const res = await api.startMDContinuous({ steps:10 });
    expect(res.completed).toBe(true);

    // Assert the first 10 requests used 500K
    expect(requests.length).toBeGreaterThanOrEqual(10);
    expect(requests.slice(0,10).every(t=>t===500)).toBe(true);

    // Cleanup fetch
    if (origFetch) global.fetch = origFetch;
  });
});
