/** @jest-environment jsdom */
// Mock-browser test: loads the viewer, loads ROY by default, auto-starts MD, and processes >= 20 frames.

import fs from 'fs';
import path from 'path';
import { getWS } from '../public/fairchem_ws_client.js';

// Minimal Babylon stubs so render/scene.js doesn't blow up
jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn)=>{}, stopRenderLoop:()=>{} },
    scene: { meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){ this._l.push(fn); } } },
    camera: { attachControl:()=>{} },
  })
}));

describe('auto MD on page load (ROY)', () => {
  test('auto-starts MD and receives >= 20 frames', async () => {
    // Ensure auto-start is allowed: do NOT set __MLIPVIEW_TEST_MODE or __MLIPVIEW_NO_AUTO_MD
    // Stub the browser WebSocket to avoid real network and trigger onopen.
    const origWS = global.WebSocket;
    let constructed = 0;
    class FakeWS {
      constructor(url){ constructed++; this.url=url; this.readyState=0; setTimeout(()=>{ this.readyState=1; this.onopen && this.onopen(); }, 0); }
      set binaryType(_){}
      send(){}
      close(){}
      onopen(){}
      onerror(){}
      onclose(){}
      set onmessage(fn){ this._onmsg = fn; }
      get onmessage(){ return this._onmsg; }
    }
    global.WebSocket = FakeWS;

    // DOM scaffold similar to index.html
    document.body.innerHTML = '';
    const canvas = document.createElement('canvas'); canvas.id='viewer'; canvas.addEventListener=()=>{}; document.body.appendChild(canvas);
    const hud = document.createElement('div'); hud.className='hud'; document.body.appendChild(hud);
    const status = document.createElement('span'); status.id='status'; hud.appendChild(status);
    const rps = document.createElement('span'); rps.id='rpsLabel'; rps.textContent='RPS: --'; hud.appendChild(rps);
    const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
    const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
    const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);

    // Serve ROY xyz via fetch
    const royPath = path.resolve(process.cwd(), 'public/molecules/roy.xyz');
    const royText = fs.readFileSync(royPath, 'utf8');
    const origFetch = global.fetch;
    global.fetch = async (url) => {
      if (typeof url === 'string' && /\/molecules\/roy\.xyz$/.test(url)) {
        return { ok:true, status:200, text: async()=> royText };
      }
      // Fallback: 404 for unexpected fetches in this test
      return { ok:false, status:404, text: async()=> 'not found' };
    };

    try {
      const { initNewViewer } = await import('../public/index.js');
      const { loadDefault } = await import('../public/util/moleculeLoader.js');

      // Initialize viewer with empty placeholders like index.html
      const api = await initNewViewer(canvas, { elements:[], positions:[], bonds:[] });
      // Expose for auto-start path that checks window.viewerApi
      window.viewerApi = api;

      // Capture outgoing WS actions using instance-level hook
      const ws = getWS();
      const sent = [];
      ws.setTestHook(msg => sent.push(msg));

      // Load ROY by default using loader (mirrors index.js behavior)
      const ld = await loadDefault(api);
      expect(ld && ld.file).toBe('molecules/roy.xyz');

      // Allow auto-start tick to schedule startMDContinuous and connection to open
      await Promise.resolve();
      await new Promise(r=>setTimeout(r,0));
      await new Promise(r=>setTimeout(r,0));

      // Wait briefly for a START_SIMULATION for MD to be sent
      let started = false;
      for(let k=0;k<10 && !started;k++){
        started = sent.some(m => m && m.simulationParams && m.simulationType != null);
        if(!started) await new Promise(r=>setTimeout(r,5));
      }
      expect(started).toBe(true);

      // Emit >= 20 frames with monotonically changing energy so the energy plot records ticks
      const N = api.state.positions.length;
      const basePos = api.state.positions.map(p=>[p.x,p.y,p.z]);
      for (let i=0;i<20;i++){
        const pos = basePos.map((p,idx)=> [p[0] + (i*0.0001) + (idx===0? i*0.00005:0), p[1], p[2]]);
        const forces = Array.from({length:N}, (_,k)=> [0.001*((k+i)%3+1), 0, 0]);
        ws.injectTestResult({ positions: pos, forces, energy: -100 - i, temperature: 300 + i*0.1 });
        await new Promise(r=>setTimeout(r,0));
      }

      // Let UI settle then assert we got at least 20 energy ticks
      await Promise.resolve();
      const ticks = api.debugEnergySeriesLength();
      expect(ticks).toBeGreaterThanOrEqual(20);
    } finally {
      global.fetch = origFetch;
      global.WebSocket = origWS;
    }
  }, 15000);
});
