/** @jest-environment jsdom */
// Integration: Only API energy responses add ticks; drag does not.

const { haveServer, getApiBase } = require('../tests/helpers/server.js');

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

describe('[integration] API-only energy ticks (live WS)', () => {
  test('drag does not tick; relax+md do', async () => {
    if (process.env.MLIPVIEW_SKIP_SERVERS === '1') { console.warn('[integration] skipping due to MLIPVIEW_SKIP_SERVERS'); return; }
    const up = await haveServer(); if(!up){ console.warn('[integration] backend not reachable; skipping'); return; }
    if (typeof global.WebSocket === 'undefined' && typeof window.WebSocket === 'undefined') { console.warn('[integration] WebSocket API not available; skipping'); return; }
  const base = getApiBase(); const u = new URL(base);
    delete window.location; window.location = { protocol: u.protocol, host: u.host, origin: u.origin, search: '', href: u.href };
    window.__MLIPVIEW_SERVER = base;
  window.__MLIPVIEW_NO_AUTO_MD = true;
  window.__MLIPVIEW_TEST_MODE = true;

    document.body.innerHTML = '';
    const canvas = document.createElement('canvas'); canvas.id='viewer'; canvas.addEventListener=()=>{}; document.body.appendChild(canvas);
    const hud = document.createElement('div'); hud.className='hud'; document.body.appendChild(hud);
    const energyWrapper=document.createElement('div'); energyWrapper.id='energyPlot'; document.body.appendChild(energyWrapper);
    const energyCanvas=document.createElement('canvas'); energyCanvas.id='energyCanvas'; energyCanvas.width=260; energyCanvas.height=80; energyCanvas.getContext=()=>({ clearRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, arc(){}, fill(){}, fillRect(){}, strokeStyle:null, lineWidth:1, fillStyle:null }); energyWrapper.appendChild(energyCanvas);
    const energyLabel=document.createElement('div'); energyLabel.id='energyLabel'; energyWrapper.appendChild(energyLabel);

    const { initNewViewer } = await import('../public/index.js');
  const api = await initNewViewer(canvas, { elements:[6,6], positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], bonds:[{i:0,j:1}] });

    // Give baseline a moment
    await new Promise(r=>setTimeout(r,100));
    const start = api.debugEnergySeriesLength();

    // Drag (no tick)
    api.state.selection = { kind:'atom', data:{ index:1 } };
    api.manipulation.beginDrag(()=>({ x:1.4,y:0,z:0 }));
    api.manipulation.updateDrag(()=>({ x:1.6,y:0,z:0 }));
    api.manipulation.endDrag();
    await new Promise(r=>setTimeout(r,200));
    const afterDrag = api.debugEnergySeriesLength();
    expect(afterDrag).toBe(start);

    // Relax step (adds tick)
  await new Promise(r=>setTimeout(r,100));
  await api.relaxStep();
    await new Promise(r=>setTimeout(r,200));
    const afterRelax = api.debugEnergySeriesLength();
    expect(afterRelax).toBeGreaterThan(afterDrag);

    // MD step (adds tick)
  await new Promise(r=>setTimeout(r,100));
  await api.mdStep();
    await new Promise(r=>setTimeout(r,200));
    const afterMd = api.debugEnergySeriesLength();
    expect(afterMd).toBeGreaterThan(afterRelax);
  }, 30000);
});
