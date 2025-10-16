import { initNewViewer } from '../public/index.js';

// Provide minimal DOM canvas for scene
function makeCanvas(){ return { getBoundingClientRect(){ return { left:0, top:0, width:800, height:600 }; }, addEventListener(){}, removeEventListener(){} }; }

// Intercept fetch calls to capture payload
let lastBody = null; let lastUrl = null;
beforeEach(()=>{ lastBody = null; lastUrl = null; global.__MLIP_API_URL = 'http://127.0.0.1:8000'; global.fetch = async (url, opts)=>{ lastUrl = url; try { lastBody = JSON.parse(opts.body); } catch { lastBody = null; } return { ok:true, status:200, json: async ()=> ({ results:{ energy: -1.23, forces:[[0,0,0]] } }) }; } });

describe('API includes cell when PBC enabled', () => {
  test('simple_calculate body has cell after toggle on', async () => {
    const canvas = makeCanvas();
    const api = await initNewViewer(canvas, { elements:['H','O','H'], positions:[[0,0,0],[0.96,0,0],[-0.24,0.93,0]], bonds:[] });
    // Enable PBC via enhanced toggle
    api.state.toggleCellVisibilityEnhanced();
    // Trigger a force compute which sends /serve/simple
    await api.ff.computeForces({ sync:true });
    expect(lastUrl).toMatch(/\/serve\/simple$/);
    expect(lastBody).toBeTruthy();
    expect(lastBody.cell).toBeTruthy();
    expect(Array.isArray(lastBody.cell)).toBe(true);
    expect(lastBody.pbc).toEqual([true,true,true]);
  });
});
