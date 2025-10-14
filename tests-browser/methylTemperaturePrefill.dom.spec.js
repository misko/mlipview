/** @jest-environment jsdom */

// Recreates loading the methyl radicals example (twomethane.xyz) and verifies
// the MD temperature defaults to 500K and the slider reflects it.

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';
import { initTemperatureSlider } from '../public/ui/temperatureSlider.js';

function makeViewer(){
  const listeners = {};
  const bus = { on:(ev,fn)=>{ (listeners[ev]||(listeners[ev]=[])).push(fn); }, emit:(ev)=>{ (listeners[ev]||[]).forEach(fn=>fn()); } };
  const state = {
    bus,
    elements: [], positions: [], bonds: [],
    cell: { a:{x:1,y:0,z:0}, b:{x:0,y:1,z:0}, c:{x:0,y:0,z:1}, enabled:false, originOffset:{x:0,y:0,z:0} },
    showCell: false,
    markCellChanged(){ bus.emit('cellChanged'); },
    markPositionsChanged(){ bus.emit('positionsChanged'); },
    markBondsChanged(){ bus.emit('bondsChanged'); },
    dynamics: {}
  };
  const viewerApi = { state, recomputeBonds: ()=>{} };
  return viewerApi;
}

beforeEach(()=>{ document.body.innerHTML = '<div id="hud"></div>'; delete window.__MLIP_TARGET_TEMPERATURE; });

function readFixture(){
  // Inline a minimal shard of twomethane.xyz header with temperature, followed by a few atoms
  return `24\ntemp=500K (methyl radicals, reactive) \nC    0.100   0.000   0.000 \nH    0.100   1.048   0.300 \nH    1.009  -0.524  -0.300 \nH   -0.809  -0.524   0.300 \nC    2.600   0.000   0.000 \nH    2.600   1.048  -0.300 \nH    3.509  -0.524   0.300 \nH    1.691  -0.524  -0.300 \nC    1.350   2.100   0.000 \nH    1.350   3.148   0.300 \nH    2.259   1.576  -0.300 \nH    0.441   1.576   0.300 \nC    1.350  -2.100   0.000 \nH    1.350  -1.052  -0.300 \nH    2.259  -2.624   0.300 \nH    0.441  -2.624  -0.300 \nC    1.350   0.000   2.100 \nH    1.350   1.048   2.400 \nH    2.259  -0.524   1.800 \nH    0.441  -0.524   2.400 \nC    1.350   0.000  -2.100 \nH    1.350   1.048  -2.400 \nH    2.259  -0.524  -1.800 \nH    0.441  -0.524  -2.400 \n`;
}

describe('methyl radicals example initializes temperature to 500K', ()=>{
  test('load XYZ -> state + slider at 500K', ()=>{
    const xyz = readFixture();
    const parsed = parseXYZ(xyz);
    expect(parsed.temperature).toBe(500);

    const viewer = makeViewer();
    const hud = document.getElementById('hud');
    // Build slider FIRST (as in app init), then apply load which dispatches a sync event
    const slider = initTemperatureSlider({ hudEl: hud, getViewer: ()=>viewer });
    applyParsedToViewer(viewer, parsed);

    // UI elements created
    const label = hud.querySelector('#tempLabel');
    const input = hud.querySelector('#mdTempSlider');
    expect(label && input).toBeTruthy();

    // Label shows 500K
    expect(label.textContent).toMatch(/500\s*K/);

    // Slider's current selection corresponds to 500K in its tick list
    // We use exposed helper to get the temperature
    expect(slider.getTemperature()).toBe(500);

    // State should carry targetTemperature as well
    expect(viewer.state.dynamics.targetTemperature).toBe(500);
    // Global mirrors 500K
    expect(window.__MLIP_TARGET_TEMPERATURE).toBe(500);
  });
});
