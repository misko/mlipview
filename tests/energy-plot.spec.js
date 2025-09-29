// Updated to mlipviewer2 implementation
let createEnergyPlot;
beforeAll(async () => {
  ({ createEnergyPlot } = await import('../public/ui/energy-plot.js'));
});

// Minimal DOM & canvas stubs (no jsdom) for Node environment.
function setupDom() {
  const elements = new Map();
  function makeCanvas(id) {
    const canvas = { id, width: 0, height: 0, getContext: () => ({
      clearRect: ()=>{}, fillRect: ()=>{}, beginPath: ()=>{}, moveTo: ()=>{}, lineTo: ()=>{}, stroke: ()=>{}, fillText: ()=>{}, measureText: (t)=>({ width: t.length * 10 }), save: ()=>{}, restore: ()=>{}, translate: ()=>{}, rotate: ()=>{}, font: '', fillStyle: '', strokeStyle: '', lineWidth: 1
    }) };
    elements.set(id, canvas);
    return canvas;
  }
  function makeDiv(id) {
    const div = { id, style: {}, appendChild: ()=>{}, setAttribute: ()=>{} };
    elements.set(id, div);
    return div;
  }
  function makeButton(id){ const b={ id, onclick:null }; elements.set(id,b); return b; }
  makeDiv('plotContainer');
  makeCanvas('energyChart');
  makeButton('clearPlot');
  global.document = { getElementById: (id) => elements.get(id) };
  global.window = { location: { search: '' } };
}

describe('energy plot step counter (unified line plot)', () => {
  it('increments steps and values with recordStep', () => {
    setupDom();
    const plot = createEnergyPlot();
    plot.addInitialDataPoint(5); // step 0
    plot.recordStep(6); // step 1
    plot.recordStep(7); // step 2
    const dbg = plot._debug();
    expect(dbg.steps).toEqual([0,1,2]);
    expect(dbg.data).toEqual([5,6,7]);
  });
});
