import { createEnergyPlot } from '../public/ui/energy-plot.js';

// Provide minimal DOM & canvas stubs without pulling full jsdom (keeps Node 18 compatibility)
function setupDom() {
  const elements = new Map();
  function makeCanvas(id) {
    const canvas = { id, getContext: () => ({}) };
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
  const charts = [];
  global.Chart = class { constructor(ctx,cfg){ this.data=cfg.data; charts.push(this);} update(){} static getChart(){ return null;} };
  return charts;
}

describe('energy plot step counter', () => {
  it('increments labels and data with recordStep', () => {
    const charts = setupDom();
    const plot = createEnergyPlot();
    plot.addInitialDataPoint(5); // step 0
    plot.recordStep(6); // step 1
    plot.recordStep(7); // step 2
    const chart = charts[0];
    expect(chart.data.labels).toEqual([0,1,2]);
    expect(chart.data.datasets[0].data).toEqual([5,6,7]);
  });
});
