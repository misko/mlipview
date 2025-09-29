import { createLinePlot } from '../plot/line-plot-core.js';
// Lightweight debug hook (set window.DEBUG_PLOT = true to enable logging if in browser)
const dbg = (...args) => { try { if (typeof window !== 'undefined' && window.DEBUG_PLOT) console.log('[plot]', ...args); } catch {} };

export function createEnergyPlot() {
  const plotContainer = (typeof document !== 'undefined') ? document.getElementById('plotContainer') : null;
  const energyCanvas = (typeof document !== 'undefined') ? document.getElementById('energyChart') : null;

  if (!plotContainer || !energyCanvas) {
    return { isEnabled: () => false, toggle: () => false, addInitialDataPoint: ()=>{}, recordStep: ()=>{}, clear: ()=>{}, destroy: ()=>{}, _debug: () => ({ data: [], steps: [] }) };
  }

  energyCanvas.width = 1024; energyCanvas.height = 512;
  let showPlot = true;
  plotContainer.style.display = 'block';
  let internalStep = 0;
  const linePlot = createLinePlot({
    width: energyCanvas.width,
    height: energyCanvas.height,
    getContext: () => energyCanvas.getContext('2d'),
    labels: { x: 'Step', y: 'Energy' },
    maxPoints: 500
  });

  function recordStep(energy) { internalStep += 1; linePlot.addPoint(energy, internalStep); }
  function addInitialDataPoint(energy) { internalStep = 0; linePlot.addPoint(energy, internalStep); }
  function toggle() { showPlot = !showPlot; plotContainer.style.display = showPlot ? 'block' : 'none'; linePlot.setVisible(showPlot); return showPlot; }
  function clear() { internalStep = 0; linePlot.reset(); }
  function destroy() {}
  linePlot.redraw();
  return { isEnabled: () => showPlot, toggle, recordStep, addInitialDataPoint, clear, destroy, _debug: () => ({ data: linePlot.data(), steps: linePlot.steps() }) };
}

