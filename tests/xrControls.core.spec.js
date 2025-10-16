/** @jest-environment jsdom */

describe('XR Controls Core Model', () => {
  test('simulation selection and forces toggle state', async () => {
    const { buildXRControlsModel } = await import('../public/vr/xr-controls-core.js');
    const viewer = {
      state: { showForces: false, toggleForceVectorsVisibility: function(){ this.showForces = !this.showForces; } },
      getMetrics: () => ({ running: null }),
      startRelaxContinuous: jest.fn(async()=>({})),
      startMDContinuous: jest.fn(async()=>({})),
      stopSimulation: jest.fn(()=>{})
    };
    const model = buildXRControlsModel({ getViewer: ()=> viewer });
    // Initial
    expect(model.simSelection()).toEqual({ relax:false, md:false, off:true });
    // Forces initially off
    expect(model.isForcesOn()).toBe(false);
    // Toggle forces -> now on
    model.toggleForces();
    expect(viewer.state.showForces).toBe(true);
    expect(model.isForcesOn()).toBe(true);
    // Start relax
    model.setSimulation('relax');
    viewer.getMetrics = () => ({ running: 'relax' });
    expect(model.simSelection()).toEqual({ relax:true, md:false, off:false });
    // Press relax again -> off
    model.setSimulation('relax');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax:false, md:false, off:true });
    // Switch to md
    model.setSimulation('md');
    viewer.getMetrics = () => ({ running: 'md' });
    expect(model.simSelection()).toEqual({ relax:false, md:true, off:false });
    // Press md again -> off
    model.setSimulation('md');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax:false, md:false, off:true });
    // Turn off
    model.setSimulation('off');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax:false, md:false, off:true });
    // ensureDefaultActive sets MD when off
    viewer.getMetrics = () => ({ running: null });
    model.ensureDefaultActive();
    expect(viewer.startMDContinuous).toHaveBeenCalled();
  });
});
