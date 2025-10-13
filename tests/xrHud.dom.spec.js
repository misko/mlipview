/** @jest-environment jsdom */

// Note: We don't instantiate Babylon in jsdom. This test validates the wiring logic
// by simulating button callbacks and ensuring correct calls into the core model.

describe('XR HUD button semantics (simulated)', () => {
  test('Buttons call core model as expected', async () => {
    const { buildXRControlsModel } = await import('../public/vr/xr-controls-core.js');
    const calls = [];
    const viewer = {
      state: { showForces: false, bus: { on: jest.fn() }, toggleForceVectorsVisibility(){ this.showForces = !this.showForces; } },
      getMetrics: () => ({ running: null }),
      startRelaxContinuous: jest.fn(async()=>{ calls.push('relax'); return {}; }),
      startMDContinuous: jest.fn(async()=>{ calls.push('md'); return {}; }),
      stopSimulation: jest.fn(()=>{ calls.push('off'); })
    };
    const model = buildXRControlsModel({ getViewer: ()=> viewer });

    // Simulate pressing Relax
    model.setSimulation('relax');
    expect(calls).toContain('relax');
    viewer.getMetrics = () => ({ running: 'relax' });
    expect(model.simSelection()).toEqual({ relax:true, md:false, off:false });

    // Simulate pressing MD
    model.setSimulation('md');
    expect(calls).toContain('md');
    viewer.getMetrics = () => ({ running: 'md' });
    expect(model.simSelection()).toEqual({ relax:false, md:true, off:false });

    // Simulate pressing Off
    model.setSimulation('off');
    expect(calls).toContain('off');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax:false, md:false, off:true });

    // Forces label flips on toggle
    expect(model.forcesLabel()).toBe('Forces On');
    model.toggleForces();
    expect(model.forcesLabel()).toBe('Forces Off');
  });
});
