/** @jest-environment jsdom */

describe('x-xr-controls-core', () => {
  test('forces toggle and simulation selection wiring', async () => {
    const { buildXRControlsModel } = await import('../public/vr/xr-controls-core.js');

    const viewer = {
      state: {
        showForces: false,
        toggleForceVectorsVisibility() {
          this.showForces = !this.showForces;
        },
      },
      getMetrics: () => ({ running: null }),
      startRelaxContinuous: jest.fn(async () => ({})),
      startMDContinuous: jest.fn(async () => ({})),
      stopSimulation: jest.fn(() => {}),
    };

    const model = buildXRControlsModel({ getViewer: () => viewer });

    expect(model.simSelection()).toEqual({ relax: false, md: false, off: true });
    expect(model.isForcesOn()).toBe(false);

    model.toggleForces();
    expect(viewer.state.showForces).toBe(true);
    expect(model.isForcesOn()).toBe(true);

    model.setSimulation('relax');
    viewer.getMetrics = () => ({ running: 'relax' });
    expect(model.simSelection()).toEqual({ relax: true, md: false, off: false });

    model.setSimulation('relax');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax: false, md: false, off: true });

    model.setSimulation('md');
    viewer.getMetrics = () => ({ running: 'md' });
    expect(model.simSelection()).toEqual({ relax: false, md: true, off: false });

    model.setSimulation('md');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax: false, md: false, off: true });

    model.setSimulation('off');
    expect(model.simSelection()).toEqual({ relax: false, md: false, off: true });

    viewer.getMetrics = () => ({ running: null });
    model.ensureDefaultActive();
    expect(viewer.startMDContinuous).toHaveBeenCalled();
  });
});

