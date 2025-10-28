/** @jest-environment jsdom */

describe('x-xr-hud-core-buttons', () => {
  test('button semantics call xr controls model hooks', async () => {
    const { buildXRControlsModel } = await import('../public/vr/xr-controls-core.js');
    const calls = [];

    const viewer = {
      state: {
        showForces: false,
        bus: { on: jest.fn() },
        toggleForceVectorsVisibility() {
          this.showForces = !this.showForces;
        },
      },
      getMetrics: () => ({ running: null }),
      startRelaxContinuous: jest.fn(async () => {
        calls.push('relax');
        return {};
      }),
      startMDContinuous: jest.fn(async () => {
        calls.push('md');
        return {};
      }),
      stopSimulation: jest.fn(() => {
        calls.push('off');
      }),
    };

    const model = buildXRControlsModel({ getViewer: () => viewer });

    model.setSimulation('relax');
    expect(calls).toContain('relax');
    viewer.getMetrics = () => ({ running: 'relax' });
    expect(model.simSelection()).toEqual({ relax: true, md: false, off: false });

    model.setSimulation('md');
    expect(calls).toContain('md');
    viewer.getMetrics = () => ({ running: 'md' });
    expect(model.simSelection()).toEqual({ relax: false, md: true, off: false });

    model.setSimulation('off');
    expect(calls).toContain('off');
    viewer.getMetrics = () => ({ running: null });
    expect(model.simSelection()).toEqual({ relax: false, md: false, off: true });

    expect(model.isForcesOn()).toBe(false);
    model.toggleForces();
    expect(model.isForcesOn()).toBe(true);
  });
});
