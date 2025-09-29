let createLinePlot;
beforeAll(async () => {
  ({ createLinePlot } = await import('../mlipviewer2/public/plot/line-plot-core.js'));
});

describe('energy plot DOM integration', () => {
  test('adds points without error', () => {
    const canvas = document.createElement('canvas');
    canvas.width = 200; canvas.height = 100;
    document.body.appendChild(canvas);
    const plot = createLinePlot(canvas, { maxPoints: 5 });
    plot.addPoint(0, 1.0);
    plot.addPoint(1, 2.0);
    expect(plot.data().length).toBe(2);
  });
});
