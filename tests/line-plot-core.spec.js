// Updated to reference mlipviewer2 unified plot core (legacy public removed)
let createLinePlot;
beforeAll(async () => {
  ({ createLinePlot } = await import('../public/plot/line-plot-core.js'));
});

// Minimal canvas context stub implementing only what core uses
function makeCtx() {
  return {
    clearRect: () => {},
    fillRect: () => {},
    beginPath: () => {},
    moveTo: () => {},
    lineTo: () => {},
    stroke: () => {},
    fillText: () => {},
    measureText: (t) => ({ width: t.length * 8 }),
    save: () => {},
    restore: () => {},
    translate: () => {},
    rotate: () => {},
    // style props
    fillStyle: '',
    strokeStyle: '',
    lineWidth: 1,
    font: '',
  };
}

function makePlot(extra = {}) {
  const ctx = makeCtx();
  const plot = createLinePlot({
    width: 400,
    height: 200,
    getContext: () => ctx,
    updateTexture: () => {},
    maxPoints: extra.maxPoints || 5,
    labels: { x: 'Step', y: 'Energy' },
  });
  return plot;
}

describe('line-plot-core', () => {
  test('auto step increments when stepIndex omitted', () => {
    const p = makePlot();
    p.addPoint(10); // step 0
    p.addPoint(11); // step 1
    p.addPoint(12); // step 2
    expect(p.steps()).toEqual([0, 1, 2]);
    expect(p.data()).toEqual([10, 11, 12]);
  });

  test('explicit step indices preserved (non-contiguous)', () => {
    const p = makePlot();
    p.addPoint(1, 5);
    p.addPoint(2, 10);
    p.addPoint(3, 20);
    expect(p.steps()).toEqual([5, 10, 20]);
    expect(p.data()).toEqual([1, 2, 3]);
  });

  test('maxPoints rollover drops oldest', () => {
    const p = makePlot({ maxPoints: 3 });
    p.addPoint(1); // 0
    p.addPoint(2); // 1
    p.addPoint(3); // 2
    p.addPoint(4); // 3 -> drop first
    expect(p.data()).toEqual([2, 3, 4]);
    expect(p.steps()).toEqual([1, 2, 3]);
  });

  test('reset clears data and steps', () => {
    const p = makePlot();
    p.addPoint(5);
    p.addPoint(6);
    p.reset();
    expect(p.data()).toEqual([]);
    expect(p.steps()).toEqual([]);
    p.addPoint(7); // should start at step 0 again
    expect(p.steps()).toEqual([0]);
  });

  test('non-finite values are ignored', () => {
    const p = makePlot();
    p.addPoint(NaN);
    p.addPoint(Infinity);
    p.addPoint(-Infinity);
    p.addPoint(1.23);
    expect(p.data()).toEqual([1.23]);
    expect(p.steps()).toEqual([0]);
  });
});
