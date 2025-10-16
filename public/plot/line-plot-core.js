// Unified lightweight line plot core used by desktop energy plot and VR HUD.
// Responsibilities: maintain (steps,data) arrays, enforce maxPoints, draw axes + line.
// Environment: assumes a 2D canvas context returned by getContext(). Optional updateTexture() for dynamic textures.

export function createLinePlot(arg1, maybeOpts) {
  // Back-compat: previously called as createLinePlot(canvasElement, { ...options })
  // New preferred signature: createLinePlot({ width, height, getContext, updateTexture, ... })
  let opts;
  let canvasRef = null;
  if (arg1 && typeof arg1 === 'object' && typeof arg1.getContext === 'function' && (arg1.nodeName === 'CANVAS' || arg1.tagName === 'CANVAS')) {
    canvasRef = arg1;
    opts = { ...(maybeOpts || {}) };
    if (!opts.getContext) {
      opts.getContext = () => canvasRef.getContext('2d');
    }
    if (!opts.width) opts.width = canvasRef.width || 512;
    if (!opts.height) opts.height = canvasRef.height || 256;
  } else {
    opts = arg1 || {};
  }

  const width = opts.width || 512;
  const height = opts.height || 256;
  // Ensure getContext is a zero-arg function returning a 2D context; bind if user passed raw method.
  let getContext = opts.getContext;
  if (canvasRef && getContext === canvasRef.getContext) {
    // Wrap to preserve canvasRef as this and default to 2d.
    getContext = () => canvasRef.getContext('2d');
  }
  if (!getContext) {
    throw new Error('createLinePlot: getContext function required (provide opts.getContext or pass a <canvas>)');
  }
  const updateTexture = opts.updateTexture || (() => {});
  const labels = opts.labels || { x: 'X', y: 'Y' };
  const maxPoints = opts.maxPoints || 200;

  const steps = [];
  const data = [];
  let autoStep = 0;
  let visible = true;

  function addPoint(value, stepIndex) {
    if (typeof value !== 'number' || !isFinite(value)) return; // ignore non-finite
    const s = (typeof stepIndex === 'number') ? stepIndex : autoStep++;
    steps.push(s);
    data.push(value);
    if (data.length > maxPoints) {
      data.shift();
      steps.shift();
    }
    if (visible) redraw();
  }

  function reset() {
    steps.length = 0;
    data.length = 0;
    autoStep = 0;
    if (visible) redraw();
  }

  function setVisible(v) {
    visible = !!v;
    if (visible) redraw();
  }

  function computeBounds() {
    if (!data.length) return { min: 0, max: 1 };
    let min = data[0];
    let max = data[0];
    for (let i = 1; i < data.length; i++) {
      const d = data[i];
      if (d < min) min = d;
      if (d > max) max = d;
    }
    if (min === max) { // expand degenerate range
      min -= 0.5; max += 0.5;
    }
    return { min, max };
  }

  function redraw() {
    const ctx = getContext();
    if (!ctx) return;

    ctx.clearRect(0,0,width,height);
    ctx.fillStyle = '#000';
    ctx.fillRect(0,0,width,height);

    // padding for axes labels
    const padL = 50;
    const padB = 40;
    const plotW = width - padL - 10;
    const plotH = height - padB - 10;
    const originX = padL;
    const originY = height - padB;

    // axes
    ctx.strokeStyle = '#ccc';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(originX, originY);
    ctx.lineTo(originX + plotW, originY);
    ctx.moveTo(originX, originY);
    ctx.lineTo(originX, originY - plotH);
    ctx.stroke();

    // labels
    ctx.fillStyle = '#ccc';
    ctx.font = '20px sans-serif';
    ctx.fillText(labels.x, originX + plotW/2 - ctx.measureText(labels.x).width/2, originY + 30);
    ctx.save();
    ctx.translate(originX - 35, originY - plotH/2);
    ctx.rotate(-Math.PI/2);
    ctx.fillText(labels.y, -ctx.measureText(labels.y).width/2, 0);
    ctx.restore();

    if (data.length) {
      const { min, max } = computeBounds();
      const range = max - min;
      ctx.strokeStyle = '#6cf';
      ctx.lineWidth = 2;
      ctx.beginPath();
      for (let i = 0; i < data.length; i++) {
        const x = originX + (i / Math.max(1, data.length - 1)) * plotW;
        const yNorm = (data[i] - min) / range; // 0..1
        const y = originY - yNorm * plotH;
        if (i === 0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
      }
      ctx.stroke();
    }

    updateTexture();
  }

  return {
    addPoint,
    reset,
    redraw,
    setVisible,
    data: () => data.slice(),
    steps: () => steps.slice(),
    setAutoStep(v) { autoStep = v|0; },
  };
}

