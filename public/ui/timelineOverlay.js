// Timeline overlay UI — hover-revealed DVR controls.

let stylesInstalled = false;

function ensureStyles() {
  if (stylesInstalled) return;
  try {
    const style = document.createElement('style');
    style.id = 'mlip-timeline-styles';
    style.textContent = `
      .mlip-timeline-rail {
        position: absolute;
        left: 50%;
        bottom: 12px;
        transform: translateX(-50%);
        width: min(720px, 80%);
        pointer-events: auto;
        background: rgba(10, 14, 22, 0.78);
        border: 1px solid rgba(255,255,255,0.12);
        border-radius: 14px;
        padding: 6px 14px 12px;
        display: flex;
        align-items: center;
        gap: 12px;
        height: 10px;
        opacity: 0.2;
        transition: height 0.18s ease, opacity 0.18s ease;
        box-shadow: 0 4px 14px rgba(0,0,0,0.35);
        backdrop-filter: blur(8px);
      }
      .mlip-timeline-rail[data-expanded="true"] {
        height: 74px;
        opacity: 1;
      }
      .mlip-timeline-inner {
        display: flex;
        align-items: center;
        gap: 12px;
        width: 100%;
        opacity: 0;
        transition: opacity 0.18s ease;
      }
      .mlip-timeline-rail[data-expanded="true"] .mlip-timeline-inner {
        opacity: 1;
      }
      .mlip-timeline-label {
        min-width: 140px;
        font-size: 12px;
        color: rgba(255,255,255,0.85);
        white-space: nowrap;
        font-variant-numeric: tabular-nums;
      }
      .mlip-timeline-slider {
        flex: 1;
        accent-color: #7fc4ff;
      }
      .mlip-timeline-buttons {
        display: inline-flex;
        gap: 6px;
      }
      .mlip-timeline-buttons button {
        background: rgba(255,255,255,0.12);
        border: 1px solid rgba(255,255,255,0.2);
        color: rgba(255,255,255,0.92);
        border-radius: 999px;
        padding: 4px 10px;
        font-size: 12px;
        cursor: pointer;
        transition: background 0.18s ease, color 0.18s ease;
      }
      .mlip-timeline-buttons button[data-active="true"] {
        background: rgba(127,196,255,0.9);
        color: #051626;
      }
      .mlip-timeline-buttons button:disabled {
        opacity: 0.35;
        cursor: default;
      }
      .mlip-timeline-hint {
        position: absolute;
        bottom: -4px;
        left: 50%;
        transform: translateX(-50%);
        width: 34px;
        height: 6px;
        border-radius: 3px;
        background: rgba(255,255,255,0.5);
        transition: opacity 0.2s ease;
        pointer-events: none;
      }
      .mlip-timeline-rail[data-expanded="true"] .mlip-timeline-hint {
        opacity: 0;
      }
      .mlip-timeline-heat {
        position: absolute;
        inset: 0;
        border-radius: inherit;
        overflow: hidden;
        pointer-events: none;
        opacity: 0.28;
      }
      .mlip-timeline-heat canvas {
        width: 100%;
        height: 100%;
        display: block;
      }
    `;
    document.head?.appendChild(style);
    stylesInstalled = true;
  } catch {
    stylesInstalled = true;
  }
}

function formatDelta(liveFrame, frame) {
  if (!liveFrame || !frame) return '—';
  const delta = Number(liveFrame.receivedAt || 0) - Number(frame.receivedAt || 0);
  if (!Number.isFinite(delta)) return '—';
  if (delta <= 80) return 'live';
  const secs = Math.max(delta / 1000, 0);
  if (secs >= 120) return `-${Math.round(secs)}s`;
  return `-${secs.toFixed(1)}s`;
}

function frameLabel(frame) {
  if (!frame) return 'Frame —';
  const tag = frame.kind ? frame.kind.toUpperCase() : '—';
  return `Frame ${frame.id ?? '?'} • ${tag}`;
}

export function installTimelineOverlay({
  container,
  controller,
  onPause,
  onPlay,
  onLive,
  onScrub,
} = {}) {
  if (!container || !controller) return null;
  ensureStyles();

  const isDebug = () => {
    try { return typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API; } catch { return false; }
  };

  const rail = document.createElement('div');
  rail.className = 'mlip-timeline-rail';
  rail.dataset.expanded = 'false';
  rail.tabIndex = 0;
  rail.setAttribute('role', 'group');
  rail.setAttribute('aria-label', 'Simulation timeline controls');

  const hint = document.createElement('div');
  hint.className = 'mlip-timeline-hint';
  rail.appendChild(hint);

  const heat = document.createElement('div');
  heat.className = 'mlip-timeline-heat';
  const heatCanvas = document.createElement('canvas');
  heat.appendChild(heatCanvas);
  rail.appendChild(heat);

  const inner = document.createElement('div');
  inner.className = 'mlip-timeline-inner';

  const label = document.createElement('span');
  label.className = 'mlip-timeline-label';
  label.textContent = 'Frame —';

  const slider = document.createElement('input');
  slider.type = 'range';
  slider.min = '0';
  slider.max = '0';
  slider.step = '1';
  slider.value = '0';
  slider.disabled = true;
  slider.className = 'mlip-timeline-slider';
  slider.setAttribute('aria-label', 'Timeline frame');
  const ticks = document.createElement('datalist');
  const listId = 'mlip-timeline-ticks';
  ticks.id = listId;
  slider.setAttribute('list', listId);

  const buttons = document.createElement('div');
  buttons.className = 'mlip-timeline-buttons';

  const playBtn = document.createElement('button');
  playBtn.type = 'button';
  playBtn.textContent = 'Play';
  playBtn.setAttribute('aria-label', 'Play from selected frame');

  const pauseBtn = document.createElement('button');
  pauseBtn.type = 'button';
  pauseBtn.textContent = 'Pause';
  pauseBtn.setAttribute('aria-label', 'Pause simulation');

  const liveBtn = document.createElement('button');
  liveBtn.type = 'button';
  liveBtn.textContent = 'Live';
  liveBtn.setAttribute('aria-label', 'Jump to live frame');

  buttons.appendChild(playBtn);
  buttons.appendChild(pauseBtn);
  buttons.appendChild(liveBtn);

  inner.appendChild(label);
  inner.appendChild(slider);
  inner.appendChild(buttons);

  rail.appendChild(inner);
  rail.appendChild(ticks);
  container.appendChild(rail);

  let collapseTimer = null;
  const expand = () => {
    rail.dataset.expanded = 'true';
    if (collapseTimer) {
      clearTimeout(collapseTimer);
      collapseTimer = null;
    }
  };
  const collapse = () => {
    if (rail.contains(document.activeElement)) return;
    rail.dataset.expanded = 'false';
  };
  const scheduleCollapse = () => {
    if (collapseTimer) clearTimeout(collapseTimer);
    collapseTimer = setTimeout(() => {
      collapseTimer = null;
      collapse();
    }, 700);
  };

  rail.addEventListener('mouseenter', expand);
  rail.addEventListener('mouseleave', scheduleCollapse);
  rail.addEventListener('focusin', expand);
  rail.addEventListener('focusout', scheduleCollapse);

  let latestSnapshot = controller.debug?.() || null;
  let bufferSnapshot = null;

  function setActiveButton(btn) {
    for (const child of buttons.children) {
      if (child instanceof HTMLElement) child.dataset.active = child === btn ? 'true' : 'false';
    }
  }

  function updateButtonStates(mode) {
    if (mode === 'playback') {
      playBtn.disabled = true;
      pauseBtn.disabled = false;
      liveBtn.disabled = false;
    } else if (mode === 'paused') {
      playBtn.disabled = false;
      pauseBtn.disabled = true;
      liveBtn.disabled = false;
    } else {
      playBtn.disabled = false;
      pauseBtn.disabled = false;
      liveBtn.disabled = false;
    }
  }

  function updateLabel(index) {
    if (!bufferSnapshot || !bufferSnapshot.frames || bufferSnapshot.frames.length === 0) {
      label.textContent = 'Frame —';
      slider.setAttribute('aria-valuetext', 'none');
      return;
    }
    const frame = bufferSnapshot.frames[index] || null;
    const liveFrame =
      bufferSnapshot.liveIndex >= 0 ? bufferSnapshot.frames[bufferSnapshot.liveIndex] : null;
    label.textContent = `${frameLabel(frame)} • ${formatDelta(liveFrame, frame)}`;
    slider.setAttribute(
      'aria-valuetext',
      `${frame ? frame.id : '—'} (${formatDelta(liveFrame, frame)})`,
    );
  }

  function drawHeat(snapshot) {
    if (!snapshot || !snapshot.frames || !snapshot.frames.length) return;
    const frames = snapshot.frames;
    try {
      heatCanvas.width = heat.clientWidth || rail.clientWidth;
      heatCanvas.height = heat.clientHeight || rail.clientHeight;
      const ctx = heatCanvas.getContext('2d');
      if (!ctx) return;
      ctx.clearRect(0, 0, heatCanvas.width, heatCanvas.height);
      const w = heatCanvas.width;
      const h = heatCanvas.height;
      const count = frames.length;
      for (let i = 0; i < count; i++) {
        const frame = frames[i];
        let color = 'rgba(120, 120, 120, 0.65)';
        if (frame.kind === 'md') color = 'rgba(255, 166, 0, 0.55)';
        else if (frame.kind === 'relax') color = 'rgba(0, 188, 255, 0.55)';
        const x0 = Math.floor((i / count) * w);
        const x1 = Math.floor(((i + 1) / count) * w);
        ctx.fillStyle = color;
        ctx.fillRect(x0, 0, Math.max(1, x1 - x0), h);
      }
    } catch {
      /* ignore */
    }
  }

  let pendingScrubIndex = null;
  let scrubScheduled = false;
  let lastInputTs = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();

  const scheduleScrub = () => {
    if (scrubScheduled) return;
    scrubScheduled = true;
    const schedule = typeof requestAnimationFrame === 'function' ? requestAnimationFrame : (fn) => setTimeout(fn, 16);
    schedule(() => {
      scrubScheduled = false;
      const idx = pendingScrubIndex;
      pendingScrubIndex = null;
      if (idx == null) return;
      controller.scrubTo?.(idx);
      setActiveButton(pauseBtn);
      if (typeof onScrub === 'function') {
        onScrub(idx);
      }
    });
  };

  const flushScrub = () => {
    if (pendingScrubIndex == null) return;
    const idx = pendingScrubIndex;
    pendingScrubIndex = null;
    scrubScheduled = false;
    const start = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    controller.scrubTo?.(idx);
    setActiveButton(pauseBtn);
    if (typeof onScrub === 'function') {
      onScrub(idx);
    }
    if (isDebug()) {
      const end = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
      try {
        console.log('[timeline][scrub-flush]', {
          index: idx,
          durationMs: Math.round((end - start) * 1000) / 1000,
          controllerMode: controller.getMode?.(),
        });
      } catch { }
    }
  };

  slider.addEventListener('input', () => {
    const value = Number(slider.value) || 0;
    updateLabel(value);
    pendingScrubIndex = value;
    if (isDebug()) {
      const now = (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
      const delta = now - lastInputTs;
      lastInputTs = now;
      try { console.log('[timeline][input]', { value, deltaMs: Math.round(delta * 100) / 100 }); } catch { }
    }
    scheduleScrub();
  });
  slider.addEventListener('change', () => {
    if (isDebug()) {
      try { console.log('[timeline][input-change]', { value: Number(slider.value) || 0 }); } catch { }
    }
    flushScrub();
  });

  playBtn.addEventListener('click', () => {
    setActiveButton(playBtn);
    onPlay?.({ index: Number(slider.value) || 0 });
  });

  pauseBtn.addEventListener('click', () => {
    setActiveButton(pauseBtn);
    onPause?.();
  });

  liveBtn.addEventListener('click', () => {
    setActiveButton(liveBtn);
    onLive?.();
  });

  controller.on?.('buffer', (snapshot) => {
    bufferSnapshot = snapshot;
    if (!snapshot || snapshot.size === 0) {
      slider.disabled = true;
      slider.max = '0';
      slider.value = '0';
      label.textContent = 'Frame —';
      return;
    }
    slider.disabled = snapshot.size <= 1;
    slider.max = String(Math.max(0, snapshot.size - 1));
    while (ticks.firstChild) ticks.removeChild(ticks.firstChild);
    const liveIndex = snapshot.liveIndex;
    for (let i = 0; i < snapshot.size; i++) {
      const opt = document.createElement('option');
      opt.value = String(i);
      const rel = i - liveIndex;
      opt.label = i === liveIndex ? '0 (live)' : String(rel);
      ticks.appendChild(opt);
    }
    if (controller.getMode?.() === 'live') {
      slider.value = String(snapshot.liveIndex);
    }
    drawHeat(snapshot);
    updateLabel(Number(slider.value) || 0);
  });

  controller.on?.('index', ({ index }) => {
    if (index == null || Number.isNaN(index)) return;
    slider.value = String(index);
    updateLabel(index);
  });

  controller.on?.('mode', ({ mode }) => {
    latestSnapshot = controller.debug?.();
    if (mode === 'live') {
      setActiveButton(liveBtn);
      updateButtonStates('live');
    } else if (mode === 'paused') {
      setActiveButton(pauseBtn);
      updateButtonStates('paused');
    } else if (mode === 'playback') {
      setActiveButton(playBtn);
      updateButtonStates('playback');
    } else if (mode === 'scrub') {
      setActiveButton(pauseBtn);
      updateButtonStates('paused');
    }
  });

  updateButtonStates(controller.getMode?.() || 'live');

  return {
    element: rail,
    destroy: () => {
      if (collapseTimer) clearTimeout(collapseTimer);
      rail.remove();
    },
    setActiveButton,
    updateLabel,
  };
}

export default { installTimelineOverlay };
