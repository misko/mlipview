// Timeline UI: hover reveal dock with discrete slider + controls.
// Responsible for rendering and DOM interactions only; actual playback handled by caller.

const STYLE_ID = 'timelineStyles';

const MODE_LIVE = 'live';
const MODE_PAUSED = 'paused';
const MODE_PLAYING = 'playing';

function injectStyles() {
  if (typeof document === 'undefined') return;
  if (document.getElementById(STYLE_ID)) return;
  const style = document.createElement('style');
  style.id = STYLE_ID;
  style.textContent = `
    #timelineDock {
      position: absolute;
      left: 0;
      right: 0;
      bottom: 0;
      pointer-events: none;
      z-index: 30;
    }
    #timelineDock[data-visible="true"] .timeline-panel {
      transform: translateY(0);
      opacity: 1;
      pointer-events: auto;
    }
    #timelineDock[data-visible="true"] .timeline-hitbox {
      opacity: 0;
    }
    #timelineDock .timeline-hitbox {
      position: absolute;
      left: 0;
      right: 0;
      bottom: 0;
      height: 18px;
      background: transparent;
      pointer-events: auto;
      cursor: ns-resize;
      transition: opacity 0.2s ease;
    }
    #timelineDock .timeline-panel {
      pointer-events: auto;
      margin: 0 auto;
      max-width: 720px;
      background: rgba(20, 24, 34, 0.9);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 10px 10px 0 0;
      padding: 10px 16px 14px;
      color: #e4ecff;
      font: 12px/1.4 system-ui, sans-serif;
      box-shadow: 0 -6px 16px rgba(0, 0, 0, 0.35);
      transform: translateY(110%);
      opacity: 0;
      transition: transform 0.22s cubic-bezier(0.4, 0, 0.2, 1),
                  opacity 0.22s ease;
    }
    #timelineDock .timeline-controls {
      display: flex;
      align-items: center;
      gap: 10px;
      margin-bottom: 6px;
    }
    #timelineDock .timeline-button {
      background: #2a3446;
      color: inherit;
      border: 1px solid rgba(255, 255, 255, 0.18);
      border-radius: 20px;
      padding: 4px 12px;
      font-size: 12px;
      font-weight: 600;
      cursor: pointer;
      transition: background 0.18s ease, border-color 0.18s ease, transform 0.12s ease;
    }
    #timelineDock .timeline-button:hover {
      background: #33425d;
      border-color: rgba(255, 255, 255, 0.25);
    }
    #timelineDock .timeline-button[data-active="true"] {
      background: #4f8bff;
      border-color: #86aefc;
      color: #0a1421;
    }
    #timelineDock .timeline-button[data-disabled="true"] {
      opacity: 0.5;
      cursor: not-allowed;
    }
    #timelineDock .timeline-slider-row {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    #timelineDock .timeline-slider {
      width: 100%;
      accent-color: #4f8bff;
    }
    #timelineDock .timeline-ticks {
      position: relative;
      width: 100%;
      height: 12px;
      display: flex;
      align-items: flex-start;
      pointer-events: none;
    }
    #timelineDock .timeline-tick {
      position: relative;
      height: 12px;
      width: 1px;
      background: rgba(255, 255, 255, 0.35);
      flex: 0 0 auto;
    }
    #timelineDock .timeline-tick[data-major="true"] {
      height: 14px;
      background: rgba(255, 255, 255, 0.55);
    }
    #timelineDock .timeline-tick-label {
      position: absolute;
      top: 16px;
      left: 50%;
      transform: translateX(-50%);
      font-size: 10px;
      color: rgba(228, 236, 255, 0.7);
      white-space: nowrap;
    }
    #timelineDock[data-mode="playing"] .timeline-button[data-action="play"],
    #timelineDock[data-mode="paused"] .timeline-button[data-action="pause"],
    #timelineDock[data-mode="live"] .timeline-button[data-action="live"] {
      background: #4f8bff;
      border-color: #86aefc;
      color: #0a1421;
    }
    #timelineDock .timeline-status {
      margin-left: auto;
      font-size: 11px;
      opacity: 0.8;
    }
  `;
  document.head.appendChild(style);
}

function clampOffset(offset, min) {
  if (!Number.isFinite(offset)) return null;
  const o = Math.floor(offset);
  if (o > -1) return -1;
  if (min != null && o < min) return min;
  return o;
}

function pickTickSpacing(capacity) {
  if (capacity >= 400) return 25;
  if (capacity >= 200) return 20;
  if (capacity >= 100) return 10;
  if (capacity >= 50) return 5;
  return 2;
}

export function installTimeline({
  host,
  capacity = 500,
  getOffsets,
  getActiveOffset,
  onRequestOffset,
  onRequestPlay,
  onRequestPause,
  onRequestLive,
}) {
  if (typeof document === 'undefined') {
    return {
      refresh() {},
      setMode() {},
      setActiveOffset() {},
      getState() {
        return { mode: MODE_LIVE, offset: -1 };
      },
    };
  }

  injectStyles();

  const parent =
    host || document.getElementById('app') || document.body || document.documentElement;

  const dock = document.createElement('div');
  dock.id = 'timelineDock';
  dock.dataset.mode = MODE_LIVE;
  dock.dataset.visible = 'false';

  const hitbox = document.createElement('div');
  hitbox.className = 'timeline-hitbox';
  hitbox.dataset.testid = 'timeline-hitbox';

  const panel = document.createElement('div');
  panel.className = 'timeline-panel';
  panel.dataset.testid = 'timeline-panel';

  const controlsRow = document.createElement('div');
  controlsRow.className = 'timeline-controls';

  const playBtn = document.createElement('button');
  playBtn.className = 'timeline-button';
  playBtn.dataset.action = 'play';
  playBtn.dataset.testid = 'timeline-play';
  playBtn.textContent = 'Play';

  const pauseBtn = document.createElement('button');
  pauseBtn.className = 'timeline-button';
  pauseBtn.dataset.action = 'pause';
  pauseBtn.dataset.testid = 'timeline-pause';
  pauseBtn.textContent = 'Pause';

  const liveBtn = document.createElement('button');
  liveBtn.className = 'timeline-button';
  liveBtn.dataset.action = 'live';
  liveBtn.dataset.testid = 'timeline-live';
  liveBtn.textContent = 'Live';

  const statusLabel = document.createElement('span');
  statusLabel.className = 'timeline-status';
  statusLabel.dataset.testid = 'timeline-status';
  statusLabel.textContent = 'Live';

  controlsRow.append(playBtn, pauseBtn, liveBtn, statusLabel);

  const sliderRow = document.createElement('div');
  sliderRow.className = 'timeline-slider-row';

  const slider = document.createElement('input');
  slider.type = 'range';
  slider.min = String(Math.min(-1, -capacity));
  slider.max = '-1';
  slider.step = '1';
  slider.value = '-1';
  slider.disabled = true;
  slider.className = 'timeline-slider';
  slider.dataset.testid = 'timeline-slider';

  const ticks = document.createElement('div');
  ticks.className = 'timeline-ticks';
  ticks.dataset.testid = 'timeline-ticks';

  sliderRow.append(slider, ticks);

  panel.append(controlsRow, sliderRow);
  dock.append(hitbox, panel);
  parent.appendChild(dock);

  let mode = MODE_LIVE;
  let activeOffset = -1;
  let minOffset = -capacity;

  const setVisible = (visible) => {
    dock.dataset.visible = visible ? 'true' : 'false';
  };

  hitbox.addEventListener('mouseenter', () => setVisible(true));
  dock.addEventListener('mouseleave', () => {
    if (mode === MODE_LIVE) setVisible(false);
  });

  function setMode(nextMode) {
    if (!nextMode) return;
    mode = nextMode;
    dock.dataset.mode = mode;
    switch (mode) {
      case MODE_PLAYING:
        statusLabel.textContent = `Playing (${activeOffset})`;
        break;
      case MODE_PAUSED:
        statusLabel.textContent = `Paused (${activeOffset})`;
        break;
      default:
        statusLabel.textContent = 'Live';
        break;
    }
    if (mode === MODE_LIVE) {
      setVisible(false);
    } else {
      setVisible(true);
    }
  }

  function setActiveOffset(offset) {
    if (!Number.isFinite(offset)) return;
    activeOffset = clampOffset(offset, minOffset) ?? activeOffset;
    slider.value = String(activeOffset);
  }

  function rebuildTicks(offsets) {
    ticks.innerHTML = '';
    if (!offsets || offsets.length === 0) return;
    const total = offsets.length;
    const cap = Math.min(capacity, total ? Math.abs(offsets[offsets.length - 1]) : capacity);
    const spacing = pickTickSpacing(capacity);
    const width = ticks.clientWidth || panel.clientWidth || 1;

    const tickOffsets = new Set([offsets[0], offsets[offsets.length - 1], -1]);
    for (let i = spacing; i <= capacity; i += spacing) {
      tickOffsets.add(-i);
    }

    const sorted = Array.from(tickOffsets)
      .filter((o) => o <= -1 && o >= offsets[offsets.length - 1])
      .sort((a, b) => a - b);

    const max = -1;
    const min = offsets[offsets.length - 1] || -1;
    const range = Math.abs(min - max) || 1;

    for (const off of sorted) {
      const tick = document.createElement('div');
      tick.className = 'timeline-tick';
      tick.dataset.offset = String(off);
      const major = off === -1 || off === min || off === Math.floor(min / 2) || off % spacing === 0;
      if (major) tick.dataset.major = 'true';

      if (width) {
        const pct = (Math.abs(off - min) / range) * 100;
        tick.style.left = `${pct}%`;
        tick.style.position = 'absolute';
      }

      if (major) {
        const label = document.createElement('div');
        label.className = 'timeline-tick-label';
        label.textContent = `${off}`;
        tick.appendChild(label);
      }
      ticks.appendChild(tick);
    }
  }

  function refresh() {
    const offsets = (typeof getOffsets === 'function' ? getOffsets() : []) || [];
    if (!offsets.length) {
      slider.disabled = true;
      slider.value = '-1';
      activeOffset = -1;
      minOffset = -capacity;
      rebuildTicks([]);
      return;
    }
    slider.disabled = false;
    minOffset = offsets[offsets.length - 1];
    slider.min = String(minOffset);
    slider.max = '-1';
    const current = typeof getActiveOffset === 'function' ? getActiveOffset() : null;
    if (Number.isFinite(current)) {
      setActiveOffset(current);
    } else if (!Number.isFinite(activeOffset)) {
      setActiveOffset(-1);
    }
    rebuildTicks(offsets);
  }

  function handleOffsetRequest(offset, source) {
    const clamped = clampOffset(offset, minOffset);
    if (clamped == null) return;
    setActiveOffset(clamped);
    if (typeof onRequestOffset === 'function') {
      onRequestOffset(clamped, { source });
    }
  }

  slider.addEventListener('change', (evt) => {
    const value = Number(evt.target.value);
    handleOffsetRequest(value, { kind: 'slider' });
  });

  slider.addEventListener('pointerup', (evt) => {
    const value = Number(slider.value);
    handleOffsetRequest(value, { kind: 'slider-pointerup' });
  });

  playBtn.addEventListener('click', () => {
    if (typeof onRequestPlay === 'function') {
      onRequestPlay(activeOffset);
    }
  });

  pauseBtn.addEventListener('click', () => {
    if (typeof onRequestPause === 'function') {
      onRequestPause(activeOffset);
    }
  });

  liveBtn.addEventListener('click', () => {
    if (typeof onRequestLive === 'function') {
      onRequestLive();
    }
  });

  const api = {
    refresh,
    setMode,
    setActiveOffset,
    getState() {
      return { mode, offset: activeOffset, visible: dock.dataset.visible === 'true' };
    },
    elements: {
      dock,
      slider,
      playBtn,
      pauseBtn,
      liveBtn,
    },
  };

  return api;
}

export default { installTimeline };
