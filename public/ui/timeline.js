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
      transform: translateY(-14px);
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
      display: flex;
      justify-content: center;
      align-items: flex-end;
      padding-bottom: 2px;
    }
    #timelineDock .timeline-hitbox::before {
      content: '';
      position: absolute;
      bottom: -18px;
      width: 140px;
      height: 70px;
      border-radius: 140px 140px 0 0;
      background: rgba(32, 40, 52, 0.55);
      border: 1px solid rgba(120, 130, 150, 0.2);
      filter: drop-shadow(0 -2px 4px rgba(0,0,0,0.35));
      transform: translateY(14px);
      transition: transform 0.2s ease, opacity 0.2s ease;
    }
    #timelineDock .timeline-hitbox::after {
      content: '\\23F1';
      position: relative;
      bottom: -16px;
      font-size: 36px;
      font-weight: 700;
      color: rgba(228, 236, 255, 0.65);
      text-shadow: 0 -1px 2px rgba(0, 0, 0, 0.45);
      transition: transform 0.2s ease, opacity 0.2s ease;
    }
    #timelineDock[data-visible="true"] .timeline-hitbox::before,
    #timelineDock[data-visible="true"] .timeline-hitbox::after {
      transform: translateY(24px);
      opacity: 0;
    }
    #timelineDock .timeline-panel {
      pointer-events: auto;
      margin: 0 auto;
      max-width: 720px;
      background: rgba(20, 24, 34, 0.9);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 12px;
      padding: 10px 16px 14px;
      color: #e4ecff;
      font: 12px/1.4 system-ui, sans-serif;
      box-shadow: 0 -6px 16px rgba(0, 0, 0, 0.35);
      transform: translateY(110%);
      opacity: 0;
      transition: transform 0.22s cubic-bezier(0.4, 0, 0.2, 1),
                  opacity 0.22s ease;
    }
    #timelineDock .timeline-main-row {
      display: flex;
      align-items: center;
      gap: 12px;
    }
    #timelineDock .timeline-controls {
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 0;
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
      flex: 1 1 auto;
      display: flex;
      flex-direction: column;
    }
    #timelineDock .timeline-slider {
      width: 100%;
      accent-color: #4f8bff;
    }
    #timelineDock[data-mode="playing"] .timeline-button[data-action="play"],
    #timelineDock[data-mode="paused"] .timeline-button[data-action="pause"],
    #timelineDock[data-mode="live"] .timeline-button[data-action="live"] {
      background: #4f8bff;
      border-color: #86aefc;
      color: #0a1421;
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
      refresh() { },
      setMode() { },
      setActiveOffset() { },
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
  playBtn.textContent = '▶';

  const pauseBtn = document.createElement('button');
  pauseBtn.className = 'timeline-button';
  pauseBtn.dataset.action = 'pause';
  pauseBtn.dataset.testid = 'timeline-pause';
  pauseBtn.textContent = '⏸';

  const liveBtn = document.createElement('button');
  liveBtn.className = 'timeline-button';
  liveBtn.dataset.action = 'live';
  liveBtn.dataset.testid = 'timeline-live';
  liveBtn.textContent = '⏭';

  controlsRow.append(playBtn, pauseBtn, liveBtn);

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

  sliderRow.append(slider);

  const mainRow = document.createElement('div');
  mainRow.className = 'timeline-main-row';
  mainRow.append(sliderRow, controlsRow);

  panel.append(mainRow);
  dock.append(hitbox, panel);
  parent.appendChild(dock);

  let mode = MODE_LIVE;
  let activeOffset = -1;
  let minOffset = -capacity;
  let pointerScrubActive = false;
  let offsetList = [];
  let pendingPointerOffset = null;
  let ignoreNextChange = false;

  const setVisible = (visible) => {
    dock.dataset.visible = visible ? 'true' : 'false';
  };

  hitbox.addEventListener('mouseenter', () => setVisible(true));
  dock.addEventListener('mouseleave', () => {
    if (mode === MODE_LIVE) setVisible(false);
  });

  const setDisabled = (btn, disabled) => {
    const on = !!disabled;
    btn.disabled = on;
    if (on) btn.setAttribute('disabled', '');
    else btn.removeAttribute('disabled');
    btn.dataset.disabled = on ? 'true' : 'false';
  };

  function updateButtons() {
    const isLive = mode === MODE_LIVE;
    const isPaused = mode === MODE_PAUSED;
    const isPlaying = mode === MODE_PLAYING;

    setDisabled(playBtn, isLive || isPlaying);
    setDisabled(pauseBtn, isPaused);
    setDisabled(liveBtn, isLive);
  }

  function setMode(nextMode) {
    if (!nextMode) return;
    mode = nextMode;
    dock.dataset.mode = mode;
    setVisible(mode !== MODE_LIVE);
    updateButtons();
  }

  function setActiveOffset(offset) {
    if (!Number.isFinite(offset)) return;
    activeOffset = clampOffset(offset, minOffset) ?? activeOffset;
    slider.value = String(activeOffset);
  }

  function refresh() {
    const offsets = (typeof getOffsets === 'function' ? getOffsets() : []) || [];
    if (!offsets.length) {
      slider.disabled = true;
      slider.value = '-1';
      activeOffset = -1;
      minOffset = -capacity;
      offsetList = [];
      return;
    }
    slider.disabled = false;
    minOffset = offsets[offsets.length - 1];
    slider.min = String(minOffset);
    slider.max = '-1';
    offsetList = offsets.slice();
    const current = typeof getActiveOffset === 'function' ? getActiveOffset() : null;
    if (Number.isFinite(current)) {
      setActiveOffset(current);
    } else if (!Number.isFinite(activeOffset)) {
      setActiveOffset(-1);
    }
  }

  function handleOffsetRequest(offset, source) {
    const clamped = clampOffset(offset, minOffset);
    if (clamped == null) return;
    setActiveOffset(clamped);
    if (typeof onRequestOffset === 'function') {
      onRequestOffset(clamped, { source });
    }
  }

  function offsetFromClientX(clientX) {
    const rect = slider.getBoundingClientRect();
    if (!rect || rect.width === 0) return Number(slider.value);
    if (!offsetList.length) return Number(slider.value);
    const clampedX = Math.min(rect.right, Math.max(rect.left, clientX));
    const ratio = (clampedX - rect.left) / rect.width;
    const idx = Math.round((offsetList.length - 1) * (1 - ratio));
    const boundedIdx = Math.max(0, Math.min(offsetList.length - 1, idx));
    return Number(offsetList[boundedIdx]);
  }

  slider.addEventListener('pointerdown', (evt) => {
    pointerScrubActive = true;
    try { slider.setPointerCapture?.(evt.pointerId); } catch { }
    if (typeof evt.preventDefault === 'function') evt.preventDefault();
    const nextOffset = offsetFromClientX(evt.clientX);
    pendingPointerOffset = nextOffset;
    slider.value = String(nextOffset);
    handleOffsetRequest(nextOffset, { kind: 'slider-pointerdown' });
  });

  slider.addEventListener('pointermove', (evt) => {
    if (!pointerScrubActive) return;
    if (typeof evt.preventDefault === 'function') evt.preventDefault();
    const nextOffset = offsetFromClientX(evt.clientX);
    pendingPointerOffset = nextOffset;
    slider.value = String(nextOffset);
    handleOffsetRequest(nextOffset, { kind: 'slider-pointermove' });
  });

  slider.addEventListener('change', (evt) => {
    if (ignoreNextChange) {
      ignoreNextChange = false;
      return;
    }
    const value = Number(evt.target.value);
    handleOffsetRequest(value, { kind: 'slider' });
  });

  const finishPointerScrub = (evt) => {
    if (!pointerScrubActive) return;
    pointerScrubActive = false;
    try { slider.releasePointerCapture?.(evt.pointerId); } catch { }
    if (typeof evt.preventDefault === 'function') evt.preventDefault();
    const value = offsetFromClientX(evt.clientX ?? evt.pageX ?? 0);
    pendingPointerOffset = null;
    slider.value = String(value);
    handleOffsetRequest(value, { kind: 'slider-pointerup' });
    ignoreNextChange = true;
  };

  slider.addEventListener('pointerup', finishPointerScrub);
  slider.addEventListener('pointercancel', () => {
    pointerScrubActive = false;
  });

  playBtn.addEventListener('click', () => {
    if (playBtn.disabled) return;
    if (typeof onRequestPlay === 'function') {
      onRequestPlay(activeOffset);
    }
  });

  pauseBtn.addEventListener('click', () => {
    if (pauseBtn.disabled) return;
    if (typeof onRequestPause === 'function') {
      onRequestPause(activeOffset);
    }
  });

  liveBtn.addEventListener('click', () => {
    if (liveBtn.disabled) return;
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

  updateButtons();

  return api;
}

export default { installTimeline };
