const STYLE_ID = 'timelineEditorPanelStyles';

function injectStyles() {
  if (typeof document === 'undefined') return;
  if (document.getElementById(STYLE_ID)) return;
  const style = document.createElement('style');
  style.id = STYLE_ID;
  style.textContent = `
    #timelineEditorPanel {
      position: absolute;
      top: 0;
      right: 0;
      bottom: 0;
      width: 360px;
      background: rgba(14, 18, 28, 0.96);
      border-left: 1px solid rgba(110, 130, 170, 0.18);
      box-shadow: -6px 0 18px rgba(0, 0, 0, 0.45);
      color: #e8f0ff;
      font: 12px/1.5 system-ui, sans-serif;
      display: flex;
      flex-direction: column;
      z-index: 50;
      padding: 16px 18px;
      gap: 16px;
      overflow-y: auto;
    }
    #timelineEditorPanel h2 {
      margin: 0 0 8px 0;
      font-size: 14px;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      color: rgba(209, 218, 240, 0.75);
    }
    #timelineEditorPanel .section {
      background: rgba(24, 32, 48, 0.75);
      border: 1px solid rgba(130, 150, 190, 0.18);
      border-radius: 12px;
      padding: 12px 14px;
      display: flex;
      flex-direction: column;
      gap: 10px;
    }
    #timelineEditorPanel label {
      display: flex;
      flex-direction: column;
      gap: 4px;
      font-weight: 600;
      color: rgba(220, 230, 255, 0.78);
    }
    #timelineEditorPanel input,
    #timelineEditorPanel select,
    #timelineEditorPanel textarea {
      background: rgba(12, 18, 32, 0.9);
      border: 1px solid rgba(110, 130, 170, 0.25);
      border-radius: 8px;
      color: #f7fbff;
      padding: 6px 8px;
      font: inherit;
    }
    #timelineEditorPanel textarea {
      min-height: 70px;
      resize: vertical;
    }
    #timelineEditorPanel button {
      appearance: none;
      border: 1px solid rgba(100, 140, 220, 0.4);
      border-radius: 20px;
      padding: 6px 12px;
      font: inherit;
      font-weight: 600;
      background: rgba(74, 114, 255, 0.2);
      color: #dce6ff;
      cursor: pointer;
      transition: background 0.15s ease, border-color 0.15s ease;
    }
    #timelineEditorPanel button:hover {
      background: rgba(74, 114, 255, 0.35);
      border-color: rgba(120, 160, 250, 0.55);
    }
    #timelineEditorPanel .secondary {
      background: rgba(18, 22, 34, 0.7);
      border-color: rgba(90, 110, 150, 0.35);
    }
    #timelineEditorPanel .danger {
      background: rgba(205, 80, 80, 0.18);
      border-color: rgba(220, 120, 120, 0.35);
      color: #f8dede;
    }
    #timelineEditorPanel .row {
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
    }
    #timelineEditorPanel .col {
      flex: 1 1 0;
      min-width: 0;
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    #timelineEditorPanel .list {
      display: flex;
      flex-direction: column;
      gap: 6px;
      max-height: 200px;
      overflow-y: auto;
    }
    #timelineEditorPanel .list-item {
      padding: 8px 10px;
      border-radius: 10px;
      background: rgba(18, 24, 40, 0.75);
      border: 1px solid transparent;
      cursor: pointer;
      transition: background 0.12s ease, border-color 0.12s ease;
    }
    #timelineEditorPanel .list-item[data-active="true"] {
      background: rgba(72, 104, 255, 0.25);
      border-color: rgba(92, 132, 255, 0.6);
    }
    #timelineEditorPanel .list-item-title {
      font-weight: 700;
      color: #f3f6ff;
    }
    #timelineEditorPanel .list-item-sub {
      font-size: 11px;
      color: rgba(200, 214, 240, 0.7);
    }
    #timelineEditorPanel .range-grid {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 6px;
    }
    #timelineEditorPanel .helper-text {
      font-size: 11px;
      color: rgba(190, 205, 230, 0.75);
    }
    #timelineEditorPanel .section-divider {
      height: 1px;
      background: rgba(100, 120, 160, 0.25);
      margin: 4px 0;
    }
    #timelineEditorPanel .action-header {
      display: flex;
      align-items: center;
      justify-content: space-between;
    }
    #timelineEditorPanel .action-header h3 {
      margin: 0;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: rgba(210, 220, 242, 0.75);
    }
    #timelineEditorPanel .action-toggle {
      display: flex;
      align-items: center;
      gap: 6px;
      font-weight: 600;
    }
    #timelineEditorPanel .action-toggle input[type="checkbox"] {
      width: 16px;
      height: 16px;
      accent-color: #5b7bff;
    }
    #timelineEditorPanel .error {
      color: #ff9c9c;
      font-size: 11px;
    }
  `;
  document.head.appendChild(style);
}

function deepClone(value) {
  try {
    return value == null ? value : JSON.parse(JSON.stringify(value));
  } catch {
    return value;
  }
}

function parseNumber(value) {
  const n = Number(value);
  return Number.isFinite(n) ? n : null;
}

function uniqueId(prefix) {
  return `${prefix}-${Date.now().toString(36)}-${Math.floor(Math.random() * 1e6)}`;
}

function ensureArray(value) {
  if (Array.isArray(value)) return value;
  if (value == null) return [];
  return [value];
}

function atomsToString(atoms) {
  return ensureArray(atoms)
    .filter((v) => Number.isInteger(v) && v >= 0)
    .join(', ');
}

function parseAtoms(text) {
  if (typeof text !== 'string') return [];
  return text
    .split(/[, ]+/)
    .map((part) => Number(part.trim()))
    .filter((n) => Number.isInteger(n) && n >= 0);
}

function formatRangeSummary(range, resolver) {
  if (!range) return 'Full timeline';
  const parts = [];
  const start = resolver(range.start || {});
  const end = resolver(range.end || {});
  if (start) parts.push(`${start}`);
  if (end) parts.push(`→ ${end}`);
  return parts.length ? parts.join(' ') : 'Full timeline';
}

function engineMessageToDraft(msg = {}) {
  const draft = {
    id: msg.id || uniqueId('cm'),
    label: typeof msg.label === 'string' && msg.label ? msg.label : (msg.id || ''),
    priority: Number.isFinite(msg.priority) ? msg.priority : 0,
    notes: typeof msg.notes === 'string' ? msg.notes : '',
    range: {
      start: {
        frameId: msg.range?.start?.frameId || '',
        frameIndex: msg.range?.start?.frameIndex != null ? String(msg.range.start.frameIndex) : '',
        offset: msg.range?.start?.offset != null ? String(msg.range.start.offset) : '',
      },
      end: {
        frameId: msg.range?.end?.frameId || '',
        frameIndex: msg.range?.end?.frameIndex != null ? String(msg.range.end.frameIndex) : '',
        offset: msg.range?.end?.offset != null ? String(msg.range.end.offset) : '',
        inclusive: msg.range?.end?.inclusive !== false,
      },
    },
    actions: {
      speed: {
        enabled: false,
        mode: 'fps',
        fps: '',
        multiplier: '',
        transitionMs: '',
        easing: '',
      },
      callout: {
        enabled: false,
        text: '',
        anchorMode: 'atom',
        anchor: {},
        offset: ['', '', ''],
        panelWidth: '',
        panelHeight: '',
        textSize: '',
        background: '',
        color: '',
      },
      opacity: {
        enabled: false,
        focusAtomsText: '',
        includeBonds: 'none',
        focusOpacityAtoms: '',
        focusOpacityBonds: '',
        backgroundAtoms: '',
        backgroundBonds: '',
        transitionMs: '',
      },
    },
  };

  ensureArray(msg.actions).forEach((action) => {
    if (!action || typeof action !== 'object') return;
    if (action.type === 'timeline.playbackSpeed') {
      draft.actions.speed.enabled = true;
      if (Number.isFinite(action.fps)) {
        draft.actions.speed.mode = 'fps';
        draft.actions.speed.fps = String(action.fps);
      } else if (Number.isFinite(action.speedMultiplier)) {
        draft.actions.speed.mode = 'multiplier';
        draft.actions.speed.multiplier = String(action.speedMultiplier);
      }
      if (Number.isFinite(action.transitionMs)) draft.actions.speed.transitionMs = String(action.transitionMs);
      if (typeof action.easing === 'string') draft.actions.speed.easing = action.easing;
    } else if (action.type === 'overlay.callout') {
      draft.actions.callout.enabled = true;
      draft.actions.callout.text = typeof action.text === 'string' ? action.text : '';
      draft.actions.callout.anchor = deepClone(action.anchor || {});
      if (typeof draft.actions.callout.anchor.mode === 'string') {
        draft.actions.callout.anchorMode = draft.actions.callout.anchor.mode;
      }
      if (Array.isArray(action.offset)) {
        draft.actions.callout.offset = action.offset.map((v) => String(v ?? ''));
      }
      if (action.panelSize) {
        if (Number.isFinite(action.panelSize.width)) draft.actions.callout.panelWidth = String(action.panelSize.width);
        if (Number.isFinite(action.panelSize.height)) draft.actions.callout.panelHeight = String(action.panelSize.height);
      }
      if (Number.isFinite(action.textSize)) draft.actions.callout.textSize = String(action.textSize);
      if (action.style) {
        if (typeof action.style.background === 'string') draft.actions.callout.background = action.style.background;
        if (typeof action.style.color === 'string') draft.actions.callout.color = action.style.color;
      }
    } else if (action.type === 'visual.opacityFocus') {
      draft.actions.opacity.enabled = true;
      draft.actions.opacity.focusAtomsText = atomsToString(action.focus?.atoms);
      draft.actions.opacity.includeBonds = action.focus?.includeBonds || 'none';
      if (Number.isFinite(action.focusOpacity?.atoms)) draft.actions.opacity.focusOpacityAtoms = String(action.focusOpacity.atoms);
      if (Number.isFinite(action.focusOpacity?.bonds)) draft.actions.opacity.focusOpacityBonds = String(action.focusOpacity.bonds);
      if (Number.isFinite(action.backgroundOpacity?.atoms)) draft.actions.opacity.backgroundAtoms = String(action.backgroundOpacity.atoms);
      if (Number.isFinite(action.backgroundOpacity?.bonds)) draft.actions.opacity.backgroundBonds = String(action.backgroundOpacity.bonds);
      if (Number.isFinite(action.transitionMs)) draft.actions.opacity.transitionMs = String(action.transitionMs);
    }
  });

  return draft;
}

function buildRangePayload(rangeDraft = {}) {
  const start = {};
  const end = {};
  const s = rangeDraft.start || {};
  const e = rangeDraft.end || {};
  if (s.frameId) start.frameId = s.frameId.trim();
  if (s.frameIndex !== '' && s.frameIndex != null) {
    const v = parseNumber(s.frameIndex);
    if (v != null) start.frameIndex = v;
  }
  if (s.offset !== '' && s.offset != null) {
    const v = parseNumber(s.offset);
    if (v != null) start.offset = v;
  }
  if (e.frameId) end.frameId = e.frameId.trim();
  if (e.frameIndex !== '' && e.frameIndex != null) {
    const v = parseNumber(e.frameIndex);
    if (v != null) end.frameIndex = v;
  }
  if (e.offset !== '' && e.offset != null) {
    const v = parseNumber(e.offset);
    if (v != null) end.offset = v;
  }
  if (typeof e.inclusive === 'boolean' && e.frameId) {
    end.inclusive = e.inclusive;
  } else if (typeof e.inclusive === 'boolean' && (e.frameIndex !== '' || e.offset !== '')) {
    end.inclusive = e.inclusive;
  }
  const hasStart = Object.keys(start).length > 0;
  const hasEnd = Object.keys(end).length > 0;
  const range = {};
  if (hasStart) range.start = start;
  if (hasEnd) range.end = end;
  return hasStart || hasEnd ? range : null;
}

function buildActionsFromDraft(draftActions = {}) {
  const actions = [];
  if (draftActions.speed?.enabled) {
    const payload = { type: 'timeline.playbackSpeed' };
    if (draftActions.speed.mode === 'multiplier') {
      const mult = parseNumber(draftActions.speed.multiplier);
      if (mult != null && mult > 0) payload.speedMultiplier = mult;
    } else {
      const fps = parseNumber(draftActions.speed.fps);
      if (fps != null && fps > 0) payload.fps = fps;
    }
    if (Number.isFinite(parseNumber(draftActions.speed.transitionMs))) {
      const t = parseNumber(draftActions.speed.transitionMs);
      if (t != null && t >= 0) payload.transitionMs = t;
    }
    if (draftActions.speed.easing) payload.easing = draftActions.speed.easing;
    if (payload.fps || payload.speedMultiplier) actions.push(payload);
  }

  if (draftActions.callout?.enabled) {
    const payload = {
      type: 'overlay.callout',
      text: draftActions.callout.text || '',
      anchor: { ...(draftActions.callout.anchor || {}) },
    };
    const mode = draftActions.callout.anchorMode || payload.anchor.mode || 'atom';
    payload.anchor.mode = mode;
    if (mode === 'atom') {
      const atoms = parseAtoms(draftActions.callout.anchor?.atoms || draftActions.callout.anchorAtomsText);
      payload.anchor.atoms = atoms;
    } else if (mode === 'bond') {
      const atoms = parseAtoms(draftActions.callout.anchor?.atoms);
      payload.anchor.atoms = atoms.slice(0, 2);
    } else if (mode === 'world' && Array.isArray(draftActions.callout.worldPosition)) {
      payload.anchor.position = draftActions.callout.worldPosition.map((v) => Number(v) || 0);
    } else if (mode === 'world' && Array.isArray(draftActions.callout.offset)) {
      payload.anchor.position = draftActions.callout.offset.map((v) => Number(v) || 0);
    }
    if (Array.isArray(draftActions.callout.offset) && draftActions.callout.offset.some((v) => v !== '' && v != null)) {
      payload.offset = draftActions.callout.offset.map((v) => Number(v) || 0);
    }
    const panelSize = {};
    const width = parseNumber(draftActions.callout.panelWidth);
    const height = parseNumber(draftActions.callout.panelHeight);
    if (width != null && width > 0) panelSize.width = width;
    if (height != null && height > 0) panelSize.height = height;
    if (Object.keys(panelSize).length) payload.panelSize = panelSize;
    const textSize = parseNumber(draftActions.callout.textSize);
    if (textSize != null && textSize > 0) payload.textSize = textSize;
    const style = {};
    if (draftActions.callout.background) style.background = draftActions.callout.background;
    if (draftActions.callout.color) style.color = draftActions.callout.color;
    if (Object.keys(style).length) payload.style = style;
    actions.push(payload);
  }

  if (draftActions.opacity?.enabled) {
    const payload = { type: 'visual.opacityFocus' };
    const focusAtoms = parseAtoms(draftActions.opacity.focusAtomsText);
    payload.focus = { atoms: focusAtoms };
    if (draftActions.opacity.includeBonds) payload.focus.includeBonds = draftActions.opacity.includeBonds;
    const focusOpacity = {};
    const foa = parseNumber(draftActions.opacity.focusOpacityAtoms);
    if (foa != null && foa >= 0) focusOpacity.atoms = foa;
    const fob = parseNumber(draftActions.opacity.focusOpacityBonds);
    if (fob != null && fob >= 0) focusOpacity.bonds = fob;
    if (Object.keys(focusOpacity).length) payload.focusOpacity = focusOpacity;
    const backgroundOpacity = {};
    const boa = parseNumber(draftActions.opacity.backgroundAtoms);
    if (boa != null && boa >= 0) backgroundOpacity.atoms = boa;
    const bob = parseNumber(draftActions.opacity.backgroundBonds);
    if (bob != null && bob >= 0) backgroundOpacity.bonds = bob;
    if (Object.keys(backgroundOpacity).length) payload.backgroundOpacity = backgroundOpacity;
    const transitionMs = parseNumber(draftActions.opacity.transitionMs);
    if (transitionMs != null && transitionMs >= 0) payload.transitionMs = transitionMs;
    actions.push(payload);
  }
  return actions;
}

export function createTimelineEditorPanel(options = {}) {
  if (typeof document === 'undefined') {
    return {
      destroy() {},
      refresh() {},
      getState() { return {}; },
      select() {},
      getCurrentDraft() { return null; },
    };
  }

  injectStyles();

  const host = options.attachTo || document.getElementById('app') || document.body;
  const panel = document.createElement('div');
  panel.id = 'timelineEditorPanel';
  host.appendChild(panel);

  const playbackSection = document.createElement('div');
  playbackSection.className = 'section';
  const playbackTitle = document.createElement('h2');
  playbackTitle.textContent = 'Playback';
  playbackSection.appendChild(playbackTitle);
  panel.appendChild(playbackSection);

  const messageSection = document.createElement('div');
  messageSection.className = 'section';
  const messageTitle = document.createElement('h2');
  messageTitle.textContent = 'Control Messages';
  messageSection.appendChild(messageTitle);
  panel.appendChild(messageSection);

  const detailSection = document.createElement('div');
  detailSection.className = 'section';
  const detailTitle = document.createElement('h2');
  detailTitle.textContent = 'Message Editor';
  detailSection.appendChild(detailTitle);
  panel.appendChild(detailSection);

  const errorBox = document.createElement('div');
  errorBox.className = 'error';
  panel.insertBefore(errorBox, panel.firstChild);

  function setError(message) {
    if (!message) {
      errorBox.textContent = '';
      errorBox.style.display = 'none';
      return;
    }
    errorBox.textContent = message;
    errorBox.style.display = 'block';
  }

  let playbackDraft = options.getPlaybackConfig ? deepClone(options.getPlaybackConfig()) : {
    defaultFps: 20,
    autoPlay: false,
    loop: false,
    loopRange: null,
    startFrame: null,
  };

  const messageMap = new Map();
  let orderedIds = [];
  let selectedId = null;

  function loadMessages() {
    const current = options.getMessages ? options.getMessages() : [];
    messageMap.clear();
    orderedIds = [];
    ensureArray(current).forEach((msg) => {
      const draft = engineMessageToDraft(msg);
      messageMap.set(draft.id, draft);
      orderedIds.push(draft.id);
    });
    if (!selectedId && orderedIds.length) {
      selectedId = orderedIds[0];
    } else if (selectedId && !messageMap.has(selectedId)) {
      selectedId = orderedIds[0] || null;
    }
  }

  loadMessages();

  const playbackForm = document.createElement('div');
  playbackForm.className = 'col';
  playbackSection.appendChild(playbackForm);

  function updatePlaybackInput(name, value) {
    if (!playbackDraft) playbackDraft = {};
    playbackDraft[name] = value;
  }

  function renderPlayback() {
    playbackForm.innerHTML = '';

    const fpsLabel = document.createElement('label');
    fpsLabel.textContent = 'Default FPS';
    const fpsInput = document.createElement('input');
    fpsInput.type = 'number';
    fpsInput.min = '1';
    fpsInput.step = '1';
    fpsInput.value = String(playbackDraft.defaultFps || 20);
    fpsInput.addEventListener('input', () => {
      const val = parseNumber(fpsInput.value);
      if (val != null && val > 0) {
        updatePlaybackInput('defaultFps', val);
      }
    });
    fpsLabel.appendChild(fpsInput);
    playbackForm.appendChild(fpsLabel);

    const togglesRow = document.createElement('div');
    togglesRow.className = 'row';
    const autoLabel = document.createElement('label');
    autoLabel.className = 'action-toggle';
    const autoInput = document.createElement('input');
    autoInput.type = 'checkbox';
    autoInput.checked = !!playbackDraft.autoPlay;
    autoInput.addEventListener('change', () => updatePlaybackInput('autoPlay', autoInput.checked));
    autoLabel.appendChild(autoInput);
    autoLabel.appendChild(document.createTextNode('Auto Play'));

    const loopLabel = document.createElement('label');
    loopLabel.className = 'action-toggle';
    const loopInput = document.createElement('input');
    loopInput.type = 'checkbox';
    loopInput.checked = !!playbackDraft.loop;
    loopInput.addEventListener('change', () => {
      updatePlaybackInput('loop', loopInput.checked);
    });
    loopLabel.appendChild(loopInput);
    loopLabel.appendChild(document.createTextNode('Loop'));

    togglesRow.appendChild(autoLabel);
    togglesRow.appendChild(loopLabel);
    playbackForm.appendChild(togglesRow);

    const startLabel = document.createElement('label');
    startLabel.textContent = 'Start Frame (optional)';
    const startRow = document.createElement('div');
    startRow.className = 'row';
    const startInput = document.createElement('input');
    startInput.placeholder = 'frame-00042';
    startInput.value = playbackDraft.startFrame?.frameId || '';
    startInput.addEventListener('input', () => {
      if (!playbackDraft.startFrame) playbackDraft.startFrame = {};
      playbackDraft.startFrame.frameId = startInput.value.trim();
    });
    const startOffsetInput = document.createElement('input');
    startOffsetInput.type = 'number';
    startOffsetInput.placeholder = 'offset';
    startOffsetInput.value = playbackDraft.startFrame?.offset != null ? String(playbackDraft.startFrame.offset) : '';
    startOffsetInput.addEventListener('input', () => {
      if (!playbackDraft.startFrame) playbackDraft.startFrame = {};
      const val = parseNumber(startOffsetInput.value);
      if (val == null) delete playbackDraft.startFrame.offset;
      else playbackDraft.startFrame.offset = val;
    });
    const startBtn = document.createElement('button');
    startBtn.className = 'secondary';
    startBtn.type = 'button';
    startBtn.textContent = 'Use current';
    startBtn.addEventListener('click', () => {
      if (!playbackDraft.startFrame) playbackDraft.startFrame = {};
      const offset = options.getCurrentOffset ? options.getCurrentOffset() : -1;
      const rid = options.offsetToFrameId ? options.offsetToFrameId(offset) : null;
      if (rid) {
        startInput.value = rid;
        playbackDraft.startFrame.frameId = rid;
      }
      if (Number.isFinite(offset)) {
        startOffsetInput.value = String(offset);
        playbackDraft.startFrame.offset = offset;
      }
    });
    startRow.appendChild(startInput);
    startRow.appendChild(startOffsetInput);
    startRow.appendChild(startBtn);
    startLabel.appendChild(startRow);
    playbackForm.appendChild(startLabel);

    const loopRangeLabel = document.createElement('label');
    loopRangeLabel.textContent = 'Loop Range (optional)';
    const loopRangeGrid = document.createElement('div');
    loopRangeGrid.className = 'range-grid';
    const loopStart = document.createElement('input');
    loopStart.placeholder = 'start frameId';
    loopStart.value = playbackDraft.loopRange?.start?.frameId || playbackDraft.loopRange?.startFrameId || '';
    loopStart.addEventListener('input', () => {
      if (!playbackDraft.loopRange) playbackDraft.loopRange = {};
      playbackDraft.loopRange.start = playbackDraft.loopRange.start || {};
      playbackDraft.loopRange.start.frameId = loopStart.value.trim();
    });
    const loopEnd = document.createElement('input');
    loopEnd.placeholder = 'end frameId';
    loopEnd.value = playbackDraft.loopRange?.end?.frameId || playbackDraft.loopRange?.endFrameId || '';
    loopEnd.addEventListener('input', () => {
      if (!playbackDraft.loopRange) playbackDraft.loopRange = {};
      playbackDraft.loopRange.end = playbackDraft.loopRange.end || {};
      playbackDraft.loopRange.end.frameId = loopEnd.value.trim();
    });
    const loopBtn = document.createElement('button');
    loopBtn.className = 'secondary';
    loopBtn.type = 'button';
    loopBtn.textContent = 'Use current';
    loopBtn.addEventListener('click', () => {
      if (!playbackDraft.loopRange) playbackDraft.loopRange = {};
      const offset = options.getCurrentOffset ? options.getCurrentOffset() : -1;
      const rid = options.offsetToFrameId ? options.offsetToFrameId(offset) : null;
      if (!playbackDraft.loopRange.start) playbackDraft.loopRange.start = {};
      if (rid) {
        loopStart.value = rid;
        playbackDraft.loopRange.start.frameId = rid;
      }
      if (!playbackDraft.loopRange.end) playbackDraft.loopRange.end = {};
      if (rid) {
        loopEnd.value = rid;
        playbackDraft.loopRange.end.frameId = rid;
      }
    });
    loopRangeGrid.appendChild(loopStart);
    loopRangeGrid.appendChild(loopEnd);
    loopRangeGrid.appendChild(loopBtn);
    loopRangeLabel.appendChild(loopRangeGrid);
    playbackForm.appendChild(loopRangeLabel);
  }

  renderPlayback();

  const listContainer = document.createElement('div');
  listContainer.className = 'list';
  messageSection.appendChild(listContainer);

  const listButtonsRow = document.createElement('div');
  listButtonsRow.className = 'row';
  const addButton = document.createElement('button');
  addButton.type = 'button';
  addButton.textContent = 'Add Message';
  addButton.addEventListener('click', () => {
    const draft = engineMessageToDraft({});
    draft.label = `Control ${orderedIds.length + 1}`;
    messageMap.set(draft.id, draft);
    orderedIds.push(draft.id);
    selectedId = draft.id;
    renderList();
    renderDetail();
  });
  const removeButton = document.createElement('button');
  removeButton.type = 'button';
  removeButton.className = 'danger';
  removeButton.textContent = 'Remove';
  removeButton.addEventListener('click', () => {
    if (!selectedId) return;
    const idx = orderedIds.indexOf(selectedId);
    messageMap.delete(selectedId);
    orderedIds = orderedIds.filter((id) => id !== selectedId);
    selectedId = orderedIds[Math.max(0, idx - 1)] || orderedIds[0] || null;
    renderList();
    renderDetail();
  });
  listButtonsRow.appendChild(addButton);
  listButtonsRow.appendChild(removeButton);
  messageSection.appendChild(listButtonsRow);

  function renderList() {
    listContainer.innerHTML = '';
    if (!orderedIds.length) {
      const empty = document.createElement('div');
      empty.className = 'helper-text';
      empty.textContent = 'No control messages. Click "Add Message" to create one.';
      listContainer.appendChild(empty);
      return;
    }
    orderedIds.forEach((id) => {
      const draft = messageMap.get(id);
      if (!draft) return;
      const item = document.createElement('div');
      item.className = 'list-item';
      item.dataset.id = id;
      if (selectedId === id) item.dataset.active = 'true';
      const title = document.createElement('div');
      title.className = 'list-item-title';
      title.textContent = draft.label || draft.id;
      const range = document.createElement('div');
      range.className = 'list-item-sub';
      const display = (ref) => {
        if (!ref) return '';
        const parts = [];
        if (ref.frameId) parts.push(ref.frameId);
        if (ref.offset !== '' && ref.offset != null) parts.push(`offset ${ref.offset}`);
        if (ref.frameIndex !== '' && ref.frameIndex != null) parts.push(`index ${ref.frameIndex}`);
        return parts.join(' · ');
      };
      range.textContent = `${display(draft.range.start)} → ${display(draft.range.end)}`
        || 'Full timeline';
      item.appendChild(title);
      item.appendChild(range);
      item.addEventListener('click', () => {
        selectedId = id;
        renderList();
        renderDetail();
      });
      listContainer.appendChild(item);
    });
  }

  const detailContainer = document.createElement('div');
  detailContainer.className = 'col';
  detailSection.appendChild(detailContainer);

  function getSelectedDraft() {
    if (!selectedId) return null;
    return messageMap.get(selectedId) || null;
  }

  function renderRangeInputs(draft, position) {
    const rangeLabel = document.createElement('label');
    rangeLabel.textContent = `${position === 'start' ? 'Start' : 'End'} reference`;
    const grid = document.createElement('div');
    grid.className = 'range-grid';

    const frameIdInput = document.createElement('input');
    frameIdInput.placeholder = 'frame-00012';
    frameIdInput.value = draft.range[position].frameId || '';
    frameIdInput.addEventListener('input', () => {
      draft.range[position].frameId = frameIdInput.value.trim();
    });

    const frameIndexInput = document.createElement('input');
    frameIndexInput.type = 'number';
    frameIndexInput.placeholder = 'index';
    frameIndexInput.value = draft.range[position].frameIndex || '';
    frameIndexInput.addEventListener('input', () => {
      draft.range[position].frameIndex = frameIndexInput.value;
    });

    const offsetInput = document.createElement('input');
    offsetInput.type = 'number';
    offsetInput.placeholder = 'offset';
    offsetInput.value = draft.range[position].offset || '';
    offsetInput.addEventListener('input', () => {
      draft.range[position].offset = offsetInput.value;
    });

    grid.appendChild(frameIdInput);
    grid.appendChild(frameIndexInput);
    grid.appendChild(offsetInput);

    const actionsRow = document.createElement('div');
    actionsRow.className = 'row';
    const useCurrent = document.createElement('button');
    useCurrent.type = 'button';
    useCurrent.className = 'secondary';
    useCurrent.textContent = 'Use current frame';
    useCurrent.addEventListener('click', () => {
      const offset = options.getCurrentOffset ? options.getCurrentOffset() : -1;
      if (Number.isFinite(offset)) {
        offsetInput.value = String(offset);
        draft.range[position].offset = String(offset);
      }
      if (options.offsetToFrameId) {
        const frameId = options.offsetToFrameId(offset);
        if (frameId) {
          frameIdInput.value = frameId;
          draft.range[position].frameId = frameId;
        }
      }
      if (options.offsetToIndex) {
        const frameIndex = options.offsetToIndex(offset);
        if (frameIndex != null) {
          frameIndexInput.value = String(frameIndex);
          draft.range[position].frameIndex = String(frameIndex);
        }
      }
    });
    actionsRow.appendChild(useCurrent);

    if (position === 'end') {
      const inclusiveToggle = document.createElement('label');
      inclusiveToggle.className = 'action-toggle';
      const inclusiveInput = document.createElement('input');
      inclusiveInput.type = 'checkbox';
      inclusiveInput.checked = draft.range.end.inclusive !== false;
      inclusiveInput.addEventListener('change', () => {
        draft.range.end.inclusive = inclusiveInput.checked;
      });
      inclusiveToggle.appendChild(inclusiveInput);
      inclusiveToggle.appendChild(document.createTextNode('Inclusive'));
      actionsRow.appendChild(inclusiveToggle);
    }

    rangeLabel.appendChild(grid);
    rangeLabel.appendChild(actionsRow);
    return rangeLabel;
  }

  function createActionToggle(labelText, checked, onChange) {
    const wrapper = document.createElement('label');
    wrapper.className = 'action-toggle';
    const input = document.createElement('input');
    input.type = 'checkbox';
    input.checked = !!checked;
    input.addEventListener('change', () => onChange(input.checked));
    wrapper.appendChild(input);
    wrapper.appendChild(document.createTextNode(labelText));
    return { wrapper, input };
  }

  function renderDetail() {
    detailContainer.innerHTML = '';
    const draft = getSelectedDraft();
    if (!draft) {
      const placeholder = document.createElement('div');
      placeholder.className = 'helper-text';
      placeholder.textContent = 'Select a message to edit its details.';
      detailContainer.appendChild(placeholder);
      return;
    }

    const metaRow = document.createElement('div');
    metaRow.className = 'row';

    const idLabel = document.createElement('label');
    idLabel.textContent = 'Message ID';
    const idField = document.createElement('input');
    idField.value = draft.id;
    idField.disabled = true;
    idLabel.appendChild(idField);
    metaRow.appendChild(idLabel);

    const labelLabel = document.createElement('label');
    labelLabel.textContent = 'Label';
    const labelField = document.createElement('input');
    labelField.value = draft.label || '';
    labelField.addEventListener('input', () => { draft.label = labelField.value; });
    labelLabel.appendChild(labelField);
    metaRow.appendChild(labelLabel);

    const priorityLabel = document.createElement('label');
    priorityLabel.textContent = 'Priority';
    const priorityField = document.createElement('input');
    priorityField.type = 'number';
    priorityField.step = '1';
    priorityField.value = String(draft.priority || 0);
    priorityField.addEventListener('input', () => {
      const val = parseNumber(priorityField.value);
      if (val != null) draft.priority = val;
    });
    priorityLabel.appendChild(priorityField);
    metaRow.appendChild(priorityLabel);

    detailContainer.appendChild(metaRow);

    const notesLabel = document.createElement('label');
    notesLabel.textContent = 'Notes';
    const notesField = document.createElement('textarea');
    notesField.value = draft.notes || '';
    notesField.addEventListener('input', () => { draft.notes = notesField.value; });
    notesLabel.appendChild(notesField);
    detailContainer.appendChild(notesLabel);

    detailContainer.appendChild(renderRangeInputs(draft, 'start'));
    detailContainer.appendChild(renderRangeInputs(draft, 'end'));

    const divider = document.createElement('div');
    divider.className = 'section-divider';
    detailContainer.appendChild(divider);

    // Playback speed action
    const speedHeader = document.createElement('div');
    speedHeader.className = 'action-header';
    const speedTitle = document.createElement('h3');
    speedTitle.textContent = 'Playback Speed';
    speedHeader.appendChild(speedTitle);
    const speedToggle = createActionToggle('Enable', draft.actions.speed.enabled, (checked) => {
      draft.actions.speed.enabled = checked;
      renderDetail();
    });
    speedHeader.appendChild(speedToggle.wrapper);
    detailContainer.appendChild(speedHeader);

    if (draft.actions.speed.enabled) {
      const modeRow = document.createElement('div');
      modeRow.className = 'row';
      ['fps', 'multiplier'].forEach((mode) => {
        const btn = document.createElement('button');
        btn.type = 'button';
        btn.className = draft.actions.speed.mode === mode ? '' : 'secondary';
        btn.textContent = mode === 'fps' ? 'Use FPS' : 'Use Multiplier';
        btn.addEventListener('click', () => {
          draft.actions.speed.mode = mode;
          renderDetail();
        });
        modeRow.appendChild(btn);
      });
      detailContainer.appendChild(modeRow);

      if (draft.actions.speed.mode === 'fps') {
        const fpsLabel = document.createElement('label');
        fpsLabel.textContent = 'Frames per second';
        const fpsInput = document.createElement('input');
        fpsInput.type = 'number';
        fpsInput.min = '1';
        fpsInput.step = '1';
        fpsInput.value = draft.actions.speed.fps || '';
        fpsInput.addEventListener('input', () => { draft.actions.speed.fps = fpsInput.value; });
        fpsLabel.appendChild(fpsInput);
        detailContainer.appendChild(fpsLabel);
      } else {
        const multLabel = document.createElement('label');
        multLabel.textContent = 'Speed multiplier';
        const multInput = document.createElement('input');
        multInput.type = 'number';
        multInput.min = '0.1';
        multInput.step = '0.1';
        multInput.value = draft.actions.speed.multiplier || '';
        multInput.addEventListener('input', () => { draft.actions.speed.multiplier = multInput.value; });
        multLabel.appendChild(multInput);
        detailContainer.appendChild(multLabel);
      }

      const transitionRow = document.createElement('div');
      transitionRow.className = 'row';
      const transLabel = document.createElement('label');
      transLabel.textContent = 'Transition (ms)';
      const transInput = document.createElement('input');
      transInput.type = 'number';
      transInput.min = '0';
      transInput.step = '10';
      transInput.value = draft.actions.speed.transitionMs || '';
      transInput.addEventListener('input', () => { draft.actions.speed.transitionMs = transInput.value; });
      transLabel.appendChild(transInput);
      transitionRow.appendChild(transLabel);

      const easingLabel = document.createElement('label');
      easingLabel.textContent = 'Easing';
      const easingSelect = document.createElement('select');
      ['', 'linear', 'ease-in', 'ease-out', 'ease-in-out'].forEach((optionValue) => {
        const opt = document.createElement('option');
        opt.value = optionValue;
        opt.textContent = optionValue || 'Default';
        if (draft.actions.speed.easing === optionValue) opt.selected = true;
        easingSelect.appendChild(opt);
      });
      easingSelect.addEventListener('change', () => { draft.actions.speed.easing = easingSelect.value; });
      easingLabel.appendChild(easingSelect);
      transitionRow.appendChild(easingLabel);
      detailContainer.appendChild(transitionRow);
    }

    const divider2 = document.createElement('div');
    divider2.className = 'section-divider';
    detailContainer.appendChild(divider2);

    // Callout action
    const calloutHeader = document.createElement('div');
    calloutHeader.className = 'action-header';
    const calloutTitle = document.createElement('h3');
    calloutTitle.textContent = 'Callout';
    calloutHeader.appendChild(calloutTitle);
    const calloutToggle = createActionToggle('Enable', draft.actions.callout.enabled, (checked) => {
      draft.actions.callout.enabled = checked;
      renderDetail();
    });
    calloutHeader.appendChild(calloutToggle.wrapper);
    detailContainer.appendChild(calloutHeader);

    if (draft.actions.callout.enabled) {
      const textLabel = document.createElement('label');
      textLabel.textContent = 'Text';
      const textArea = document.createElement('textarea');
      textArea.value = draft.actions.callout.text || '';
      textArea.addEventListener('input', () => { draft.actions.callout.text = textArea.value; });
      textLabel.appendChild(textArea);
      detailContainer.appendChild(textLabel);

      const modeLabel = document.createElement('label');
      modeLabel.textContent = 'Anchor mode';
      const modeSelect = document.createElement('select');
      ['atom', 'bond', 'world'].forEach((mode) => {
        const opt = document.createElement('option');
        opt.value = mode;
        opt.textContent = mode.charAt(0).toUpperCase() + mode.slice(1);
        if (draft.actions.callout.anchorMode === mode) opt.selected = true;
        modeSelect.appendChild(opt);
      });
      modeSelect.addEventListener('change', () => {
        draft.actions.callout.anchorMode = modeSelect.value;
        renderDetail();
      });
      modeLabel.appendChild(modeSelect);
      detailContainer.appendChild(modeLabel);

      const anchorMode = draft.actions.callout.anchorMode || 'atom';
      if (anchorMode === 'atom' || anchorMode === 'bond') {
        const atomsLabel = document.createElement('label');
        atomsLabel.textContent = anchorMode === 'atom' ? 'Atom Index' : 'Bond atoms';
        const atomsInput = document.createElement('input');
        atomsInput.placeholder = anchorMode === 'atom' ? 'e.g. 5' : 'e.g. 2, 15';
        const currentAtoms = draft.actions.callout.anchor?.atoms;
        atomsInput.value = atomsToString(currentAtoms);
        atomsInput.addEventListener('input', () => {
          draft.actions.callout.anchor = draft.actions.callout.anchor || {};
          draft.actions.callout.anchor.atoms = parseAtoms(atomsInput.value);
        });
        atomsLabel.appendChild(atomsInput);
        const useSelection = document.createElement('button');
        useSelection.type = 'button';
        useSelection.className = 'secondary';
        useSelection.textContent = 'Use selection';
        useSelection.addEventListener('click', () => {
          if (!options.getSelectionSnapshot) return;
          const sel = options.getSelectionSnapshot();
          if (anchorMode === 'atom' && sel?.kind === 'atom') {
            atomsInput.value = String(sel.index);
            draft.actions.callout.anchor = { mode: 'atom', atoms: [sel.index] };
          } else if (anchorMode === 'bond' && sel?.kind === 'bond') {
            const atoms = ensureArray(sel.atoms).slice(0, 2);
            atomsInput.value = atomsToString(atoms);
            draft.actions.callout.anchor = { mode: 'bond', atoms };
          }
        });
        atomsLabel.appendChild(useSelection);
        detailContainer.appendChild(atomsLabel);
      } else if (anchorMode === 'world') {
        const posLabel = document.createElement('label');
        posLabel.textContent = 'World position (x, y, z)';
        const posRow = document.createElement('div');
        posRow.className = 'row';
        const worldPosition = draft.actions.callout.worldPosition || draft.actions.callout.anchor?.position || [0, 0, 0];
        for (let i = 0; i < 3; i += 1) {
          const input = document.createElement('input');
          input.type = 'number';
          input.step = '0.01';
          input.value = String(worldPosition[i] ?? '');
          input.addEventListener('input', () => {
            draft.actions.callout.worldPosition = draft.actions.callout.worldPosition || [0, 0, 0];
            draft.actions.callout.worldPosition[i] = Number(input.value) || 0;
          });
          posRow.appendChild(input);
        }
        const useCamera = document.createElement('button');
        useCamera.type = 'button';
        useCamera.className = 'secondary';
        useCamera.textContent = 'Use camera target';
        if (typeof options.getCameraTarget === 'function') {
          useCamera.addEventListener('click', () => {
            const target = options.getCameraTarget();
            if (Array.isArray(target) && target.length === 3) {
              posRow.querySelectorAll('input').forEach((input, index) => {
                input.value = String(target[index] || 0);
              });
              draft.actions.callout.worldPosition = target.slice(0, 3);
            }
          });
        } else {
          useCamera.disabled = true;
        }
        posRow.appendChild(useCamera);
        posLabel.appendChild(posRow);
        detailContainer.appendChild(posLabel);
      }

      const offsetLabel = document.createElement('label');
      offsetLabel.textContent = 'Offset (optional)';
      const offsetRow = document.createElement('div');
      offsetRow.className = 'row';
      for (let i = 0; i < 3; i += 1) {
        const input = document.createElement('input');
        input.type = 'number';
        input.step = '0.01';
        input.placeholder = i === 0 ? 'x' : i === 1 ? 'y' : 'z';
        input.value = draft.actions.callout.offset[i] || '';
        input.addEventListener('input', () => {
          draft.actions.callout.offset[i] = input.value;
        });
        offsetRow.appendChild(input);
      }
      offsetLabel.appendChild(offsetRow);
      detailContainer.appendChild(offsetLabel);

      const panelRow = document.createElement('div');
      panelRow.className = 'row';
      const widthLabel = document.createElement('label');
      widthLabel.textContent = 'Panel width';
      const widthInput = document.createElement('input');
      widthInput.type = 'number';
      widthInput.step = '0.1';
      widthInput.value = draft.actions.callout.panelWidth || '';
      widthInput.addEventListener('input', () => { draft.actions.callout.panelWidth = widthInput.value; });
      widthLabel.appendChild(widthInput);
      panelRow.appendChild(widthLabel);

      const heightLabel = document.createElement('label');
      heightLabel.textContent = 'Panel height';
      const heightInput = document.createElement('input');
      heightInput.type = 'number';
      heightInput.step = '0.1';
      heightInput.value = draft.actions.callout.panelHeight || '';
      heightInput.addEventListener('input', () => { draft.actions.callout.panelHeight = heightInput.value; });
      heightLabel.appendChild(heightInput);
      panelRow.appendChild(heightLabel);
      detailContainer.appendChild(panelRow);

      const styleRow = document.createElement('div');
      styleRow.className = 'row';
      const textSizeLabel = document.createElement('label');
      textSizeLabel.textContent = 'Text size';
      const textSizeInput = document.createElement('input');
      textSizeInput.type = 'number';
      textSizeInput.step = '0.1';
      textSizeInput.value = draft.actions.callout.textSize || '';
      textSizeInput.addEventListener('input', () => { draft.actions.callout.textSize = textSizeInput.value; });
      textSizeLabel.appendChild(textSizeInput);
      styleRow.appendChild(textSizeLabel);

      const bgLabel = document.createElement('label');
      bgLabel.textContent = 'Background';
      const bgInput = document.createElement('input');
      bgInput.type = 'text';
      bgInput.placeholder = 'rgba(...)';
      bgInput.value = draft.actions.callout.background || '';
      bgInput.addEventListener('input', () => { draft.actions.callout.background = bgInput.value; });
      bgLabel.appendChild(bgInput);
      styleRow.appendChild(bgLabel);

      const colorLabel = document.createElement('label');
      colorLabel.textContent = 'Text color';
      const colorInput = document.createElement('input');
      colorInput.type = 'text';
      colorInput.placeholder = '#ffffff';
      colorInput.value = draft.actions.callout.color || '';
      colorInput.addEventListener('input', () => { draft.actions.callout.color = colorInput.value; });
      colorLabel.appendChild(colorInput);
      styleRow.appendChild(colorLabel);

      detailContainer.appendChild(styleRow);
    }

    const divider3 = document.createElement('div');
    divider3.className = 'section-divider';
    detailContainer.appendChild(divider3);

    // Opacity action
    const opacityHeader = document.createElement('div');
    opacityHeader.className = 'action-header';
    const opacityTitle = document.createElement('h3');
    opacityTitle.textContent = 'Opacity Focus';
    opacityHeader.appendChild(opacityTitle);
    const opacityToggle = createActionToggle('Enable', draft.actions.opacity.enabled, (checked) => {
      draft.actions.opacity.enabled = checked;
      renderDetail();
    });
    opacityHeader.appendChild(opacityToggle.wrapper);
    detailContainer.appendChild(opacityHeader);

    if (draft.actions.opacity.enabled) {
      const atomsLabel = document.createElement('label');
      atomsLabel.textContent = 'Focus atoms';
      const atomsInput = document.createElement('input');
      atomsInput.placeholder = 'e.g. 1, 2, 3';
      atomsInput.value = draft.actions.opacity.focusAtomsText || '';
      atomsInput.addEventListener('input', () => { draft.actions.opacity.focusAtomsText = atomsInput.value; });
      atomsLabel.appendChild(atomsInput);
      const useSelection = document.createElement('button');
      useSelection.type = 'button';
      useSelection.className = 'secondary';
      useSelection.textContent = 'Use selection';
      useSelection.addEventListener('click', () => {
        if (!options.getSelectionSnapshot) return;
        const sel = options.getSelectionSnapshot();
        if (sel?.kind === 'atom') {
          atomsInput.value = String(sel.index);
          draft.actions.opacity.focusAtomsText = String(sel.index);
        } else if (sel?.kind === 'bond') {
          const atoms = ensureArray(sel.atoms).slice(0, 2);
          atomsInput.value = atomsToString(atoms);
          draft.actions.opacity.focusAtomsText = atomsToString(atoms);
        }
      });
      atomsLabel.appendChild(useSelection);
      detailContainer.appendChild(atomsLabel);

      const includeLabel = document.createElement('label');
      includeLabel.textContent = 'Include bonds';
      const includeSelect = document.createElement('select');
      ['none', 'connected', 'exact'].forEach((mode) => {
        const opt = document.createElement('option');
        opt.value = mode;
        opt.textContent = mode;
        if (draft.actions.opacity.includeBonds === mode) opt.selected = true;
        includeSelect.appendChild(opt);
      });
      includeSelect.addEventListener('change', () => { draft.actions.opacity.includeBonds = includeSelect.value; });
      includeLabel.appendChild(includeSelect);
      detailContainer.appendChild(includeLabel);

      const opacityRow = document.createElement('div');
      opacityRow.className = 'row';
      const focusAtomsLabel = document.createElement('label');
      focusAtomsLabel.textContent = 'Focus atom opacity';
      const focusAtomsInput = document.createElement('input');
      focusAtomsInput.type = 'number';
      focusAtomsInput.step = '0.05';
      focusAtomsInput.min = '0';
      focusAtomsInput.max = '1';
      focusAtomsInput.value = draft.actions.opacity.focusOpacityAtoms || '';
      focusAtomsInput.addEventListener('input', () => { draft.actions.opacity.focusOpacityAtoms = focusAtomsInput.value; });
      focusAtomsLabel.appendChild(focusAtomsInput);
      opacityRow.appendChild(focusAtomsLabel);

      const focusBondsLabel = document.createElement('label');
      focusBondsLabel.textContent = 'Focus bond opacity';
      const focusBondsInput = document.createElement('input');
      focusBondsInput.type = 'number';
      focusBondsInput.step = '0.05';
      focusBondsInput.min = '0';
      focusBondsInput.max = '1';
      focusBondsInput.value = draft.actions.opacity.focusOpacityBonds || '';
      focusBondsInput.addEventListener('input', () => { draft.actions.opacity.focusOpacityBonds = focusBondsInput.value; });
      focusBondsLabel.appendChild(focusBondsInput);
      opacityRow.appendChild(focusBondsLabel);
      detailContainer.appendChild(opacityRow);

      const backgroundRow = document.createElement('div');
      backgroundRow.className = 'row';
      const bgAtomsLabel = document.createElement('label');
      bgAtomsLabel.textContent = 'Background atom opacity';
      const bgAtomsInput = document.createElement('input');
      bgAtomsInput.type = 'number';
      bgAtomsInput.step = '0.05';
      bgAtomsInput.min = '0';
      bgAtomsInput.max = '1';
      bgAtomsInput.value = draft.actions.opacity.backgroundAtoms || '';
      bgAtomsInput.addEventListener('input', () => { draft.actions.opacity.backgroundAtoms = bgAtomsInput.value; });
      bgAtomsLabel.appendChild(bgAtomsInput);
      backgroundRow.appendChild(bgAtomsLabel);

      const bgBondsLabel = document.createElement('label');
      bgBondsLabel.textContent = 'Background bond opacity';
      const bgBondsInput = document.createElement('input');
      bgBondsInput.type = 'number';
      bgBondsInput.step = '0.05';
      bgBondsInput.min = '0';
      bgBondsInput.max = '1';
      bgBondsInput.value = draft.actions.opacity.backgroundBonds || '';
      bgBondsInput.addEventListener('input', () => { draft.actions.opacity.backgroundBonds = bgBondsInput.value; });
      bgBondsLabel.appendChild(bgBondsInput);
      backgroundRow.appendChild(bgBondsLabel);
      detailContainer.appendChild(backgroundRow);

      const transitionLabel = document.createElement('label');
      transitionLabel.textContent = 'Transition (ms)';
      const transitionInput = document.createElement('input');
      transitionInput.type = 'number';
      transitionInput.min = '0';
      transitionInput.step = '10';
      transitionInput.value = draft.actions.opacity.transitionMs || '';
      transitionInput.addEventListener('input', () => { draft.actions.opacity.transitionMs = transitionInput.value; });
      transitionLabel.appendChild(transitionInput);
      detailContainer.appendChild(transitionLabel);
    }

    const saveRow = document.createElement('div');
    saveRow.className = 'row';
    const saveButton = document.createElement('button');
    saveButton.type = 'button';
    saveButton.textContent = 'Save Changes';
    saveButton.addEventListener('click', () => {
      const result = commitChanges();
      if (result.ok) {
        setError('');
      } else {
        setError(result.error || 'Validation failed');
      }
    });
    saveRow.appendChild(saveButton);

    const refreshButton = document.createElement('button');
    refreshButton.type = 'button';
    refreshButton.className = 'secondary';
    refreshButton.textContent = 'Reload state';
    refreshButton.addEventListener('click', () => {
      loadMessages();
      playbackDraft = options.getPlaybackConfig ? deepClone(options.getPlaybackConfig()) : playbackDraft;
      renderPlayback();
      renderList();
      renderDetail();
      setError('');
    });
    saveRow.appendChild(refreshButton);

    detailContainer.appendChild(saveRow);
  }

  renderList();
  renderDetail();

  function commitChanges() {
    const allDrafts = orderedIds.map((id) => messageMap.get(id)).filter(Boolean);
    const serialized = [];
    for (const draft of allDrafts) {
      const actions = buildActionsFromDraft(draft.actions);
      if (!actions.length) {
        return { ok: false, error: `Message "${draft.label || draft.id}" has no enabled actions.` };
      }
      const payload = {
        id: draft.id,
        label: draft.label,
        priority: draft.priority,
        notes: draft.notes,
        actions,
      };
      const rangePayload = buildRangePayload(draft.range);
      if (rangePayload) payload.range = rangePayload;
      serialized.push(payload);
    }

    try {
      if (options.setMessages) options.setMessages(serialized);
      if (options.refreshControlEngine) options.refreshControlEngine();
      const playbackPayload = deepClone(playbackDraft);
      if (options.setPlaybackConfig) options.setPlaybackConfig(playbackPayload);
      if (options.onCommit) options.onCommit({
        messages: serialized,
        playback: playbackPayload,
      });
      return { ok: true };
    } catch (err) {
      return { ok: false, error: err?.message || String(err) };
    }
  }

  return {
    destroy() {
      try {
        if (panel.parentNode) panel.parentNode.removeChild(panel);
      } catch { /* noop */ }
    },
    refresh() {
      loadMessages();
      playbackDraft = options.getPlaybackConfig ? deepClone(options.getPlaybackConfig()) : playbackDraft;
      renderPlayback();
      renderList();
      renderDetail();
    },
    getState() {
      return {
        selectedId,
        messageCount: orderedIds.length,
        playback: deepClone(playbackDraft),
      };
    },
    select(id) {
      if (!messageMap.has(id)) return false;
      selectedId = id;
      renderList();
      renderDetail();
      return true;
    },
    getCurrentDraft() {
      return deepClone(getSelectedDraft());
    },
  };
}

export default { createTimelineEditorPanel };
