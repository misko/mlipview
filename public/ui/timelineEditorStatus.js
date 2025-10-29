const STYLE_ID = 'timelineEditorStatusStyles';

function injectStyles() {
  if (typeof document === 'undefined') return;
  if (document.getElementById(STYLE_ID)) return;
  const style = document.createElement('style');
  style.id = STYLE_ID;
  style.textContent = `
    #timelineEditorStatus {
      position: absolute;
      left: 50%;
      bottom: 12px;
      transform: translateX(-50%);
      display: flex;
      align-items: center;
      gap: 24px;
      padding: 8px 18px;
      background: rgba(16, 20, 32, 0.88);
      border: 1px solid rgba(120, 150, 200, 0.18);
      border-radius: 999px;
      color: #e8f0ff;
      font: 12px/1.5 system-ui, sans-serif;
      pointer-events: none;
      z-index: 40;
      backdrop-filter: blur(6px);
    }
    #timelineEditorStatus .status-label {
      font-weight: 600;
      letter-spacing: .04em;
      text-transform: uppercase;
      color: rgba(208, 216, 232, 0.75);
    }
    #timelineEditorStatus .status-value {
      font-weight: 600;
      color: #fdfdfd;
    }
    #timelineEditorStatus .status-block {
      display: flex;
      flex-direction: column;
      gap: 2px;
      text-align: center;
      min-width: 120px;
    }
  `;
  document.head.appendChild(style);
}

export function createTimelineEditorStatus({ attachTo } = {}) {
  if (typeof document === 'undefined') {
    return {
      update() {},
      destroy() {},
      getElement() {
        return null;
      },
    };
  }

  injectStyles();

  const host = attachTo || document.getElementById('app') || document.body;
  const bar = document.createElement('div');
  bar.id = 'timelineEditorStatus';

  const frameBlock = document.createElement('div');
  frameBlock.className = 'status-block';
  const frameLabel = document.createElement('div');
  frameLabel.className = 'status-label';
  frameLabel.textContent = 'Frame';
  const frameValue = document.createElement('div');
  frameValue.className = 'status-value';
  frameValue.textContent = 'Live';
  frameBlock.appendChild(frameLabel);
  frameBlock.appendChild(frameValue);

  const selectionBlock = document.createElement('div');
  selectionBlock.className = 'status-block';
  const selectionLabel = document.createElement('div');
  selectionLabel.className = 'status-label';
  selectionLabel.textContent = 'Selection';
  const selectionValue = document.createElement('div');
  selectionValue.className = 'status-value';
  selectionValue.textContent = 'None';
  selectionBlock.appendChild(selectionLabel);
  selectionBlock.appendChild(selectionValue);

  bar.appendChild(frameBlock);
  bar.appendChild(selectionBlock);
  host.appendChild(bar);

  function formatFrameDisplay(info = {}) {
    if (info.isLive) return 'Live';
    const parts = [];
    if (info.frameId) parts.push(info.frameId);
    if (Number.isFinite(info.offset)) parts.push(`offset ${info.offset}`);
    if (Number.isFinite(info.frameIndex)) parts.push(`index ${info.frameIndex}`);
    return parts.length ? parts.join(' · ') : '—';
  }

  function formatSelectionDisplay(sel = null) {
    if (!sel) return 'None';
    if (sel.kind === 'atom') {
      const coords = Array.isArray(sel.position)
        ? sel.position.map((v) => Number(v).toFixed(2)).join(', ')
        : '';
      return `Atom ${sel.index}${sel.element ? ` (${sel.element})` : ''}${coords ? ` · [${coords}]` : ''}`;
    }
    if (sel.kind === 'bond') {
      const atoms = sel.atoms || [];
      const label = atoms.length === 2 ? `${atoms[0]}–${atoms[1]}` : 'Bond';
      const elements = sel.elements ? ` (${sel.elements.join('–')})` : '';
      return `Bond ${label}${elements}`;
    }
    return 'None';
  }

  return {
    update({ frame, selection } = {}) {
      frameValue.textContent = formatFrameDisplay(frame);
      selectionValue.textContent = formatSelectionDisplay(selection);
    },
    destroy() {
      try {
        if (bar.parentNode) bar.parentNode.removeChild(bar);
      } catch { /* noop */ }
    },
    getElement() {
      return bar;
    },
  };
}

export default { createTimelineEditorStatus };
