// Callout overlay rendered above the molecule timeline playback.

function toVector3(arr) {
  if (!Array.isArray(arr)) return null;
  if (arr.length < 3) return null;
  return { x: Number(arr[0]) || 0, y: Number(arr[1]) || 0, z: Number(arr[2]) || 0 };
}

export function createCalloutLayer({
  scene = null,
  babylon = (typeof BABYLON !== 'undefined' ? BABYLON : null),
  renderingGroupId = 3,
  getAtomPositionWorld = () => null,
  getBondMidpointWorld = () => null,
} = {}) {
  const state = {
    current: null,
    mesh: null,
    adt: null,
    textBlock: null,
  };

  function disposeVisual() {
    try { state.adt?.dispose?.(); } catch { }
    try { state.mesh?.dispose?.(); } catch { }
    state.mesh = null;
    state.adt = null;
    state.textBlock = null;
  }

  function ensureVisual(cfg) {
    if (!babylon || !scene) return;
    const width = Number(cfg.panelSize?.width) > 0 ? Number(cfg.panelSize.width) : 1.4;
    const height = Number(cfg.panelSize?.height) > 0 ? Number(cfg.panelSize.height) : 0.7;
    if (!state.mesh) {
      const plane = babylon.MeshBuilder?.CreatePlane
        ? babylon.MeshBuilder.CreatePlane('timeline-callout', { width: 1, height: 1 }, scene)
        : null;
      if (!plane) return;
      plane.billboardMode = babylon.Mesh.BILLBOARDMODE_ALL;
      plane.isPickable = false;
      plane.renderingGroupId = renderingGroupId;
      state.mesh = plane;
    }
    const baseWidth = 1;
    const baseHeight = 1;
    state.mesh.scaling.x = width / baseWidth;
    state.mesh.scaling.y = height / baseHeight;
    if (!state.adt && babylon.GUI?.AdvancedDynamicTexture) {
      state.adt = babylon.GUI.AdvancedDynamicTexture.CreateForMesh(state.mesh, 1024, 512, false);
      state.adt.wrap = false;
      const rect = new babylon.GUI.Rectangle('timeline-callout-rect');
      rect.thickness = 0;
      rect.background = 'rgba(10,14,24,0.85)';
      rect.alpha = 1;
      rect.width = 1;
      rect.height = 1;
      state.adt.addControl(rect);
      const tb = new babylon.GUI.TextBlock('timeline-callout-text');
      tb.textHorizontalAlignment = babylon.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
      tb.textVerticalAlignment = babylon.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
      tb.color = '#f0f7ff';
      tb.textWrapping = true;
      tb.fontFamily = 'system-ui, sans-serif';
      rect.addControl(tb);
      state.textBlock = tb;
    }
    if (state.textBlock) {
      state.textBlock.text = cfg.text || '';
      const size = Number(cfg.textSize);
      state.textBlock.fontSize = Number.isFinite(size) && size > 0 ? size * 48 : 42;
      if (cfg.style && typeof cfg.style === 'object') {
        if (cfg.style.color) state.textBlock.color = cfg.style.color;
        if (cfg.style.fontFamily) state.textBlock.fontFamily = cfg.style.fontFamily;
      }
    }
    if (state.adt && cfg.style && typeof cfg.style === 'object') {
      const rect = state.adt.rootContainer?.children?.[0];
      if (rect) {
        if (cfg.style.background) rect.background = cfg.style.background;
        if (cfg.style.opacity != null) rect.alpha = Number(cfg.style.opacity);
        if (cfg.style.cornerRadius != null) rect.cornerRadius = Number(cfg.style.cornerRadius);
      }
    }
    if (typeof state.mesh.setEnabled === 'function') {
      state.mesh.setEnabled(true);
    } else {
      state.mesh.isVisible = true;
    }
  }

  function computeAnchorPosition(cfg) {
    if (!cfg || !cfg.anchor) return null;
    const anchor = cfg.anchor;
    const offset = Array.isArray(cfg.offset) ? cfg.offset : Array.isArray(anchor.offset) ? anchor.offset : null;
    let base = null;
    switch (anchor.mode) {
      case 'world':
      case 'xyz': {
        base = toVector3(anchor.position || anchor.xyz || anchor.value);
        break;
      }
      case 'atom': {
        const idx = Number.isInteger(anchor.atom) ? anchor.atom : Array.isArray(anchor.atoms) ? anchor.atoms[0] : null;
        if (Number.isInteger(idx)) base = getAtomPositionWorld(idx);
        break;
      }
      case 'bond': {
        let atoms = Array.isArray(anchor.atoms) ? anchor.atoms : null;
        if (!atoms && Number.isInteger(anchor.i) && Number.isInteger(anchor.j)) {
          atoms = [anchor.i, anchor.j];
        }
        if (atoms && atoms.length >= 2) {
          base = getBondMidpointWorld(atoms[0], atoms[1]);
        }
        break;
      }
      default:
        break;
    }
    if (!base) return null;
    if (offset && offset.length === 3) {
      return {
        x: base.x + (Number(offset[0]) || 0),
        y: base.y + (Number(offset[1]) || 0),
        z: base.z + (Number(offset[2]) || 0),
      };
    }
    return base;
  }

  function updatePosition(cfg) {
    if (!state.mesh || !cfg) return;
    const pos = computeAnchorPosition(cfg);
    if (!pos) return;
    state.mesh.position.x = pos.x;
    state.mesh.position.y = pos.y;
    state.mesh.position.z = pos.z;
  }

  function clear() {
    state.current = null;
    if (state.mesh) {
      if (typeof state.mesh.setEnabled === 'function') state.mesh.setEnabled(false);
      else state.mesh.isVisible = false;
    }
  }

  return {
    show(config) {
      state.current = config ? { ...config } : null;
      if (!state.current) {
        clear();
        return;
      }
      ensureVisual(state.current);
      updatePosition(state.current);
    },
    update() {
      if (!state.current) return;
      updatePosition(state.current);
    },
    hide: clear,
    dispose: disposeVisual,
    getState: () => ({
      visible: !!state.current,
      config: state.current ? { ...state.current } : null,
    }),
  };
}

export default { createCalloutLayer };
