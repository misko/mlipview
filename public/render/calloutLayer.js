// Callout overlay rendered above the molecule timeline playback.

function toVector3(arr) {
  if (!Array.isArray(arr) || arr.length < 3) return null;
  return {
    x: Number(arr[0]) || 0,
    y: Number(arr[1]) || 0,
    z: Number(arr[2]) || 0,
  };
}

export function createCalloutLayer({
  scene = null,
  babylon = (typeof BABYLON !== 'undefined' ? BABYLON : null),
  renderingGroupId = 3,
  getAtomPositionWorld = () => null,
  getBondMidpointWorld = () => null,
} = {}) {
  const entries = new Map();

  function disposeEntry(entry) {
    if (!entry) return;
    try { entry.adt?.dispose?.(); } catch {}
    try { entry.mesh?.dispose?.(); } catch {}
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
        if (!atoms && Number.isInteger(anchor.i) && Number.isInteger(anchor.j)) atoms = [anchor.i, anchor.j];
        if (atoms && atoms.length >= 2) base = getBondMidpointWorld(atoms[0], atoms[1]);
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

  function ensureVisual(entry, cfg) {
    if (!babylon || !scene) return;
    const width = Number(cfg.panelSize?.width) > 0 ? Number(cfg.panelSize.width) : 1.4;
    const height = Number(cfg.panelSize?.height) > 0 ? Number(cfg.panelSize.height) : 0.7;
    if (!entry.mesh) {
      const plane = babylon.MeshBuilder?.CreatePlane
        ? babylon.MeshBuilder.CreatePlane(`timeline-callout-${entry.key}`, { width: 1, height: 1 }, scene)
        : null;
      if (!plane) return;
      plane.billboardMode = babylon.Mesh.BILLBOARDMODE_ALL;
      plane.isPickable = false;
      plane.renderingGroupId = renderingGroupId;
      entry.mesh = plane;
    }
    entry.mesh.scaling.x = width;
    entry.mesh.scaling.y = height;
    if (!entry.adt && babylon.GUI?.AdvancedDynamicTexture) {
      entry.adt = babylon.GUI.AdvancedDynamicTexture.CreateForMesh(entry.mesh, 1024, 512, false);
      entry.adt.wrap = false;
      const rect = new babylon.GUI.Rectangle(`timeline-callout-rect-${entry.key}`);
      rect.thickness = 0;
      rect.background = 'rgba(10,14,24,0.85)';
      rect.alpha = 1;
      rect.width = 1;
      rect.height = 1;
      entry.adt.addControl(rect);
      const tb = new babylon.GUI.TextBlock(`timeline-callout-text-${entry.key}`);
      tb.textHorizontalAlignment = babylon.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
      tb.textVerticalAlignment = babylon.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
      tb.color = '#f0f7ff';
      tb.textWrapping = true;
      tb.fontFamily = 'system-ui, sans-serif';
      rect.addControl(tb);
      entry.textBlock = tb;
    }
    if (entry.textBlock) {
      entry.textBlock.text = cfg.text || '';
      const size = Number(cfg.textSize);
      entry.textBlock.fontSize = Number.isFinite(size) && size > 0 ? size * 48 : 42;
      if (cfg.style && typeof cfg.style === 'object') {
        if (cfg.style.color) entry.textBlock.color = cfg.style.color;
        if (cfg.style.fontFamily) entry.textBlock.fontFamily = cfg.style.fontFamily;
      }
    }
    if (entry.adt && cfg.style && typeof cfg.style === 'object') {
      const rect = entry.adt.rootContainer?.children?.[0];
      if (rect) {
        if (cfg.style.background) rect.background = cfg.style.background;
        if (cfg.style.opacity != null) rect.alpha = Number(cfg.style.opacity);
        if (cfg.style.cornerRadius != null) rect.cornerRadius = Number(cfg.style.cornerRadius);
      }
    }
    if (typeof entry.mesh.setEnabled === 'function') entry.mesh.setEnabled(true);
    else entry.mesh.isVisible = true;
  }

  function updatePosition(entry) {
    if (!entry?.mesh || !entry.config) return;
    const pos = computeAnchorPosition(entry.config);
    if (!pos) return;
    entry.mesh.position.x = pos.x;
    entry.mesh.position.y = pos.y;
    entry.mesh.position.z = pos.z;
  }

  function showAll(configs = []) {
    const keep = new Set();
    configs.forEach((cfg, index) => {
      const key = cfg.key || cfg.sourceId || `callout-${index}`;
      let entry = entries.get(key);
      if (!entry) {
        entry = { key, mesh: null, adt: null, textBlock: null, config: null };
        entries.set(key, entry);
      }
      entry.config = { ...cfg };
      ensureVisual(entry, entry.config);
      updatePosition(entry);
      keep.add(key);
    });
    for (const [key, entry] of entries.entries()) {
      if (!keep.has(key)) {
        disposeEntry(entry);
        entries.delete(key);
      }
    }
  }

  function hideAll() {
    for (const [key, entry] of entries.entries()) {
      disposeEntry(entry);
      entries.delete(key);
    }
  }

  return {
    show(config) {
      showAll(config ? [config] : []);
    },
    showAll,
    update() {
      for (const entry of entries.values()) {
        updatePosition(entry);
      }
    },
    hide: hideAll,
    dispose() {
      for (const entry of entries.values()) disposeEntry(entry);
      entries.clear();
    },
    getState: () => ({
      visible: entries.size > 0,
      count: entries.size,
      configs: Array.from(entries.values()).map((entry) => ({ ...entry.config })),
    }),
  };
}

export default { createCalloutLayer };
