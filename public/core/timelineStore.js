// Timeline frame ring buffer
// Stores the last N frames with flattened payload snapshots for playback.

const DEFAULT_CAPACITY = 500;

function nowMs() {
  try {
    return typeof performance !== 'undefined' && performance.now
      ? performance.now()
      : Date.now();
  } catch {
    return Date.now();
  }
}

function toFloat32Triples(triples) {
  if (!triples) return null;
  if (ArrayBuffer.isView(triples) && triples.length) {
    const buf = new Float32Array(triples.length);
    for (let i = 0; i < triples.length; i++) {
      buf[i] = Number(triples[i]) || 0;
    }
    return buf;
  }
  if (!Array.isArray(triples) || !triples.length) return null;
  const out = new Float32Array(triples.length * 3);
  let k = 0;
  for (let i = 0; i < triples.length; i++) {
    const t = triples[i] || [0, 0, 0];
    out[k++] = Number(t[0]) || 0;
    out[k++] = Number(t[1]) || 0;
    out[k++] = Number(t[2]) || 0;
  }
  return out;
}

function vecToArray(vec) {
  if (!vec) return [0, 0, 0];
  if (Array.isArray(vec)) {
    return [Number(vec[0]) || 0, Number(vec[1]) || 0, Number(vec[2]) || 0];
  }
  return [Number(vec.x) || 0, Number(vec.y) || 0, Number(vec.z) || 0];
}

function matrixRow(row) {
  return [
    Number(row?.[0]) || 0,
    Number(row?.[1]) || 0,
    Number(row?.[2]) || 0,
  ];
}

export function expandFloat32Triples(buffer) {
  if (!buffer || !buffer.length) return [];
  const triples = new Array(Math.floor(buffer.length / 3));
  for (let i = 0, k = 0; i < triples.length; i++) {
    triples[i] = [buffer[k++], buffer[k++], buffer[k++]];
  }
  return triples;
}

export function createTimelineStore({ capacity = DEFAULT_CAPACITY } = {}) {
  const cap = Math.max(1, Number(capacity) || DEFAULT_CAPACITY);
  const frames = new Array(cap);
  const subs = new Set();
  let cursor = 0; // next insertion slot
  let size = 0;
  let nextId = 1;
  let version = 0;

  function notify(type = 'update', payload) {
    const snapshot = getSnapshot();
    snapshot.type = type;
    snapshot.meta = payload || null;
    for (const fn of subs) {
      try {
        fn(snapshot);
      } catch {
        /* ignore */
      }
    }
  }

  function add(frame) {
    if (!frame || !frame.payload) return null;
    const id = frame.id != null ? frame.id : nextId++;
    if (id >= nextId) nextId = id + 1;
    const prepared = {
      id,
      timestamp: Number(frame.timestamp) || nowMs(),
      receivedAt: Number(frame.receivedAt) || Date.now(),
      kind: frame.kind || 'idle',
      simStep: Number.isFinite(frame.simStep) ? frame.simStep : undefined,
      userInteractionCount: Number.isFinite(frame.userInteractionCount)
        ? frame.userInteractionCount
        : undefined,
      stateVersion: {
        userInteraction: Number(frame.stateVersion?.userInteraction) || 0,
        totalInteraction: Number(frame.stateVersion?.totalInteraction) || 0,
      },
      payload: {
        positions: toFloat32Triples(frame.payload.positions),
        velocities: toFloat32Triples(frame.payload.velocities),
        forces: toFloat32Triples(frame.payload.forces),
        energy:
          typeof frame.payload.energy === 'number'
            ? frame.payload.energy
            : undefined,
        temperature:
          typeof frame.payload.temperature === 'number'
            ? frame.payload.temperature
            : undefined,
        bonds: Array.isArray(frame.payload.bonds)
          ? frame.payload.bonds.map((b) => ({
              i: Number.isInteger(b.i) ? b.i | 0 : Number(b.i) || 0,
              j: Number.isInteger(b.j) ? b.j | 0 : Number(b.j) || 0,
              length: typeof b.length === 'number' ? b.length : 0,
              opacity: typeof b.opacity === 'number' ? b.opacity : 1,
              crossing: !!b.crossing,
            }))
          : undefined,
        cell: (() => {
          const cell = frame.payload.cell;
          if (!cell) return undefined;
          const rows = Array.isArray(cell)
            ? cell
            : Array.isArray(cell.matrix)
              ? cell.matrix
              : null;
          if (!rows || rows.length !== 3) return undefined;
          return {
            a: matrixRow(rows[0]),
            b: matrixRow(rows[1]),
            c: matrixRow(rows[2]),
            originOffset: vecToArray(
              cell.originOffset || cell.origin || [0, 0, 0],
            ),
            enabled: !!cell.enabled,
          };
        })(),
        stress: Array.isArray(frame.payload.stress)
          ? frame.payload.stress.slice(0, 9).map((v) => Number(v) || 0)
          : undefined,
      },
    };
    const evicted = frames[cursor];
    frames[cursor] = prepared;
    cursor = (cursor + 1) % cap;
    if (size < cap) size++;
    version++;
    notify('add', { added: prepared, evicted });
    return prepared;
  }

  function getSnapshot() {
    const out = [];
    for (let i = 0; i < size; i++) {
      const idx = (cursor - size + i + cap) % cap;
      const frame = frames[idx];
      if (frame) out.push(frame);
    }
    return {
      frames: out,
      capacity: cap,
      size,
      liveIndex: size ? size - 1 : -1,
      version,
    };
  }

  function get(index) {
    if (!size) return null;
    if (index < 0 || index >= size) return null;
    const idx = (cursor - size + index + cap) % cap;
    return frames[idx] || null;
  }

  function clear() {
    for (let i = 0; i < cap; i++) frames[i] = undefined;
    cursor = 0;
    size = 0;
    version++;
    notify('clear');
  }

  function liveFrame() {
    if (!size) return null;
    return get(size - 1);
  }

  function subscribe(fn) {
    if (typeof fn !== 'function') return () => {};
    subs.add(fn);
    try {
      fn(getSnapshot());
    } catch {
      /* ignore */
    }
    return () => subs.delete(fn);
  }

  function getInfo() {
    return {
      size,
      capacity: cap,
      liveIndex: size ? size - 1 : -1,
      version,
    };
  }

  return {
    add,
    get,
    liveFrame,
    clear,
    getSnapshot,
    getInfo,
    subscribe,
    capacity: cap,
  };
}

export default { createTimelineStore, expandFloat32Triples };
