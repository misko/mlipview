// Ring buffer for storing recent simulation frames with cloned payloads.
// Stores up to `capacity` frames; offsets are negative (e.g. -1 latest, -capacity oldest).

const DEFAULT_CAPACITY = 500;

const KIND_IDLE = 'idle';
const KIND_MD = 'md';
const KIND_RELAX = 'relax';

function isFiniteNumber(v) {
  return typeof v === 'number' && Number.isFinite(v);
}

function toFloat32Array(triples) {
  if (!Array.isArray(triples)) return null;
  const n = triples.length;
  if (!n) return new Float32Array(0);
  const out = new Float32Array(n * 3);
  let k = 0;
  for (let i = 0; i < n; i++) {
    const t = triples[i] || [0, 0, 0];
    out[k++] = +t[0] || 0;
    out[k++] = +t[1] || 0;
    out[k++] = +t[2] || 0;
  }
  return out;
}

function fromFloat32Array(arr) {
  if (!(arr instanceof Float32Array) && !(arr instanceof Float64Array)) {
    if (Array.isArray(arr)) {
      // Already triple array
      return arr.map((t) => {
        if (!t || t.length !== 3) return [0, 0, 0];
        return [Number(t[0]) || 0, Number(t[1]) || 0, Number(t[2]) || 0];
      });
    }
    return [];
  }
  const n = Math.floor(arr.length / 3);
  const out = new Array(n);
  let k = 0;
  for (let i = 0; i < n; i++) {
    out[i] = [arr[k++] || 0, arr[k++] || 0, arr[k++] || 0];
  }
  return out;
}

function cloneScalars(arr) {
  if (!Array.isArray(arr)) return [];
  return arr.map((v) => [Number(v[0]) || 0, Number(v[1]) || 0, Number(v[2]) || 0]);
}

function toSignature(f32) {
  if (!(f32 instanceof Float32Array)) {
    const copy = toFloat32Array(f32);
    if (!copy) return 0;
    f32 = copy;
  }
  // Simple FNV-1a hash over float bytes (converted to int via view)
  const view = new DataView(f32.buffer, f32.byteOffset, f32.byteLength);
  let hash = 0x811c9dc5;
  for (let i = 0; i < view.byteLength; i++) {
    hash ^= view.getUint8(i);
    hash = (hash * 0x01000193) >>> 0;
  }
  return hash >>> 0;
}

function normalizeKind(kind) {
  if (kind === KIND_MD || kind === KIND_RELAX || kind === KIND_IDLE) return kind;
  return KIND_IDLE;
}

export function createFrameBuffer({ capacity = DEFAULT_CAPACITY } = {}) {
  const cap = Math.max(1, Math.floor(capacity) || DEFAULT_CAPACITY);
  const store = new Array(cap);
  let head = -1; // points at latest index
  let count = 0;
  let nextId = 1;

  function record(kind, frame, meta = {}) {
    if (!frame || !Array.isArray(frame.positions) || frame.positions.length === 0) {
      return null;
    }

    const pos = toFloat32Array(frame.positions);
    if (!pos) return null;

    const velocities = Array.isArray(frame.velocities) ? toFloat32Array(frame.velocities) : null;
    const forces = Array.isArray(frame.forces) ? toFloat32Array(frame.forces) : null;

    head = (head + 1) % cap;
    count = Math.min(cap, count + 1);

    const id = nextId++;
    const simStep = isFiniteNumber(frame.simStep) ? frame.simStep : isFiniteNumber(frame.sim_step) ? frame.sim_step : null;
    const userInteractionCount = isFiniteNumber(frame.userInteractionCount)
      ? frame.userInteractionCount
      : isFiniteNumber(frame.userInteractionCount ?? frame.user_interaction_count)
        ? (frame.userInteractionCount ?? frame.user_interaction_count)
        : null;
    const seq = isFiniteNumber(frame.seq) ? frame.seq : 0;
    const energyIndex = Number.isFinite(meta.energyIndex) ? (meta.energyIndex | 0) : null;

    store[head] = {
      id,
      kind: normalizeKind(kind),
      timestamp: Date.now(),
      seq,
      simStep,
      userInteractionCount,
      positions: pos,
      velocities,
      forces,
      energy: isFiniteNumber(frame.energy) ? frame.energy : null,
      temperature: isFiniteNumber(frame.temperature) ? frame.temperature : null,
      stress: Array.isArray(frame.stress) ? cloneScalars(frame.stress) : null,
      signature: toSignature(pos),
      energyIndex,
    };

    return { id, offset: -1 };
  }

  function size() {
    return count;
  }

  function indexForOffset(offset) {
    if (offset === 0) return -1;
    if (!count) return -1;
    const off = Math.abs(Math.floor(offset));
    if (off < 1 || off > count) return -1;
    let idx = head - (off - 1);
    if (idx < 0) idx += cap;
    return idx;
  }

  function getByOffset(offset) {
    const idx = indexForOffset(offset);
    if (idx < 0) return null;
    const entry = store[idx];
    if (!entry) return null;
    return {
      id: entry.id,
      kind: entry.kind,
      seq: entry.seq,
      simStep: entry.simStep,
      userInteractionCount: entry.userInteractionCount,
      timestamp: entry.timestamp,
      energy: entry.energy,
      temperature: entry.temperature,
      stress: entry.stress ? entry.stress.map((row) => row.slice()) : null,
      signature: entry.signature,
      energyIndex: entry.energyIndex != null ? entry.energyIndex : null,
      positions: fromFloat32Array(entry.positions),
      velocities: entry.velocities ? fromFloat32Array(entry.velocities) : null,
      forces: entry.forces ? fromFloat32Array(entry.forces) : null,
    };
  }

  function listOffsets() {
    const out = [];
    for (let i = 1; i <= count; i++) out.push(-i);
    return out;
  }

  function latestOffset() {
    return count ? -1 : null;
  }

  function getSignature(offset) {
    const idx = indexForOffset(offset);
    if (idx < 0) return null;
    return store[idx]?.signature ?? null;
  }

  function stats() {
    return {
      capacity: cap,
      size: count,
      latestOffset: latestOffset(),
    };
  }

  return {
    record,
    size,
    stats,
    listOffsets,
    latestOffset,
    getByOffset,
    getSignature,
  };
}

export default { createFrameBuffer };
