// Thin WebSocket client for fairchem_local_server2 /ws protobuf stream
// Import modern ESM protobuf stubs directly
import { __count } from './util/funcCount.js';
import { create, toBinary, fromBinary } from '@bufbuild/protobuf';
import {
  // Core messages (oneof + flat arrays)
  ClientActionSchema,
  ClientAction_UserInteractionSchema,
  ClientAction_StartSchema,
  ClientAction_Start_SimType,
  ClientAction_StopSchema,
  ClientAction_PingSchema,
  ServerResultSchema,
  // Supporting messages
  Mat3Schema,
  SimulationParamsSchema,
} from '/proto/fairchem_local_server2/session_pb.js';

const __hasPB = !!ClientActionSchema && !!ServerResultSchema;
async function __ensurePB() {
  if (!__hasPB) throw new Error('[WS][protobuf-missing] ESM stubs not found');
  return true;
}
function __n(v) {
  return typeof v === 'bigint' ? Number(v) : typeof v === 'number' ? v : 0;
}
// Flatten [[x,y,z], ...] -> [x0,y0,z0, x1,y1,z1, ...]
function __triplesToFlat3(triples) {
  try {
    if (!Array.isArray(triples)) return undefined;
    const out = new Array(triples.length * 3);
    let j = 0;
    for (let i = 0; i < triples.length; i++) {
      const p = triples[i] || [0, 0, 0];
      out[j++] = +p[0];
      out[j++] = +p[1];
      out[j++] = +p[2];
    }
    return out;
  } catch {
    return undefined;
  }
}
// Expand flat [x0,y0,z0, ...] -> [[x,y,z], ...]
function __flat3ToTriples(arr) {
  try {
    if (!Array.isArray(arr)) return undefined;
    const n = Math.floor(arr.length / 3);
    const out = new Array(n);
    for (let i = 0, k = 0; i < n; i++) out[i] = [arr[k++], arr[k++], arr[k++]];
    return out;
  } catch {
    return undefined;
  }
}
function __wsDebugOn() {
  try {
    if (typeof window !== 'undefined') {
      // URL param wsDebug=1|true toggles
      try {
        const q = new URLSearchParams(window.location.search || '');
        if (
          q.get('wsDebug') === '1' ||
          q.get('wsDebug') === 'true' ||
          q.get('debug') === '1' ||
          q.get('debug') === 'true'
        )
          window.__MLIPVIEW_DEBUG_WS = true;
      } catch { }
      // localStorage key 'mlip.wsDebug' toggles
      try {
        if (window.localStorage && window.localStorage.getItem('mlip.wsDebug') === '1')
          window.__MLIPVIEW_DEBUG_WS = true;
      } catch { }
      // Expose a global toggle for convenience
      try {
        if (typeof window.__WS_DEBUG_ENABLE__ !== 'function') {
          window.__WS_DEBUG_ENABLE__ = (on) => {
            try {
              window.__MLIPVIEW_DEBUG_WS = !!on;
              if (on && window.localStorage) window.localStorage.setItem('mlip.wsDebug', '1');
              else if (window.localStorage) window.localStorage.removeItem('mlip.wsDebug');
              console.log('[WS] debug set to', !!on);
            } catch { }
          };
        }
      } catch { }
      return !!(window.__MLIPVIEW_DEBUG_WS || window.__MLIPVIEW_DEBUG_API);
    }
  } catch { }
  try {
    return !!(typeof process !== 'undefined' && process.env && process.env.JEST_WORKER_ID);
  } catch { }
  return false;
}
function __wsLog(...a) {
  try {
    if (__wsDebugOn()) console.log(...a);
  } catch { }
}
function __wsWarn(...a) {
  try {
    if (__wsDebugOn()) console.warn(...a);
  } catch { }
}
function __wsErr(...a) {
  try {
    console.error(...a);
  } catch { }
}

function resolveWsBase() {
  if (typeof window === 'undefined') return 'ws://127.0.0.1:8000';
  const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
  const host = window.__MLIPVIEW_SERVER
    ? new URL(window.__MLIPVIEW_SERVER).host
    : window.location.host;
  return `${proto}//${host}`;
}

let __singleton = null;
let __connectPromise = null; // guard parallel connects

export function createFairchemWS() {
  __count('ws#createFairchemWS');
  let ws = null;
  let seq = 0;
  let clientAck = 0;
  let lastCounters = { userInteractionCount: 0, simStep: 0 };
  let dragLock = null; // { index, pos }
  // Optional per-instance test hook; prefer this over any global fallbacks
  let __testHook = null;
  const listeners = new Set();
  function onFrame(fn) {
    listeners.add(fn);
    return () => listeners.delete(fn);
  }
  async function connect() {
    await __ensurePB();
    const base = resolveWsBase();
    const url = base.replace(/\/$/, '') + '/ws';
    __wsLog('[WS][connect]', url);
    ws = new WebSocket(url);
    return new Promise((resolve, reject) => {
      ws.binaryType = 'arraybuffer';
      ws.onopen = () => {
        __wsLog('[WS][open]', { readyState: ws?.readyState });
        try {
          if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function')
            window.__WS_ON_STATE__({ type: 'open', url });
        } catch { }
        resolve(true);
      };
      ws.onerror = (e) => {
        __wsErr('[WS][error]', e?.message || e);
        try {
          if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function')
            window.__WS_ON_STATE__({ type: 'error', error: e?.message || String(e) });
        } catch { }
        reject(e);
      };
      ws.onclose = (ev) => {
        __wsWarn('[WS][close]', { code: ev?.code, reason: ev?.reason });
        try {
          if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function')
            window.__WS_ON_STATE__({ type: 'close', code: ev?.code, reason: ev?.reason });
        } catch { }
      };

      // --- small local helpers for this connect() scope ---
      const __bytesFromWSData = (data) => {
        if (data instanceof ArrayBuffer) return new Uint8Array(data);
        if (typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView && ArrayBuffer.isView(data)) {
          return new Uint8Array(data.buffer, data.byteOffset || 0, data.byteLength || 0);
        }
        if (typeof Uint8Array !== 'undefined' && data instanceof Uint8Array) return data;
        if (typeof Buffer !== 'undefined' && data && Buffer.isBuffer(data)) {
          return new Uint8Array(data.buffer, data.byteOffset || 0, data.byteLength || data.length || 0);
        }
        if (data && typeof data.arrayBuffer === 'function') {
          // Blob (browser)
          return data.arrayBuffer().then((ab) => new Uint8Array(ab));
        }
        const t = data && data.constructor && data.constructor.name ? data.constructor.name : typeof data;
        throw new Error('[WS] Expected binary protobuf frame, got: ' + t);
      };
      const __mat3ToRows = (m) =>
        Array.isArray(m) && m.length === 9
          ? [
            [m[0], m[1], m[2]],
            [m[3], m[4], m[5]],
            [m[6], m[7], m[8]],
          ]
          : undefined;

      ws.onmessage = async (ev) => {
        // Protobuf-only decode (but tolerate Blob and JSON-delivered-as-binary)
        try {
          // ✅ New: if the browser delivered a TEXT frame, ignore it explicitly.
          if (typeof ev.data === 'string') {
            __wsErr('[WS] Received TEXT frame; protobuf-only server. Ignoring.');
            return;
          }

          // Normalize to Uint8Array (handles ArrayBuffer/Blob/Buffer/typed array views)
          const maybeBytes = __bytesFromWSData(ev.data);
          const bytes = maybeBytes instanceof Promise ? await maybeBytes : maybeBytes;


          const r = fromBinary(ServerResultSchema, bytes);
          const out = {
            seq: __n(r.seq),
            client_seq: __n(r.clientSeq),
            userInteractionCount: __n(r.userInteractionCount),
            simStep: __n(r.simStep),
          };

          const which = r.payload?.case;
          if (which === 'frame') {
            const fr = r.payload.value;
            if (fr.energy != null) out.energy = fr.energy;

            if (Array.isArray(fr.positions)) out.positions = __flat3ToTriples(fr.positions) || [];
            if (Array.isArray(fr.forces)) out.forces = __flat3ToTriples(fr.forces) || [];
            if (Array.isArray(fr.velocities)) out.velocities = __flat3ToTriples(fr.velocities) || [];

            if (fr.cell && Array.isArray(fr.cell.m)) out.cell = __mat3ToRows(fr.cell.m);
            if (fr.stress && Array.isArray(fr.stress.m)) out.stress = __mat3ToRows(fr.stress.m);

            // Apply drag-lock override for the currently dragged atom, if any
            try {
              if (dragLock && out.positions && out.positions[dragLock.index]) {
                out.positions[dragLock.index] = dragLock.pos || out.positions[dragLock.index];
              }
            } catch { }
          } else if (which === 'notice') {
            const n = r.payload.value;
            if (typeof n.message === 'string') {
              // Normalize variants like "WAITING_FOR_ACKS" or "WAITING_FOR_ACK seq=123"
              const m = n.message.toUpperCase();
              out.message = m.includes('WAITING_FOR_ACK') ? 'WAITING_FOR_ACK' : n.message;
            }
            if (typeof n.simulationStopped === 'boolean') out.simulationStopped = !!n.simulationStopped;
          } else {
            // Unknown payload -> ignore
            return;
          }

          __wsLog('[WS][rx]', {
            seq: out.seq,
            client_seq: out.client_seq,
            uic: out.userInteractionCount,
            simStep: out.simStep,
            have: {
              pos: !!out.positions,
              forces: !!out.forces,
              vel: !!out.velocities,
              cell: !!out.cell,
              energy: typeof out.energy === 'number',
            },
            simulationStopped: !!out.simulationStopped,
          });

          try {
            const dbgApi = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API;
            if ((dbgApi || __wsDebugOn()) && Array.isArray(out.positions)) {
              console.log('[WS][rx][positions]', out.positions);
            }
          } catch { }

          for (const fn of listeners) {
            try { fn(out, lastCounters); } catch { }
          }
        } catch (e) {
          try { __wsErr('[WS][decode-error]', e?.message || e); } catch { }
        }
      };
    });
  }

  function waitForEnergy({ timeoutMs = 10000 } = {}) {
    return new Promise((resolve, reject) => {
      let off = null;
      let done = false;
      let to = null;
      function finish(v) {
        if (done) return;
        done = true;
        try {
          off && off();
        } catch { }
        if (to)
          try {
            clearTimeout(to);
          } catch { }
        resolve(v);
      }
      function fail(e) {
        if (done) return;
        done = true;
        try {
          off && off();
        } catch { }
        if (to)
          try {
            clearTimeout(to);
          } catch { }
        reject(e);
      }
      off = onFrame((res) => {
        try {
          if (!res || typeof res !== 'object') return;
          const energy = res.energy != null ? res.energy : undefined;
          const forces = Array.isArray(res.forces) ? res.forces : [];
          // Resolve on the first frame carrying energy (idle compute or simulation)
          if (energy != null) {
            try {
              if (typeof res.seq === 'number') ack(res.seq);
            } catch { }
            finish({
              energy,
              forces,
              userInteractionCount: res.userInteractionCount | 0,
              simStep: res.simStep | 0,
            });
          }
        } catch { }
      });
      to = setTimeout(() => fail(new Error('waitForEnergy timeout')), timeoutMs | 0);
    });
  }
  function sendBytes(buf) {
    if (ws && ws.readyState === 1) {
      try {
        __wsLog('[WS][tx-bytes]', { len: buf?.byteLength || buf?.length || 0 });
      } catch { }
      try {
        ws.send(buf); // << sends protobuf bytes only; never JSON
      } catch (e) {
        __wsErr('[WS][send-error]', e?.message || e);
      }
    } else {
      __wsWarn('[WS][send-skipped:not-open]', { readyState: ws?.readyState });
    }
  }
  function close() {
    try {
      ws && ws.close();
    } catch { }
  }
  function setCounters({ userInteractionCount, simStep }) {
    if (Number.isFinite(userInteractionCount))
      lastCounters.userInteractionCount = userInteractionCount | 0;
    if (Number.isFinite(simStep)) lastCounters.simStep = simStep | 0;
  }
  function nextSeq() {
    seq = (seq | 0) + 1;
    return seq;
  }
  // Test hook: prefer an instance-level hook set via setTestHook(fn),
  // then fall back to a single global location on globalThis for legacy tests.
  function __notifyTestHook(msg, explicitKind, meta) {
    try {
      const sinks = [];
      try {
        if (__testHook && typeof __testHook === 'function') sinks.push(__testHook);
      } catch { }
      try {
        if (
          typeof globalThis !== 'undefined' &&
          typeof globalThis.__WS_TEST_HOOK__ === 'function'
        )
          sinks.push(globalThis.__WS_TEST_HOOK__);
      } catch { }
      try {
        if (typeof window !== 'undefined' && typeof window.__WS_TEST_HOOK__ === 'function')
          sinks.push(window.__WS_TEST_HOOK__);
      } catch { }
      if (!sinks.length) return;
      // Project a minimal JSON snapshot (type, simulation params, counters)
      const kind =
        explicitKind ||
        (msg && msg.userInteraction
          ? 'USER_INTERACTION'
          : msg && msg.start
            ? 'START_SIMULATION'
            : msg && msg.stop
              ? 'STOP_SIMULATION'
              : msg && msg.ping
                ? 'PING'
                : undefined);
      const out = {
        seq: msg.seq | 0,
        type: kind,
        userInteractionCount: msg.userInteractionCount | 0,
        simStep: msg.simStep | 0,
        // Include simulation params when present
        simulationType: meta?.simulationType ?? msg.start?.simulationType,
        simulationParams:
          meta?.simulationParams ??
          (msg.start?.simulationParams
            ? {
              calculator: msg.start.simulationParams.calculator || '',
              temperature: msg.start.simulationParams.temperature,
              timestepFs: msg.start.simulationParams.timestepFs,
              friction: msg.start.simulationParams.friction,
              fmax: msg.start.simulationParams.fmax,
              maxStep: msg.start.simulationParams.maxStep,
              optimizer: msg.start.simulationParams.optimizer || undefined,
            }
            : undefined),
        positionsCount:
          meta?.positionsLenTriples ??
          (Array.isArray(msg.userInteraction?.positions)
            ? Math.floor(msg.userInteraction.positions.length / 3)
            : undefined),
        velocitiesCount:
          meta?.velocitiesLenTriples ??
          (Array.isArray(msg.userInteraction?.velocities)
            ? Math.floor(msg.userInteraction.velocities.length / 3)
            : undefined),
        hasCell: meta?.hasCell ?? !!msg.userInteraction?.cell,
      };
      try {
        const dbg = (typeof process !== 'undefined' && process.env && process.env.JEST_WORKER_ID) ||
          (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_WS);
        if (dbg) console.log('[WS][TEST_HOOK]', out);
      } catch { }
      for (const fn of sinks) {
        try {
          fn(out);
        } catch { }
      }
    } catch { }
  }

  // Allow tests to inject a hook without using globals
  function setTestHook(fn) {
    try {
      __testHook = typeof fn === 'function' ? fn : null;
    } catch {
      __testHook = null;
    }
  }
  function setAck(s) {
    clientAck = Math.max(clientAck, s | 0);
  }
  function getState() {
    return { seq, clientAck, ...lastCounters, connected: !!ws && ws.readyState === 1 };
  }
  async function ensureConnected() {
    if (ws && ws.readyState === 1) return true;
    if (__connectPromise) {
      try {
        await __connectPromise;
        return true;
      } catch {
        /* fallthrough */
      }
    }
    __connectPromise = connect().finally(() => {
      __connectPromise = null;
    });
    await __connectPromise;
    return true;
  }

  // INIT removed: delegate to USER_INTERACTION for initialization.
  function initSystem({ atomic_numbers, positions, velocities, cell }) {
    __wsWarn('[WS] initSystem deprecated; sending USER_INTERACTION for init');
    return userInteraction({ atomic_numbers, positions, velocities, cell });
  }

  // USER_INTERACTION: allow partial state updates (positions, velocities, or cell)
  function userInteraction({ atomic_numbers, positions, velocities, cell, dragLockIndex } = {}) {
    const payload = {};
    if (Array.isArray(atomic_numbers)) payload.atomicNumbers = atomic_numbers.map((z) => z | 0);
    if (Array.isArray(positions)) payload.positions = __triplesToFlat3(positions);
    if (Array.isArray(velocities)) payload.velocities = __triplesToFlat3(velocities);
    if (cell && Array.isArray(cell) && cell.length === 3)
      payload.cell = {
        m: [
          +cell[0][0],
          +cell[0][1],
          +cell[0][2],
          +cell[1][0],
          +cell[1][1],
          +cell[1][2],
          +cell[2][0],
          +cell[2][1],
          +cell[2][2],
        ],
      };
    // Build message with counters set at creation time to ensure encoder includes them
    const msg = create(ClientActionSchema, {
      seq: nextSeq(),
      schemaVersion: 1,
      userInteractionCount: Number.isFinite(lastCounters.userInteractionCount)
        ? (lastCounters.userInteractionCount | 0)
        : undefined,
      simStep: Number.isFinite(lastCounters.simStep) ? (lastCounters.simStep | 0) : undefined,
      payload: { case: 'userInteraction', value: payload },
    });
    __notifyTestHook(msg, 'USER_INTERACTION', {
      positionsLenTriples: Array.isArray(positions) ? positions.length : undefined,
      velocitiesLenTriples: Array.isArray(velocities) ? velocities.length : undefined,
      hasCell: !!payload.cell,
    });
    try {
      __wsLog('[WS][tx][USER_INTERACTION]', {
        seq: msg.seq | 0,
        uic: msg.userInteractionCount | 0,
        atoms: Array.isArray(payload.atomicNumbers) ? payload.atomicNumbers.length : 0,
        positions: Array.isArray(payload.positions) ? payload.positions.length : 0,
        velocities: Array.isArray(payload.velocities) ? payload.velocities.length : 0,
        hasCell: !!msg.userInteraction?.cell,
      });
      // When debug is enabled via ?debug=1 or wsDebug=1, print the positions being sent
      try {
        const dbgApi = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API;
        if ((dbgApi || __wsDebugOn()) && Array.isArray(positions)) {
          console.log('[WS][tx][USER_INTERACTION][positions]', positions);
          console.log('[WS][tx][USER_INTERACTION][flat3.len]',
            Array.isArray(payload.positions) ? payload.positions.length : -1);
        }
      } catch { }
    } catch { }
    const buf = toBinary(ClientActionSchema, msg);
    sendBytes(buf);
    try {
      if (Number.isInteger(dragLockIndex) && Array.isArray(positions)) {
        dragLock = { index: dragLockIndex | 0, pos: positions[dragLockIndex | 0] };
      }
    } catch { }
  }

  function beginDrag(index) {
    try {
      dragLock = { index: index | 0, pos: null };
    } catch { }
  }
  function endDrag() {
    try {
      dragLock = null;
    } catch { }
  }

  // Subscribe to decoded ServerResult frames (convenience over raw onFrame)
  function onResult(fn) {
    // Directly invoke the provided callback with decoded frames.
    // Do NOT re-broadcast via window.__ON_WS_RESULT__ here to avoid recursion,
    // since window.__ON_WS_RESULT__ itself fans out to these listeners in tests.
    return onFrame((decoded) => {
      try {
        if (!decoded || typeof decoded !== 'object') return;
        fn(decoded);
      } catch { }
    });
  }

  // Start simulation (md or relax) with SimulationParams
  function startSimulation({ type, params }) {
    const t =
      String(type || 'md').toLowerCase() === 'md'
        ? ClientAction_Start_SimType.MD
        : ClientAction_Start_SimType.RELAX;
    const sp = params
      ? {
        calculator: params.calculator || '',
        temperature: typeof params.temperature === 'number' ? +params.temperature : undefined,
        timestepFs: typeof params.timestep_fs === 'number' ? +params.timestep_fs : undefined,
        friction: typeof params.friction === 'number' ? +params.friction : undefined,
        fmax: typeof params.fmax === 'number' ? +params.fmax : undefined,
        maxStep: typeof params.max_step === 'number' ? +params.max_step : undefined,
        optimizer: params.optimizer || undefined,
      }
      : undefined;
    const msg = create(ClientActionSchema, {
      seq: nextSeq(),
      schemaVersion: 1,
      payload: {
        case: 'start',
        value: {
          simulationType: t,
          simulationParams: sp,
        },
      },
      userInteractionCount: lastCounters.userInteractionCount | 0,
      simStep: lastCounters.simStep | 0,
    });
    __notifyTestHook(msg, 'START_SIMULATION', {
      simulationType: t,
      simulationParams: sp
        ? {
          calculator: sp.calculator || '',
          temperature: sp.temperature,
          timestepFs: sp.timestepFs,
          friction: sp.friction,
          fmax: sp.fmax,
          maxStep: sp.maxStep,
          optimizer: sp.optimizer || undefined,
        }
        : undefined,
    });
    try {
      __wsLog('[WS][tx][START_SIMULATION]', {
        seq: msg.seq | 0,
        simType: t === ClientAction_Start_SimType.MD ? 'MD' : 'RELAX',
        uic: msg.userInteractionCount | 0,
        simStep: msg.simStep | 0,
        params: {
          calculator: sp?.calculator,
          temperature: sp?.temperature,
          timestepFs: sp?.timestepFs,
          friction: sp?.friction,
          fmax: sp?.fmax,
          maxStep: sp?.maxStep,
          optimizer: sp?.optimizer,
        },
      });
    } catch { }
    const buf = toBinary(ClientActionSchema, msg);
    sendBytes(buf);
  }

  function stopSimulation() {
    const msg = create(ClientActionSchema, {
      seq: nextSeq(),
      schemaVersion: 1,
      payload: { case: 'stop', value: {} },
    });
    __notifyTestHook(msg, 'STOP_SIMULATION');
    try {
      __wsLog('[WS][tx][STOP_SIMULATION]', { seq: msg.seq | 0 });
    } catch { }
    const buf = toBinary(ClientActionSchema, msg);
    sendBytes(buf);
  }

  function ack(seqNum) {
    try {
      const msg = create(ClientActionSchema, {
        seq: nextSeq(),
        schemaVersion: 1,
        ack: seqNum | 0,
        payload: { case: 'ping', value: {} },
      });
      __notifyTestHook(msg, 'PING');
      const buf = toBinary(ClientActionSchema, msg);
      sendBytes(buf);
      __wsLog('[WS][tx][ACK]', { seq: msg.seq | 0, ack: seqNum | 0 });
      setAck(seqNum | 0);
    } catch (e) {
      try {
        console.error('[WS][ack-error]', e?.message || e);
      } catch { }
    }
  }

  // One-shot request for a single simulation step (md or relax). Starts, waits 1 frame, stops, resolves.
  async function requestSingleStep({ type, params }) {
    if (!ws || ws.readyState !== 1) {
      __wsWarn('[WS][single-step][not-connected]');
      throw new Error('ws not connected');
    }
    return new Promise((resolve, reject) => {
      let off = null;
      let doneCalled = false;
      function done(v) {
        if (doneCalled) return;
        doneCalled = true;
        try {
          off && off();
        } catch { }
        try {
          stopSimulation();
        } catch { }
        resolve(v);
      }
      off = onFrame((res) => {
        try {
          if (res && typeof res === 'object') {
            __wsLog('[WS][rx][single-step]', {
              hasPos: !!res.positions,
              hasForces: !!res.forces,
              energy: res.energy,
            });
            try {
              if (typeof res.seq === 'number') ack(res.seq);
            } catch { }
            done(res);
          }
        } catch { }
      });
      try {
        startSimulation({ type, params });
      } catch (e) {
        try {
          off && off();
        } catch { }
        reject(e);
      }
      setTimeout(() => {
        if (!doneCalled) {
          __wsWarn('[WS][single-step][timeout] after 15s');
          try {
            off && off();
          } catch { }
          reject(new Error('single_step timeout'));
        }
      }, 15000);
    });
  }
  // Allow tests to inject decoded frames directly into listeners, avoiding globals
  function injectTestResult(obj) {
    try {
      for (const fn of listeners) {
        try {
          fn(obj, lastCounters);
        } catch { }
      }
    } catch { }
  }
  const api = {
    connect,
    ensureConnected,
    sendBytes,
    close,
    onFrame,
    onResult,
    setCounters,
    nextSeq,
    setAck,
    ack,
    getState,
    requestSingleStep,
    waitForEnergy,
    initSystem,
    userInteraction,
    beginDrag,
    endDrag,
    startSimulation,
    stopSimulation,
    setTestHook,
    injectTestResult,
  };
  try {
    if (typeof window !== 'undefined') {
      window.__fairchem_ws__ = api;
      // Toggle WS debug at runtime from tests or devtools
      if (typeof window.__WS_DEBUG_ENABLE__ !== 'function') {
        window.__WS_DEBUG_ENABLE__ = (on) => {
          try {
            window.__MLIPVIEW_DEBUG_WS = !!on;
            if (on && window.localStorage) window.localStorage.setItem('mlip.wsDebug', '1');
            else if (window.localStorage) window.localStorage.removeItem('mlip.wsDebug');
            console.log('[WS] debug set to', !!on);
          } catch { }
        };
      }
      // Provide a fan-out for tests running in jsdom that call window.__ON_WS_RESULT__
      if (typeof window.__ON_WS_RESULT__ !== 'function') {
        window.__ON_WS_RESULT__ = (obj) => {
          try {
            for (const fn of listeners) {
              try {
                fn(obj, lastCounters);
              } catch { }
            }
          } catch { }
        };
      }
    }
    // Also expose test hook in Node/Jest environments without window
    try {
      if (typeof globalThis !== 'undefined') {
        globalThis.__fairchem_ws__ = globalThis.__fairchem_ws__ || api;
        if (typeof globalThis.__ON_WS_RESULT__ !== 'function') {
          globalThis.__ON_WS_RESULT__ = (obj) => {
            try {
              for (const fn of listeners) {
                try {
                  fn(obj, lastCounters);
                } catch { }
              }
            } catch { }
          };
        }
      }
    } catch { }
  } catch { }
  return api;
}

export function getWS() {
  if (!__singleton) __singleton = createFairchemWS();
  return __singleton;
}

export default { createFairchemWS };
