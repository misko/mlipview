// fairchem_ws_client.js — thin WebSocket client for fairchem_local_server2
// Refactored: single subscribe API, batched ACKs, no dragLock mutation.

import { __count } from './util/funcCount.js';
import { create, toBinary, fromBinary } from '@bufbuild/protobuf';
import {
  ClientActionSchema,
  ClientAction_Start_SimType,
  ServerResultSchema,
} from '/proto/fairchem_local_server2/session_pb.js';

const __hasPB = !!ClientActionSchema && !!ServerResultSchema;
function __assertPB() { if (!__hasPB) throw new Error('[WS][protobuf-missing] ESM stubs not found'); }

let __DEBUG_CACHE = null;
function __wsDebugOn() {
  if (__DEBUG_CACHE != null) return __DEBUG_CACHE;
  let on = false;
  try {
    if (typeof window !== 'undefined') {
      const q = new URLSearchParams(window.location.search || '');
      on =
        q.get('wsDebug') === '1' ||
        q.get('wsDebug') === 'true' ||
        q.get('debug') === '1' ||
        q.get('debug') === 'true' ||
        (window.localStorage && window.localStorage.getItem('mlip.wsDebug') === '1') ||
        !!window.__MLIPVIEW_DEBUG_WS ||
        !!window.__MLIPVIEW_DEBUG_API;
    }
  } catch { }
  try { on = on || !!(typeof process !== 'undefined' && process.env && process.env.JEST_WORKER_ID); } catch { }
  __DEBUG_CACHE = !!on;
  return __DEBUG_CACHE;
}
const __log = (...a) => { if (__wsDebugOn()) try { console.log(...a); } catch { } };
const __warn = (...a) => { if (__wsDebugOn()) try { console.warn(...a); } catch { } };
const __err = (...a) => { try { console.error(...a); } catch { } };

function __n(v) { return typeof v === 'bigint' ? Number(v) : typeof v === 'number' ? v : 0; }

function __triplesToFlat3(triples) {
  if (!Array.isArray(triples)) return undefined;
  const out = new Array(triples.length * 3);
  let j = 0;
  for (let i = 0; i < triples.length; i++) {
    const p = triples[i] || [0, 0, 0];
    out[j++] = +p[0]; out[j++] = +p[1]; out[j++] = +p[2];
  }
  return out;
}

function __flat3ToTriples(arr) {
  if (!Array.isArray(arr)) return undefined;
  const n = Math.floor(arr.length / 3);
  const out = new Array(n);
  for (let i = 0, k = 0; i < n; i++) out[i] = [arr[k++], arr[k++], arr[k++]];
  return out;
}

function __mat3ToRows(m) {
  return Array.isArray(m) && m.length === 9
    ? [[m[0], m[1], m[2]], [m[3], m[4], m[5]], [m[6], m[7], m[8]]]
    : undefined;
}

function resolveWsBase() {
  if (typeof window === 'undefined') return 'ws://127.0.0.1:8000';
  const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
  const host = window.__MLIPVIEW_SERVER ? new URL(window.__MLIPVIEW_SERVER).host : window.location.host;
  return `${proto}//${host}`;
}

async function toBytes(data) {
  if (data instanceof ArrayBuffer) return new Uint8Array(data);
  if (typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView && ArrayBuffer.isView(data)) {
    return new Uint8Array(data.buffer, data.byteOffset || 0, data.byteLength || 0);
  }
  if (typeof Uint8Array !== 'undefined' && data instanceof Uint8Array) return data;
  if (typeof Buffer !== 'undefined' && Buffer.isBuffer?.(data)) {
    return new Uint8Array(data.buffer, data.byteOffset, data.byteLength);
  }
  if (data && typeof data.arrayBuffer === 'function') {
    const ab = await data.arrayBuffer();
    return new Uint8Array(ab);
  }
  const t = data && data.constructor && data.constructor.name ? data.constructor.name : typeof data;
  throw new Error('[WS] Expected binary protobuf frame, got: ' + t);
}

/* -------------------------------------------------------------------------------------------------
 * Singleton
 * ------------------------------------------------------------------------------------------------- */

let __singleton = null;
let __connectPromise = null;

/* -------------------------------------------------------------------------------------------------
 * Factory
 * ------------------------------------------------------------------------------------------------- */

/**
 * @returns {{
 *  connect: () => Promise<boolean>,
 *  ensureConnected: () => Promise<boolean>,
 *  sendBytes: (buf: Uint8Array) => void,
 *  close: () => void,
 *  onFrame: (fn: (res: any, counters: any)=>void) => () => void,
 *  onResult: (fn: (res: any)=>void) => () => void,
 *  setCounters: ({userInteractionCount?:number, simStep?:number}) => void,
 *  nextSeq: () => number,
 *  setAck: (s:number) => void,
 *  ack: (seq:number) => void,
 *  getState: () => { seq:number, clientAck:number, userInteractionCount:number, simStep:number, connected:boolean },
 *  requestSingleStep: ({type:'md'|'relax', params?:any}) => Promise<any>,
 *  waitForEnergy: ({timeoutMs?:number}?) => Promise<any>,
 *  initSystem: (args:any)=>void,
 *  userInteraction: (args:any)=>void,
 *  beginDrag: (index:number)=>void,
 *  endDrag: ()=>void,
 *  startSimulation: ({type:'md'|'relax', params?:any})=>void,
 *  stopSimulation: ()=>void,
 *  setTestHook: (fn:Function)=>void,
 *  injectTestResult: (obj:any)=>void,
 * }}
 */
export function createFairchemWS() {
  __count('ws#createFairchemWS');
  __assertPB();

  /** @type {WebSocket|null} */
  let ws = null;
  let seq = 0;
  let clientAck = 0;
  const lastCounters = { userInteractionCount: 0, simStep: 0 };

  // Instance-level test hook
  let __testHook = null;

  /** @type {Set<Function>} */
  const listeners = new Set();

  function subscribe(fn) {
    listeners.add(fn);
    return () => listeners.delete(fn);
  }
  // Compat aliases
  const onFrame = (fn) => subscribe((res) => { try { fn(res, lastCounters); } catch { } });
  const onResult = (fn) => subscribe((res) => { try { fn(res); } catch { } });

  function fanout(obj) {
    for (const fn of listeners) {
      try { fn(obj); } catch { }
    }
  }

  function nextSeq() { seq = (seq | 0) + 1; return seq; }
  function setAck(s) { clientAck = Math.max(clientAck, s | 0); }
  function getState() { return { seq, clientAck, ...lastCounters, connected: !!ws && ws.readyState === 1 }; }
  function setCounters({ userInteractionCount, simStep }) {
    if (Number.isFinite(userInteractionCount)) lastCounters.userInteractionCount = userInteractionCount | 0;
    if (Number.isFinite(simStep)) lastCounters.simStep = simStep | 0;
  }

  function __notifyTestHook(msg, explicitKind, meta) {
    try {
      const sinks = [];
      if (typeof __testHook === 'function') sinks.push(__testHook);
      if (typeof globalThis !== 'undefined' && typeof globalThis.__WS_TEST_HOOK__ === 'function') sinks.push(globalThis.__WS_TEST_HOOK__);
      if (typeof window !== 'undefined' && typeof window.__WS_TEST_HOOK__ === 'function') sinks.push(window.__WS_TEST_HOOK__);
      if (!sinks.length) return;

      const kind =
        explicitKind ||
        (msg && msg.userInteraction ? 'USER_INTERACTION'
          : msg && msg.start ? 'START_SIMULATION'
            : msg && msg.stop ? 'STOP_SIMULATION'
              : msg && msg.ping ? 'PING' : undefined);

      const out = {
        seq: msg.seq | 0,
        type: kind,
        userInteractionCount: msg.userInteractionCount | 0,
        simStep: msg.simStep | 0,
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
      if (__wsDebugOn()) console.log('[WS][TEST_HOOK]', out);
      for (const fn of sinks) { try { fn(out); } catch { } }
    } catch { }
  }

  function sendBytes(buf) {
    if (ws && ws.readyState === 1) {
      __log('[WS][tx-bytes]', { len: buf?.byteLength || buf?.length || 0 });
      try { ws.send(buf); } catch (e) { __err('[WS][send-error]', e?.message || e); }
    } else {
      __warn('[WS][send-skipped:not-open]', { readyState: ws?.readyState });
    }
  }

  let wantConnected = false;

  function close() { wantConnected = false; try { ws && ws.close(); } catch { } }

  async function connect() {
    __assertPB();
    const base = resolveWsBase();
    const url = base.replace(/\/$/, '') + '/ws';
    __log('[WS][connect]', url);

    ws = new WebSocket(url);
    ws.binaryType = 'arraybuffer';

    return new Promise((resolve, reject) => {
      ws.onopen = () => {
        wantConnected = true;
        __log('[WS][open]', { readyState: ws?.readyState });
        try { if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function') window.__WS_ON_STATE__({ type: 'open', url }); } catch { }
        resolve(true);
      };
      ws.onerror = (e) => {
        __err('[WS][error]', e?.message || e);
        try { if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function') window.__WS_ON_STATE__({ type: 'error', error: e?.message || String(e) }); } catch { }
        reject(e);
      };
      ws.onclose = (ev) => {
        __warn('[WS][close]', { code: ev?.code, reason: ev?.reason });
        try { if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__ === 'function') window.__WS_ON_STATE__({ type: 'close', code: ev?.code, reason: ev?.reason }); } catch { }
        if (wantConnected) setTimeout(() => { if (!__connectPromise) ensureConnected().catch(() => { }); }, 1000);
      };
      ws.onmessage = async (ev) => {
        try {
          if (typeof ev.data === 'string') { __err('[WS] Received TEXT frame; protobuf-only server. Ignoring.'); return; }
          const bytes = await toBytes(ev.data);
          let i = 0; while (i < bytes.length && bytes[i] <= 32) i++;
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
            if (typeof fr.temperature === 'number') out.temperature = fr.temperature;
            if (typeof fr.kinetic === 'number') out.kinetic = fr.kinetic;
            if (Array.isArray(fr.positions)) out.positions = __flat3ToTriples(fr.positions) || [];
            if (Array.isArray(fr.forces)) out.forces = __flat3ToTriples(fr.forces) || [];
            if (Array.isArray(fr.velocities)) out.velocities = __flat3ToTriples(fr.velocities) || [];
            if (fr.cell && Array.isArray(fr.cell.m)) out.cell = __mat3ToRows(fr.cell.m);
            if (fr.stress && Array.isArray(fr.stress.m)) out.stress = __mat3ToRows(fr.stress.m);
          } else if (which === 'notice') {
            const n = r.payload.value;
            if (typeof n.message === 'string') {
              const upper = n.message.toUpperCase();
              out.message = upper.includes('WAITING_FOR_ACK') ? 'WAITING_FOR_ACK' : n.message;
            }
            if (typeof n.simulationStopped === 'boolean') out.simulationStopped = !!n.simulationStopped;
            console.log('[WS][notice]', out.message || '(no message)', { simulation_stopped: !!out.simulationStopped, seq: out.seq });
          } else {
            __err('[WS] Unknown payload case:', which);
            return;
          }

          __log('[WS][rx]', {
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
            message: out.message || null,
          });

          try {
            const dbgApi = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API;
            if ((dbgApi || __wsDebugOn()) && Array.isArray(out.positions)) {
              console.log('[WS][rx][positions]', out.positions);
              console.log('[WS][rx][velocities]', out.velocities);
            }
          } catch { }

          fanout(out);
          try { if (typeof out.seq === 'number') maybeAck(out.seq); } catch { }
        } catch (e) {
          __err('[WS] onmessage decode error:', e?.message || String(e));
        }
      };
    });
  }

  async function ensureConnected() {
    if (ws && ws.readyState === 1) return true;
    if (__connectPromise) { await __connectPromise; return true; }
    __connectPromise = connect().finally(() => { __connectPromise = null; });
    await __connectPromise;
    return true;
  }

  function withCounters(payload) {
    return create(ClientActionSchema, {
      seq: nextSeq(),
      schemaVersion: 1,
      userInteractionCount: Number.isFinite(lastCounters.userInteractionCount) ? (lastCounters.userInteractionCount | 0) : undefined,
      simStep: Number.isFinite(lastCounters.simStep) ? (lastCounters.simStep | 0) : undefined,
      ...payload,
    });
  }

  function initSystem({ atomic_numbers, positions, velocities, cell }) {
    __warn('[WS] initSystem deprecated; sending USER_INTERACTION for init');
    return userInteraction({ atomic_numbers, positions, velocities, cell });
  }

  function userInteraction({ atomic_numbers, positions, velocities, cell } = {}) {
    const payload = {};
    if (Array.isArray(atomic_numbers)) payload.atomicNumbers = atomic_numbers.map((z) => z | 0);
    if (Array.isArray(positions)) payload.positions = __triplesToFlat3(positions);
    if (Array.isArray(velocities)) payload.velocities = __triplesToFlat3(velocities);
    if (cell && Array.isArray(cell) && cell.length === 3) {
      payload.cell = {
        m: [
          +cell[0][0], +cell[0][1], +cell[0][2],
          +cell[1][0], +cell[1][1], +cell[1][2],
          +cell[2][0], +cell[2][1], +cell[2][2],
        ]
      };
    }

    const msg = withCounters({ payload: { case: 'userInteraction', value: payload } });

    __notifyTestHook(msg, 'USER_INTERACTION', {
      positionsLenTriples: Array.isArray(positions) ? positions.length : undefined,
      velocitiesLenTriples: Array.isArray(velocities) ? velocities.length : undefined,
      hasCell: !!payload.cell,
    });

    __log('[WS][tx][USER_INTERACTION]', {
      seq: msg.seq | 0,
      uic: msg.userInteractionCount | 0,
      atoms: Array.isArray(payload.atomicNumbers) ? payload.atomicNumbers.length : 0,
      positions: Array.isArray(payload.positions) ? payload.positions.length : 0,
      velocities: Array.isArray(payload.velocities) ? payload.velocities.length : 0,
      hasCell: !!msg.userInteraction?.cell,
    });

    try {
      const dbgApi = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_API;
      if ((dbgApi || __wsDebugOn()) && Array.isArray(positions)) {
        console.log('[WS][tx][USER_INTERACTION][positions]', positions);
        console.log('[WS][tx][USER_INTERACTION][flat3.len]',
          Array.isArray(payload.positions) ? payload.positions.length : -1);
      }
    } catch { }

    sendBytes(toBinary(ClientActionSchema, msg));
  }

  // No-ops retained for compatibility
  function beginDrag(_) { }
  function endDrag() { }

  function startSimulation({ type, params }) {
    const t = String(type || 'md').toLowerCase() === 'md'
      ? ClientAction_Start_SimType.MD
      : ClientAction_Start_SimType.RELAX;

    const sp = params ? {
      calculator: params.calculator || 'uma',
      temperature: typeof params.temperature === 'number' ? +params.temperature : undefined,
      timestepFs: typeof params.timestep_fs === 'number' ? +params.timestep_fs : undefined,
      friction: typeof params.friction === 'number' ? +params.friction : undefined,
      fmax: typeof params.fmax === 'number' ? +params.fmax : undefined,
      maxStep: typeof params.max_step === 'number' ? +params.max_step : undefined,
      optimizer: params.optimizer || undefined,
    } : undefined;

    const msg = withCounters({ payload: { case: 'start', value: { simulationType: t, simulationParams: sp } } });

    __notifyTestHook(msg, 'START_SIMULATION', {
      simulationType: t,
      simulationParams: sp ? {
        calculator: sp.calculator || '',
        temperature: sp.temperature,
        timestepFs: sp.timestepFs,
        friction: sp.friction,
        fmax: sp.fmax,
        maxStep: sp.maxStep,
        optimizer: sp.optimizer || undefined,
      } : undefined,
    });

    __log('[WS][tx][START_SIMULATION]', {
      seq: msg.seq | 0,
      simType: t === ClientAction_Start_SimType.MD ? 'MD' : 'RELAX',
      uic: msg.userInteractionCount | 0,
      simStep: msg.simStep | 0,
      params: sp,
    });

    sendBytes(toBinary(ClientActionSchema, msg));
  }

  function stopSimulation() {
    const msg = create(ClientActionSchema, {
      seq: nextSeq(),
      schemaVersion: 1,
      payload: { case: 'stop', value: {} },
    });
    __notifyTestHook(msg, 'STOP_SIMULATION');
    __log('[WS][tx][STOP_SIMULATION]', { seq: msg.seq | 0 });
    sendBytes(toBinary(ClientActionSchema, msg));
  }

  // Batched ACK (one per frame)
  let lastAckSent = 0, highestSeq = 0, ackScheduled = false;
  function flushAck() {
    ackScheduled = false;
    const n = highestSeq | 0;
    if (n <= lastAckSent) return;
    const msg = create(ClientActionSchema, { seq: nextSeq(), schemaVersion: 1, ack: n, payload: { case: 'ping', value: {} } });
    __notifyTestHook(msg, 'PING');
    sendBytes(toBinary(ClientActionSchema, msg));
    __log('[WS][tx][ACK]', { seq: msg.seq | 0, ack: n | 0 });
    lastAckSent = n;
    setAck(n);
  }
  function maybeAck(seqNum) {
    highestSeq = Math.max(highestSeq, seqNum | 0);
    if (!ackScheduled) {
      ackScheduled = true;
      (typeof requestAnimationFrame === 'function')
        ? requestAnimationFrame(flushAck)
        : setTimeout(flushAck, 16);
    }
  }
  function ack(seqNum) { maybeAck(seqNum); }

  function waitForEnergy({ timeoutMs = 10000 } = {}) {
    return new Promise((resolve, reject) => {
      let cleared = false;
      const off = subscribe((res) => {
        try {
          if (!res || typeof res !== 'object') return;
          const energy = res.energy != null ? res.energy : undefined;
          const forces = Array.isArray(res.forces) ? res.forces : [];
          if (energy != null) {
            if (typeof res.seq === 'number') ack(res.seq);
            if (!cleared) {
              cleared = true; off(); clearTimeout(to); resolve({
                energy, forces,
                userInteractionCount: res.userInteractionCount | 0,
                simStep: res.simStep | 0,
              });
            }
          }
        } catch { }
      });
      const to = setTimeout(() => {
        if (!cleared) { cleared = true; off(); reject(new Error('waitForEnergy timeout')); }
      }, timeoutMs | 0);
    });
  }

  async function requestSingleStep({ type, params }) {
    if (!ws || ws.readyState !== 1) {
      __warn('[WS][single-step][not-connected]');
      throw new Error('ws not connected');
    }
    return new Promise((resolve, reject) => {
      let done = false;
      const off = subscribe((res) => {
        try {
          if (res && typeof res === 'object') {
            __log('[WS][rx][single-step]', { hasPos: !!res.positions, hasForces: !!res.forces, energy: res.energy });
            if (typeof res.seq === 'number') ack(res.seq);
            if (!done) { done = true; off(); try { stopSimulation(); } catch { } resolve(res); }
          }
        } catch { }
      });
      try { startSimulation({ type, params }); }
      catch (e) { off(); return reject(e); }

      setTimeout(() => {
        if (!done) { done = true; off(); __warn('[WS][single-step][timeout] after 15s'); reject(new Error('single_step timeout')); }
      }, 15000);
    });
  }

  function setTestHook(fn) { __testHook = typeof fn === 'function' ? fn : null; }
  function injectTestResult(obj) { fanout(obj); }

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
      if (typeof window.__ON_WS_RESULT__ !== 'function') {
        window.__ON_WS_RESULT__ = (obj) => { fanout(obj); };
      }
    }
    if (typeof globalThis !== 'undefined') {
      globalThis.__fairchem_ws__ = globalThis.__fairchem_ws__ || api;
      if (typeof globalThis.__ON_WS_RESULT__ !== 'function') {
        globalThis.__ON_WS_RESULT__ = (obj) => { fanout(obj); };
      }
    }
  } catch { }

  return api;
}

export function getWS() {
  if (!__singleton) __singleton = createFairchemWS();
  return __singleton;
}

export default { createFairchemWS };