// Thin WebSocket client for fairchem_local_server2 /ws protobuf stream
// Import modern ESM protobuf stubs directly
import { __count } from './util/funcCount.js';
import { create, toBinary, fromBinary } from '@bufbuild/protobuf';
import {
		ClientActionSchema,
		ClientAction_Type,
		ClientAction_SimType,
		ServerResultSchema,
		Vec3Schema,
		Mat3Schema,
		SimulationParamsSchema,
} from '/proto/fairchem_local_server2/session_pb.js';
const __hasPB = !!ClientActionSchema && !!ServerResultSchema;
async function __ensurePB(){ if(!__hasPB) throw new Error('[WS][protobuf-missing] ESM stubs not found'); return true; }
function __n(v){ return (typeof v==='bigint') ? Number(v) : (typeof v==='number' ? v : 0); }
function __wsDebugOn(){
	try {
		if (typeof window !== 'undefined') {
			// URL param wsDebug=1|true toggles
			try {
				const q = new URLSearchParams(window.location.search||'');
				if (q.get('wsDebug') === '1' || q.get('wsDebug') === 'true') window.__MLIPVIEW_DEBUG_WS = true;
			} catch {}
			// localStorage key 'mlip.wsDebug' toggles
			try { if (window.localStorage && window.localStorage.getItem('mlip.wsDebug') === '1') window.__MLIPVIEW_DEBUG_WS = true; } catch {}
			return !!(window.__MLIPVIEW_DEBUG_WS || window.__MLIPVIEW_DEBUG_API);
		}
	} catch {}
	try { return !!(typeof process!=='undefined' && process.env && process.env.JEST_WORKER_ID); } catch {}
	return false;
}
function __wsLog(...a){ try { if(__wsDebugOn()) console.log(...a); } catch {} }
function __wsWarn(...a){ try { if(__wsDebugOn()) console.warn(...a); } catch {} }
function __wsErr(...a){ try { console.error(...a); } catch {} }

function resolveWsBase(){
		if (typeof window === 'undefined') return 'ws://127.0.0.1:8000';
	const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
	const host = window.__MLIPVIEW_SERVER ? new URL(window.__MLIPVIEW_SERVER).host : window.location.host;
	return `${proto}//${host}`;
}

let __singleton = null;
let __connectPromise = null; // guard parallel connects

export function createFairchemWS(){
	__count('ws#createFairchemWS');
	let ws = null;
	let seq = 0;
	let clientAck = 0;
	let lastCounters = { userInteractionCount: 0, simStep: 0 };
	let dragLock = null; // { index, pos }
	const listeners = new Set();
	function onFrame(fn){ listeners.add(fn); return ()=>listeners.delete(fn); }

	async function connect(){
		await __ensurePB();
		const base = resolveWsBase();
		const url = base.replace(/\/$/, '') + '/ws';
		__wsLog('[WS][connect]', url);
		ws = new WebSocket(url);
			return new Promise((resolve, reject)=>{
			ws.binaryType = 'arraybuffer';
			ws.onopen = ()=> {
				__wsLog('[WS][open]', { readyState: ws?.readyState });
				try { if (typeof window !== 'undefined' && typeof window.__WS_ON_STATE__==='function') window.__WS_ON_STATE__({ type:'open', url }); } catch {}
				resolve(true);
			};
			ws.onerror = (e)=> { __wsErr('[WS][error]', e?.message||e); try { if (typeof window!=='undefined' && typeof window.__WS_ON_STATE__==='function') window.__WS_ON_STATE__({ type:'error', error: (e?.message||String(e)) }); } catch {} reject(e); };
			ws.onclose = (ev)=>{ __wsWarn('[WS][close]', { code: ev?.code, reason: ev?.reason }); try { if (typeof window!=='undefined' && typeof window.__WS_ON_STATE__==='function') window.__WS_ON_STATE__({ type:'close', code: ev?.code, reason: ev?.reason }); } catch {} };
						ws.onmessage = (ev)=>{
								// Decode protobuf ServerResult frames; also accept JSON
								try {
												// Text JSON path (handle primitive string and String object)
												const isStringLike = (typeof ev.data === 'string') || (!!ev.data && ev.data.constructor && ev.data.constructor.name === 'String');
												if (isStringLike){
													try {
														const s = String(ev.data);
														const obj = JSON.parse(s);
														const out = {
															seq: __n(obj.seq),
															client_seq: __n(obj.client_seq||obj.clientSeq),
															userInteractionCount: __n(obj.userInteractionCount||obj.user_interaction_count),
															simStep: __n(obj.simStep||obj.sim_step),
															message: obj.message,
														};
														if (typeof obj.energy === 'number') out.energy = obj.energy;
														if (Array.isArray(obj.positions)) out.positions = obj.positions;
														if (Array.isArray(obj.forces)) out.forces = obj.forces;
														if (Array.isArray(obj.velocities)) out.velocities = obj.velocities;
														if (Array.isArray(obj.cell)) out.cell = obj.cell;
														__wsLog('[WS][rx]', { seq: out.seq, client_seq: out.client_seq, uic: out.userInteractionCount, simStep: out.simStep, have: { pos: !!out.positions, forces: !!out.forces, vel: !!out.velocities, cell: !!out.cell, energy: typeof out.energy==='number' }, message: out.message });
														try { if (out && out.message === 'WAITING_FOR_ACK') __wsWarn('[WS][rx][WAITING_FOR_ACK]', { seq: out.seq, clientAck }); } catch {}
														for(const fn of listeners){ try { fn(out, lastCounters); } catch{} }
														return;
													} catch(_jerr) {
														// Non-JSON string (e.g., 'pong') — ignore quietly
														try { const s = String(ev.data).trim(); if (s === 'pong') return; } catch {}
														// For any other text that isn't valid JSON, ignore and return (server may send metrics or keepalive)
														return;
													}
												}
												let bytes;
												// Accept ArrayBuffer (browser)
												if (ev.data instanceof ArrayBuffer) {
													bytes = new Uint8Array(ev.data);
												// Accept any ArrayBufferView (Uint8Array, DataView, Buffer-backed typed arrays)
												} else if (typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView && ArrayBuffer.isView(ev.data)) {
													const v = ev.data; bytes = new Uint8Array(v.buffer, v.byteOffset||0, v.byteLength||0);
												// Accept Uint8Array explicitly (redundant with ArrayBufferView but safe)
												} else if (typeof Uint8Array !== 'undefined' && ev.data instanceof Uint8Array) {
													bytes = ev.data;
												// Accept Node Buffer (ws)
												} else if (typeof Buffer !== 'undefined' && ev.data && Buffer.isBuffer(ev.data)) {
													bytes = new Uint8Array(ev.data.buffer, ev.data.byteOffset||0, ev.data.byteLength||ev.data.length||0);
												// Blob fallback (browser)
												} else if (ev.data && ev.data.arrayBuffer && typeof ev.data.arrayBuffer === 'function') {
													return ev.data.arrayBuffer().then(ab=>{ ws.onmessage({ data: ab }); }).catch(()=>{});
												} else {
													const t = (ev && ev.data) ? (ev.data.constructor && ev.data.constructor.name ? ev.data.constructor.name : typeof ev.data) : String(ev && ev.data);
													throw new Error('[WS] Expected binary protobuf frame, got unsupported message type: ' + t);
												}
												// Some Node WS stacks deliver text frames as Buffer/Uint8Array; detect JSON and parse
												try {
													if (bytes && bytes.length) {
														let i=0; while(i<bytes.length && bytes[i]<=32) i++;
														const b0 = bytes[i];
														if (b0===123 /*'{'*/ || b0===91 /*'['*/){
															const td = new TextDecoder('utf-8');
															const s = td.decode(bytes);
															// Re-dispatch to the JSON path
															try { ws.onmessage({ data: s }); return; } catch {}
														}
													}
												} catch {}
									const r = fromBinary(ServerResultSchema, bytes);
									const out = {
										seq: __n(r.seq),
										client_seq: __n(r.clientSeq),
										userInteractionCount: __n(r.userInteractionCount),
										simStep: __n(r.simStep),
										message: r.message,
									};
									if (r.energy != null) out.energy = r.energy;
									if (Array.isArray(r.positions)) out.positions = r.positions.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									// Apply dragLock override for the dragged atom if present
									try {
										if (dragLock && out.positions && out.positions[dragLock.index]) {
											out.positions[dragLock.index] = dragLock.pos || out.positions[dragLock.index];
										}
									} catch{}
									if (Array.isArray(r.forces)) out.forces = r.forces.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									if (Array.isArray(r.velocities)) out.velocities = r.velocities.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									if (r.cell && Array.isArray(r.cell.m)){
										const m = r.cell.m;
										if (m.length===9) out.cell = [ [m[0],m[1],m[2]],[m[3],m[4],m[5]],[m[6],m[7],m[8]] ];
									}
									__wsLog('[WS][rx]', { seq: out.seq, client_seq: out.client_seq, uic: out.userInteractionCount, simStep: out.simStep, have: { pos: !!out.positions, forces: !!out.forces, vel: !!out.velocities, cell: !!out.cell, energy: typeof out.energy==='number' } });
									for(const fn of listeners){ try { fn(out, lastCounters); } catch{} }
								} catch (e) {
						// Surface decode errors to console but do not crash the socket
								try { __wsErr('[WS][decode-error]', e?.message||e); } catch {}
					}
				};
		});
	}
						function waitForEnergy({ timeoutMs=5000 }={}){
			return new Promise((resolve, reject)=>{
				let off = null; let done=false; let to=null;
				function finish(v){ if(done) return; done=true; try{ off&&off(); }catch{} if(to) try{ clearTimeout(to); }catch{} resolve(v); }
				function fail(e){ if(done) return; done=true; try{ off&&off(); }catch{} if(to) try{ clearTimeout(to); }catch{} reject(e); }
				off = onFrame((res)=>{
					try {
						if (!res || typeof res !== 'object') return;
						const energy = (res.energy!=null ? res.energy : undefined);
						const forces = Array.isArray(res.forces)? res.forces : [];
						// Only resolve when we have a real energy-bearing frame or an explicit simple_calculate message.
						// This avoids prematurely resolving on INIT_SYSTEM 'initialized' frames that may carry zero forces.
						const isSimple = (res && typeof res.message === 'string' && res.message.toLowerCase() === 'simple_calculate');
						if (energy != null || isSimple) finish({ energy, forces, userInteractionCount: res.userInteractionCount|0, simStep: res.simStep|0 });
					} catch {}
				});
				to = setTimeout(()=> fail(new Error('waitForEnergy timeout')), timeoutMs|0);
			});
		}
							function sendBytes(buf){ if(ws && ws.readyState===1) { try { __wsLog('[WS][tx-bytes]', { len: buf?.byteLength||buf?.length||0 }); } catch{} try { ws.send(buf); } catch(e){ __wsErr('[WS][send-error]', e?.message||e); } } else { __wsWarn('[WS][send-skipped:not-open]', { readyState: ws?.readyState }); } }
	function close(){ try { ws && ws.close(); } catch{} }
	function setCounters({ userInteractionCount, simStep }){
		if(Number.isFinite(userInteractionCount)) lastCounters.userInteractionCount = userInteractionCount|0;
		if(Number.isFinite(simStep)) lastCounters.simStep = simStep|0;
	}
	function nextSeq(){ seq = (seq|0) + 1; return seq; }
	// Test hook: if window.__WS_TEST_HOOK__ is a function, call with a shallow JSON of the outgoing message
	function __notifyTestHook(msg){
		try {
			const sink = (typeof window !== 'undefined' && typeof window.__WS_TEST_HOOK__ === 'function')
				? window.__WS_TEST_HOOK__
				: (typeof globalThis !== 'undefined' && typeof globalThis.__WS_TEST_HOOK__ === 'function')
					? globalThis.__WS_TEST_HOOK__
					: (typeof global !== 'undefined' && typeof global.__WS_TEST_HOOK__ === 'function')
						? global.__WS_TEST_HOOK__
						: null;
			if (!sink) return;
			// Project a minimal JSON snapshot (type, simulation params, counters)
			const out = {
				seq: msg.seq|0,
				type: msg.type,
				userInteractionCount: msg.userInteractionCount|0,
				simStep: msg.simStep|0,
				// Include simulation params when present
				simulationType: msg.simulationType,
				simulationParams: msg.simulationParams ? {
					calculator: msg.simulationParams.calculator || '',
					temperature: msg.simulationParams.temperature,
					timestepFs: msg.simulationParams.timestepFs,
					friction: msg.simulationParams.friction,
					fmax: msg.simulationParams.fmax,
					maxStep: msg.simulationParams.maxStep,
					optimizer: msg.simulationParams.optimizer || undefined,
				} : undefined,
				positionsCount: Array.isArray(msg.positions) ? msg.positions.length : undefined,
				velocitiesCount: Array.isArray(msg.velocities) ? msg.velocities.length : undefined,
				hasCell: !!msg.cell,
			};
			sink(out);
		} catch {}
	}
	function setAck(s){ clientAck = Math.max(clientAck, s|0); }
	function getState(){ return { seq, clientAck, ...lastCounters, connected: !!ws && ws.readyState===1 }; }
		async function ensureConnected(){
			if (ws && ws.readyState === 1) return true;
			if (__connectPromise) { try { await __connectPromise; return true; } catch { /* fallthrough */ } }
			__connectPromise = connect().finally(()=>{ __connectPromise = null; });
			await __connectPromise; return true;
		}

				// INIT removed: delegate to USER_INTERACTION for initialization.
			function initSystem({ atomic_numbers, positions, velocities, cell }){
				__wsWarn('[WS] initSystem deprecated; sending USER_INTERACTION for init');
				return userInteraction({ atomic_numbers, positions, velocities, cell });
			}

				// USER_INTERACTION: allow partial state updates (positions, velocities, or cell)
				function userInteraction({ atomic_numbers, positions, velocities, cell, dragLockIndex } = {}){
			const out = { seq: nextSeq(), type: ClientAction_Type.USER_INTERACTION };
			if (Array.isArray(atomic_numbers)) out.atomicNumbers = atomic_numbers.map(z=> z|0);
			if (Array.isArray(positions)) out.positions = positions.map(p=> create(Vec3Schema, { v: [+p[0],+p[1],+p[2]] }));
			if (Array.isArray(velocities)) out.velocities = velocities.map(p=> create(Vec3Schema, { v: [+p[0],+p[1],+p[2]] }));
			if (cell && Array.isArray(cell) && cell.length===3) out.cell = create(Mat3Schema, { m: [
				+cell[0][0], +cell[0][1], +cell[0][2],
				+cell[1][0], +cell[1][1], +cell[1][2],
				+cell[2][0], +cell[2][1], +cell[2][2],
			] });
			const msg = create(ClientActionSchema, out);
			__notifyTestHook(msg);
			try { __wsLog('[WS][tx][USER_INTERACTION]', { seq: msg.seq|0, atoms: msg.atomicNumbers?.length|0, positions: msg.positions?.length|0, velocities: msg.velocities?.length|0, hasCell: !!msg.cell }); } catch {}
			const buf = toBinary(ClientActionSchema, msg);
			sendBytes(buf);
					try {
						if (Number.isInteger(dragLockIndex) && Array.isArray(positions)) {
							dragLock = { index: dragLockIndex|0, pos: positions[dragLockIndex|0] };
						}
					} catch {}
		}

				function beginDrag(index){ try { dragLock = { index: index|0, pos: null }; } catch {} }
				function endDrag(){ try { dragLock = null; } catch {} }

				// Subscribe to decoded ServerResult frames (convenience over raw onFrame)
				function onResult(fn){
					// Directly invoke the provided callback with decoded frames.
					// Do NOT re-broadcast via window.__ON_WS_RESULT__ here to avoid recursion,
					// since window.__ON_WS_RESULT__ itself fans out to these listeners in tests.
					return onFrame((decoded)=>{
						try {
							if (!decoded || typeof decoded !== 'object') return;
							fn(decoded);
						} catch {}
					});
				}

			// Start simulation (md or relax) with SimulationParams
			function startSimulation({ type, params }){
				const t = String(type||'md').toLowerCase()==='md' ? ClientAction_SimType.MD : ClientAction_SimType.RELAX;
				const sp = params ? create(SimulationParamsSchema, {
					calculator: params.calculator || '',
					temperature: typeof params.temperature==='number' ? +params.temperature : undefined,
					timestepFs: typeof params.timestep_fs==='number' ? +params.timestep_fs : undefined,
					friction: typeof params.friction==='number' ? +params.friction : undefined,
					fmax: typeof params.fmax==='number' ? +params.fmax : undefined,
					maxStep: typeof params.max_step==='number' ? +params.max_step : undefined,
					optimizer: params.optimizer || undefined,
				}) : undefined;
				const msg = create(ClientActionSchema, {
					seq: nextSeq(),
					type: ClientAction_Type.START_SIMULATION,
					simulationType: t,
					simulationParams: sp,
					userInteractionCount: lastCounters.userInteractionCount|0,
					simStep: lastCounters.simStep|0,
				});
				__notifyTestHook(msg);
				try { __wsLog('[WS][tx][START_SIMULATION]', { seq: msg.seq|0, simType: (t===ClientAction_SimType.MD?'MD':'RELAX'), uic: msg.userInteractionCount|0, simStep: msg.simStep|0, params: { calculator: sp?.calculator, temperature: sp?.temperature, timestepFs: sp?.timestepFs, friction: sp?.friction, fmax: sp?.fmax, maxStep: sp?.maxStep, optimizer: sp?.optimizer } }); } catch {}
				const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
			}

			function stopSimulation(){
				const msg = create(ClientActionSchema, { seq: nextSeq(), type: ClientAction_Type.STOP_SIMULATION });
				__notifyTestHook(msg);
				try { __wsLog('[WS][tx][STOP_SIMULATION]', { seq: msg.seq|0 }); } catch {}
				const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
			}

				function ack(seqNum){
					try {
						const msg = create(ClientActionSchema, {
							seq: nextSeq(),
							type: ClientAction_Type.PING,
							ack: seqNum|0,
						});
						const buf = toBinary(ClientActionSchema, msg);
						sendBytes(buf);
						__wsLog('[WS][tx][ACK]', { seq: msg.seq|0, ack: seqNum|0 });
						setAck(seqNum|0);
					} catch (e) { try { console.error('[WS][ack-error]', e?.message||e); } catch {} }
				}
		// High-level helper: request one-shot simple calculation via WS.
			// Encodes a ClientAction SIMPLE_CALCULATE referencing current session state.
			// Decodes first matching ServerResult with energy/forces and resolves.
			// Returns { energy, forces, userInteractionCount, simStep } so caller can gate staleness.
			async function requestSimpleCalculate(){
			if(!ws || ws.readyState!==1) { __wsWarn('[WS][simple][not-connected]'); throw new Error('ws not connected'); }
			return new Promise((resolve, reject)=>{
				let off = null; let to=null;
				const done = (v)=>{ try{ off && off(); }catch{} if(to) try{ clearTimeout(to); }catch{} resolve(v); };
				const fail = (e)=>{ try{ off && off(); }catch{} if(to) try{ clearTimeout(to); }catch{} reject(e); };
				off = onFrame((res)=>{
					try {
						if (!res || typeof res !== 'object') return;
						const energy = (res.energy!=null ? res.energy : undefined);
						const forces = Array.isArray(res.forces)? res.forces : [];
						const isSimple = (res && typeof res.message === 'string' && res.message.toLowerCase() === 'simple_calculate');
						if (energy != null || isSimple) done({ energy, forces, userInteractionCount: res.userInteractionCount|0, simStep: res.simStep|0 });
					} catch(e){ /* ignore non-matching frames */ }
				});
					try {
						const msg = create(ClientActionSchema, {
							seq: nextSeq(),
							ack: clientAck||undefined,
							type: ClientAction_Type.SIMPLE_CALCULATE,
							userInteractionCount: lastCounters.userInteractionCount|0,
							simStep: lastCounters.simStep|0,
						});
						__notifyTestHook(msg);
						try { __wsLog('[WS][tx][SIMPLE_CALCULATE]', { seq: msg.seq|0, uic: msg.userInteractionCount|0, simStep: msg.simStep|0, ack: msg.ack? (msg.ack|0) : undefined, state: getState() }); } catch {}
						const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
					} catch(e){ fail(e); }
				// optional timeout
				to = setTimeout(()=> { __wsWarn('[WS][simple][timeout] after 5s'); fail(new Error('simple_calculate timeout')); }, 5000);
			});
		}

		// One-shot request for a single simulation step (md or relax). Starts, waits 1 frame, stops, resolves.
		async function requestSingleStep({ type, params }){
			if(!ws || ws.readyState!==1) { __wsWarn('[WS][single-step][not-connected]'); throw new Error('ws not connected'); }
			return new Promise((resolve, reject)=>{
				let off = null; let doneCalled=false;
				function done(v){ if(doneCalled) return; doneCalled=true; try{ off && off(); }catch{} try{ stopSimulation(); }catch{} resolve(v); }
				off = onFrame((res)=>{ try { if(res && typeof res==='object'){ __wsLog('[WS][rx][single-step]', { hasPos: !!res.positions, hasForces: !!res.forces, energy: res.energy }); done(res); } } catch{} });
				try { startSimulation({ type, params }); } catch(e){ try{ off&&off(); }catch{} reject(e); }
				setTimeout(()=>{ if(!doneCalled) { __wsWarn('[WS][single-step][timeout] after 5s'); try{ off&&off(); }catch{} reject(new Error('single_step timeout')); } }, 5000);
			});
		}
		    	const api = { connect, ensureConnected, sendBytes, close, onFrame, onResult, setCounters, nextSeq, setAck, ack, getState, requestSimpleCalculate, requestSingleStep, waitForEnergy, initSystem, userInteraction, beginDrag, endDrag, startSimulation, stopSimulation };
		            try {
				    if (typeof window !== 'undefined') {
		            window.__fairchem_ws__ = api;
		            // Test bridge: allow tests to inject decoded frames without protobuf by calling window.__ON_WS_RESULT__(obj)
		            if (typeof window.__ON_WS_RESULT__ !== 'function') {
		              window.__ON_WS_RESULT__ = (obj)=>{ try { for(const fn of listeners){ try { fn(obj, lastCounters); } catch{} } } catch{} };
		            }
		            // Toggle WS debug at runtime from tests or devtools
		            if (typeof window.__WS_DEBUG_ENABLE__ !== 'function') {
		              window.__WS_DEBUG_ENABLE__ = (on)=>{ try { window.__MLIPVIEW_DEBUG_WS = !!on; if (on && window.localStorage) window.localStorage.setItem('mlip.wsDebug','1'); else if (window.localStorage) window.localStorage.removeItem('mlip.wsDebug'); console.log('[WS] debug set to', !!on); } catch {} };
		            }
		          }
		          // Also expose test hook in Node/Jest environments without window
		          try {
		            if (typeof globalThis !== 'undefined') {
		              globalThis.__fairchem_ws__ = globalThis.__fairchem_ws__ || api;
		              if (typeof globalThis.__ON_WS_RESULT__ !== 'function') {
		                globalThis.__ON_WS_RESULT__ = (obj)=>{ try { for(const fn of listeners){ try { fn(obj, lastCounters); } catch{} } } catch{} };
		              }
		            }
		            if (typeof global !== 'undefined') {
		              global.__fairchem_ws__ = global.__fairchem_ws__ || api;
		              if (typeof global.__ON_WS_RESULT__ !== 'function') {
		                global.__ON_WS_RESULT__ = (obj)=>{ try { for(const fn of listeners){ try { fn(obj, lastCounters); } catch{} } } catch{} };
		              }
		            }
		          } catch {}
		          } catch {}
		            return api;
}

		export function getWS(){ if(!__singleton) __singleton = createFairchemWS(); return __singleton; }

export default { createFairchemWS };
