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

function resolveWsBase(){
		if (typeof window === 'undefined') return 'ws://127.0.0.1:8000';
	const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
	const host = window.__MLIPVIEW_SERVER ? new URL(window.__MLIPVIEW_SERVER).host : window.location.host;
	return `${proto}//${host}`;
}

let __singleton = null;

export function createFairchemWS(){
	__count('ws#createFairchemWS');
	let ws = null;
	let seq = 0;
	let clientAck = 0;
	let lastCounters = { userInteractionCount: 0, simStep: 0 };
	const listeners = new Set();
	function onFrame(fn){ listeners.add(fn); return ()=>listeners.delete(fn); }

	async function connect(){
		await __ensurePB();
		const base = resolveWsBase();
		const url = base.replace(/\/$/, '') + '/ws';
		ws = new WebSocket(url);
			return new Promise((resolve, reject)=>{
			ws.binaryType = 'arraybuffer';
			ws.onopen = ()=> resolve(true);
			ws.onerror = (e)=> reject(e);
			ws.onclose = ()=>{};
			            ws.onmessage = (ev)=>{
					// Decode protobuf ServerResult frames only (no JSON fallback)
								try {
									if (!(ev.data instanceof ArrayBuffer)) {
										throw new Error('[WS] Expected binary protobuf frame, got non-binary message');
									}
									const bytes = new Uint8Array(ev.data);
									const r = fromBinary(ServerResultSchema, bytes);
									const out = {
										seq: r.seq,
										client_seq: r.clientSeq,
										userInteractionCount: r.userInteractionCount,
										simStep: r.simStep,
										message: r.message,
									};
									if (r.energy != null) out.energy = r.energy;
									if (Array.isArray(r.positions)) out.positions = r.positions.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									if (Array.isArray(r.forces)) out.forces = r.forces.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									if (Array.isArray(r.velocities)) out.velocities = r.velocities.map(v=> v && Array.isArray(v.v) ? v.v : undefined).filter(Boolean);
									if (r.cell && Array.isArray(r.cell.m)){
										const m = r.cell.m;
										if (m.length===9) out.cell = [ [m[0],m[1],m[2]],[m[3],m[4],m[5]],[m[6],m[7],m[8]] ];
									}
									for(const fn of listeners){ try { fn(out, lastCounters); } catch{} }
								} catch (e) {
						// Surface decode errors to console but do not crash the socket
						try { console.error('[WS][decode-error]', e?.message||e); } catch {}
					}
				};
		});
	}
					function sendBytes(buf){ if(ws && ws.readyState===1) ws.send(buf); }
	function close(){ try { ws && ws.close(); } catch{} }
	function setCounters({ userInteractionCount, simStep }){
		if(Number.isFinite(userInteractionCount)) lastCounters.userInteractionCount = userInteractionCount|0;
		if(Number.isFinite(simStep)) lastCounters.simStep = simStep|0;
	}
	function nextSeq(){ seq = (seq|0) + 1; return seq; }
	function setAck(s){ clientAck = Math.max(clientAck, s|0); }
	function getState(){ return { seq, clientAck, ...lastCounters, connected: !!ws && ws.readyState===1 }; }
		async function ensureConnected(){ if(ws && ws.readyState===1) return true; await connect(); return true; }

				// INIT_SYSTEM: send atoms and initial positions (and optional velocities/cell)
		function initSystem({ atomic_numbers, positions, velocities, cell }){
			const msg = create(ClientActionSchema, {
				seq: nextSeq(),
				type: ClientAction_Type.INIT_SYSTEM,
				atomicNumbers: Array.isArray(atomic_numbers) ? atomic_numbers.map(z=> z|0) : [],
				positions: Array.isArray(positions) ? positions.map(p=> create(Vec3Schema, { v: [+p[0],+p[1],+p[2]] })) : [],
				velocities: Array.isArray(velocities) ? velocities.map(p=> create(Vec3Schema, { v: [+p[0],+p[1],+p[2]] })) : [],
				cell: (cell && Array.isArray(cell) && cell.length===3) ? create(Mat3Schema, { m: [
					+cell[0][0], +cell[0][1], +cell[0][2],
					+cell[1][0], +cell[1][1], +cell[1][2],
					+cell[2][0], +cell[2][1], +cell[2][2],
				] }) : undefined,
			});
			const buf = toBinary(ClientActionSchema, msg);
			sendBytes(buf);
		}

				// UPDATE_POSITIONS: keep server-side state synchronized for one-shot calculations
		function updatePositions(positions){
			const msg = create(ClientActionSchema, {
				seq: nextSeq(),
				type: ClientAction_Type.UPDATE_POSITIONS,
				positions: Array.isArray(positions) ? positions.map(p=> create(Vec3Schema, { v: [+p[0],+p[1],+p[2]] })) : [],
			});
			const buf = toBinary(ClientActionSchema, msg);
			sendBytes(buf);
		}

				// Subscribe to decoded ServerResult frames (convenience over raw onFrame)
			function onResult(fn){
					// With binary-only decode in ws.onmessage, listeners already receive decoded objects
					return onFrame((decoded)=>{
						try {
							if (!decoded || typeof decoded !== 'object') return;
							try { if (typeof window !== 'undefined') { window.__ON_WS_RESULT__ && window.__ON_WS_RESULT__(decoded); } } catch {}
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
				const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
			}

			function stopSimulation(){
				const msg = create(ClientActionSchema, { seq: nextSeq(), type: ClientAction_Type.STOP_SIMULATION });
				const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
			}

				function ack(seqNum){
					try {
						const msg = new __pb.ClientAction();
						msg.setSeq(nextSeq());
						msg.setType(__pb.ClientAction.Type.PING);
						if (typeof msg.setAck==='function') msg.setAck(seqNum|0);
						const buf = msg.serializeBinary(); sendBytes(buf);
					} catch {}
				}
		// High-level helper: request one-shot simple calculation via WS.
			// Encodes a ClientAction SIMPLE_CALCULATE referencing current session state.
			// Decodes first matching ServerResult with energy/forces and resolves.
			async function requestSimpleCalculate(){
			if(!ws || ws.readyState!==1) throw new Error('ws not connected');
			return new Promise((resolve, reject)=>{
				let off = null;
				const done = (v)=>{ try{ off && off(); }catch{} resolve(v); };
				const fail = (e)=>{ try{ off && off(); }catch{} reject(e); };
				off = onFrame((res)=>{
					try {
						if (!res || typeof res !== 'object') return;
						const energy = (res.energy!=null ? res.energy : undefined);
						const forces = Array.isArray(res.forces)? res.forces : [];
						if (energy != null || (forces && forces.length)) done({ energy, forces });
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
						const buf = toBinary(ClientActionSchema, msg); sendBytes(buf);
					} catch(e){ fail(e); }
				// optional timeout
				setTimeout(()=> fail(new Error('simple_calculate timeout')), 5000);
			});
		}
	    		const api = { connect, ensureConnected, sendBytes, close, onFrame, onResult, setCounters, nextSeq, setAck, ack, getState, requestSimpleCalculate, initSystem, updatePositions, startSimulation, stopSimulation };
		            try { if (typeof window !== 'undefined') { window.__fairchem_ws__ = api; } } catch {}
		            return api;
}

		export function getWS(){ if(!__singleton) __singleton = createFairchemWS(); return __singleton; }

export default { createFairchemWS };
