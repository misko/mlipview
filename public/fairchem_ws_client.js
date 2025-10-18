// Thin WebSocket client for fairchem_local_server2 /ws protobuf stream
// Non-breaking addition: not yet integrated into index.js loops; provides a
// small API for future migration.

import { __count } from './util/funcCount.js';
// Protobuf stubs will be loaded on-demand in the browser by injecting scripts.
let __pb = (globalThis && globalThis.proto && globalThis.proto.fairchem && globalThis.proto.fairchem.session) || {};
let __hasPB = !!(__pb && __pb.ClientAction && __pb.ServerResult);
async function __ensurePB(){
	if (__hasPB) return true;
	if (typeof window === 'undefined' || typeof document === 'undefined') {
		throw new Error('[WS][protobuf-missing] session_pb.js not loaded in non-browser environment');
	}
	// Provide a minimal require shim for generated stubs
	if (!('require' in window)) {
		window.require = function(name){
			if (name === 'google-protobuf') {
				if (window.jspb) return window.jspb;
				// Fallback: some builds attach exports under goog.jspb
				if (window.goog && window.goog.jspb) return window.goog.jspb;
				throw new Error('[WS] google-protobuf not loaded');
			}
			throw new Error('Unsupported module '+name);
		};
	}
	// 1) google-protobuf runtime
	if (!window.goog) {
			await new Promise((resolve, reject)=>{
				const s = document.createElement('script');
				// Load the browser-ready Closure runtime (classic script, defines window.goog/jspb)
				s.src = '/vendor/google-protobuf/google-protobuf.js';
			s.async = false;
			s.onload = ()=> resolve(true);
			s.onerror = ()=> reject(new Error('Failed to load google-protobuf runtime'));
			document.head.appendChild(s);
		});
	}
	// 2) generated session_pb.js (commonjs style expects require)
	if (!(globalThis && globalThis.proto && globalThis.proto.fairchem && globalThis.proto.fairchem.session)) {
		await new Promise((resolve, reject)=>{
			const s = document.createElement('script');
			s.src = '/proto/fairchem_local_server2/session_pb.js';
			s.async = false;
			s.onload = ()=> resolve(true);
			s.onerror = ()=> reject(new Error('Failed to load session_pb.js'));
			document.head.appendChild(s);
		});
	}
	__pb = (globalThis && globalThis.proto && globalThis.proto.fairchem && globalThis.proto.fairchem.session) || {};
	__hasPB = !!(__pb && __pb.ClientAction && __pb.ServerResult);
	if (!__hasPB) throw new Error('[WS][protobuf-missing] Unable to initialize protobuf stubs');
	return true;
}

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
						const r = __pb.ServerResult.deserializeBinary(bytes);
						const out = {
							seq: r.getSeq && r.getSeq(),
							client_seq: r.getClientSeq && r.getClientSeq(),
							userInteractionCount: (typeof r.getUserInteractionCount==='function')? r.getUserInteractionCount(): undefined,
							simStep: (typeof r.getSimStep==='function')? r.getSimStep(): undefined,
							message: r.getMessage && r.getMessage(),
						};
						if (typeof r.hasEnergy==='function' && r.hasEnergy()) out.energy = r.getEnergy();
						if (r.getPositionsList) out.positions = r.getPositionsList().map(v=> v.getVList());
						if (r.getForcesList) out.forces = r.getForcesList().map(v=> v.getVList());
						if (r.getVelocitiesList) out.velocities = r.getVelocitiesList().map(v=> v.getVList());
						if (r.getCell && r.getCell()){
							const m = r.getCell().getMList();
							if (Array.isArray(m) && m.length===9) out.cell = [ [m[0],m[1],m[2]],[m[3],m[4],m[5]],[m[6],m[7],m[8]] ];
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
			const msg = new __pb.ClientAction();
			msg.setSeq(nextSeq());
			msg.setType(__pb.ClientAction.Type.INIT_SYSTEM);
			if (Array.isArray(atomic_numbers)) {
				for (const z of atomic_numbers) msg.addAtomicNumbers(z|0);
			}
			if (Array.isArray(positions)) {
				for (const p of positions) {
					const v = new __pb.Vec3(); v.setVList([+p[0], +p[1], +p[2]]); msg.addPositions(v);
				}
			}
			if (Array.isArray(velocities)) {
				for (const p of velocities) {
					const v = new __pb.Vec3(); v.setVList([+p[0], +p[1], +p[2]]); msg.addVelocities(v);
				}
			}
			if (cell && Array.isArray(cell) && cell.length===3) {
				const m = new __pb.Mat3();
				const flat = [
					+cell[0][0], +cell[0][1], +cell[0][2],
					+cell[1][0], +cell[1][1], +cell[1][2],
					+cell[2][0], +cell[2][1], +cell[2][2],
				];
				m.setMList(flat);
				if (typeof msg.setCell === 'function') msg.setCell(m);
			}
			const buf = msg.serializeBinary();
			sendBytes(buf);
		}

				// UPDATE_POSITIONS: keep server-side state synchronized for one-shot calculations
		function updatePositions(positions){
			const msg = new __pb.ClientAction();
			msg.setSeq(nextSeq());
			msg.setType(__pb.ClientAction.Type.UPDATE_POSITIONS);
			if (Array.isArray(positions)) {
				for (const p of positions) {
					const v = new __pb.Vec3(); v.setVList([+p[0], +p[1], +p[2]]); msg.addPositions(v);
				}
			}
			const buf = msg.serializeBinary();
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
				const msg = new __pb.ClientAction();
				msg.setSeq(nextSeq());
				msg.setType(__pb.ClientAction.Type.START_SIMULATION);
				const t = String(type||'md').toLowerCase()==='md' ? __pb.ClientAction.SimType.MD : __pb.ClientAction.SimType.RELAX;
				if (typeof msg.setSimulationType==='function') msg.setSimulationType(t);
				if (params){
					const sp = new __pb.SimulationParams();
					if (params.calculator) sp.setCalculator(String(params.calculator));
					if (typeof params.temperature==='number') sp.setTemperature(+params.temperature);
					if (typeof params.timestep_fs==='number') sp.setTimestepFs(+params.timestep_fs);
					if (typeof params.friction==='number') sp.setFriction(+params.friction);
					if (typeof params.fmax==='number') sp.setFmax(+params.fmax);
					if (typeof params.max_step==='number') sp.setMaxStep(+params.max_step);
					if (params.optimizer) sp.setOptimizer(String(params.optimizer));
					if (typeof msg.setSimulationParams==='function') msg.setSimulationParams(sp);
				}
				if (typeof msg.setUserInteractionCount==='function') msg.setUserInteractionCount(lastCounters.userInteractionCount|0);
				if (typeof msg.setSimStep==='function') msg.setSimStep(lastCounters.simStep|0);
				const buf = msg.serializeBinary(); sendBytes(buf);
			}

			function stopSimulation(){
				const msg = new __pb.ClientAction();
				msg.setSeq(nextSeq());
				msg.setType(__pb.ClientAction.Type.STOP_SIMULATION);
				const buf = msg.serializeBinary(); sendBytes(buf);
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
						const msg = new __pb.ClientAction(); msg.setSeq(nextSeq()); if (clientAck) msg.setAck(clientAck); msg.setType(__pb.ClientAction.Type.SIMPLE_CALCULATE); if (typeof msg.setUserInteractionCount === 'function') { msg.setUserInteractionCount(lastCounters.userInteractionCount|0); } if (typeof msg.setSimStep === 'function') { msg.setSimStep(lastCounters.simStep|0); } const buf = msg.serializeBinary(); sendBytes(buf);
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
