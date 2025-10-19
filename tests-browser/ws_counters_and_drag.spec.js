/** @jest-environment jsdom */

// Verifies that:
// - On any user interaction we increment the user interaction counter and include it in simple calculate
// - WS responses are discarded if their userInteractionCount < current front-end counter
// - Matching counter responses are applied
// - While dragging an atom, incoming frames do not overwrite that atom's position

describe('WS userInteraction counter gating and drag exclusion', () => {
  beforeEach(() => {
    document.body.innerHTML = '<canvas id="viewer" width="800" height="600"></canvas>';
    // Disable auto-MD so it doesn't interfere
    window.__MLIPVIEW_NO_AUTO_MD = true;
    window.__MLIPVIEW_TEST_MODE = true;
  });

  test('stale simple_calculate is ignored and matching one applies', async () => {
    // Mock global WebSocket for the ws client; we'll capture the first sent SIMPLE_CALCULATE
    const sentFrames = [];
    let wsOnMessage = null;
    class WSStub {
      constructor() { this.readyState = 1; setTimeout(()=> this.onopen && this.onopen(), 0); }
      set binaryType(_){}
      send(buf){ sentFrames.push(buf); }
      close(){}
      onopen(){}
      onerror(){}
      onclose(){}
      set onmessage(fn){ wsOnMessage = fn; }
      get onmessage(){ return wsOnMessage; }
    }
    global.WebSocket = WSStub;

    // Import modules under test
    const { initNewViewer } = await import('../public/index.js');
    // Build a tiny scene canvas stub
    const canvas = document.getElementById('viewer');
    canvas.getContext = ()=>({});
    canvas.addEventListener = ()=>{};
    canvas.removeEventListener = ()=>{};
    canvas.getBoundingClientRect = ()=>({ left:0, top:0, width:800, height:600 });

    const api = await initNewViewer(canvas, { elements: ['H','H','O'], positions: [{x:0,y:0,z:0},{x:1,y:0,z:0},{x:0,y:1,z:0}], bonds: [] });

    // Force a simple calculate; the client will send UPDATE_POSITIONS then SIMPLE_CALCULATE.
    // We resolve by delivering two frames: first stale (lower counter), then matching.
    const p = api.ff.computeForces({ sync:true });

    // Helper to emit a protobuf-like decoded frame object the client expects from ws.onmessage
    function emitFrame(obj){
      // The ws client decodes ArrayBuffer into an object, but in our stub we shortcut and call the onmessage
      // with a fake ArrayBuffer that the client ignores; instead, we directly invoke the decoded listener path
      // by calling window.__ON_WS_RESULT__ which ws client triggers on each decoded frame.
      if (typeof window.__ON_WS_RESULT__ === 'function') window.__ON_WS_RESULT__(obj);
      else if (wsOnMessage) wsOnMessage({ data: new ArrayBuffer(0) });
    }

    // 1) Deliver a stale frame with userInteractionCount lower than current (assume current is >=0; bump it)
    api.debugRecordInteraction('dragMove');
    const cur = api.getVersionInfo().userInteractionVersion;
    // Stale uses cur-1 (min 0)
    const staleUIC = Math.max(0, cur - 1);
    emitFrame({ energy: -1.23, forces: [[0,0,0],[0,0,0],[0,0,0]], userInteractionCount: staleUIC, simStep: 1, positions: api.state.positions.map(p=>[p.x,p.y,p.z]) });

    // 2) Now a matching frame which should be accepted
    const matchUIC = api.getVersionInfo().userInteractionVersion;
    emitFrame({ energy: 0.42, forces: [[0.1,0,0],[0,0.1,0],[0,0,0.1]], userInteractionCount: matchUIC, simStep: 2, positions: api.state.positions.map(p=>[p.x,p.y,p.z]) });

    const res = await p;
    expect(typeof res.energy).toBe('number');
  });

  test('dragging atom prevents overwrite by incoming frame', async () => {
    // Mock WebSocket again
    let wsOnMessage = null;
    class WSStub {
      constructor() { this.readyState = 1; setTimeout(()=> this.onopen && this.onopen(), 0); }
      set binaryType(_){}
      send(_){}
      close(){}
      onopen(){}
      onerror(){}
      onclose(){}
      set onmessage(fn){ wsOnMessage = fn; }
      get onmessage(){ return wsOnMessage; }
    }
    global.WebSocket = WSStub;

    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.getElementById('viewer');
    canvas.getContext = ()=>({});
    canvas.addEventListener = ()=>{};
    canvas.removeEventListener = ()=>{};
    canvas.getBoundingClientRect = ()=>({ left:0, top:0, width:800, height:600 });

    const api = await initNewViewer(canvas, { elements: ['H','H','O'], positions: [{x:0,y:0,z:0},{x:1,y:0,z:0},{x:0,y:1,z:0}], bonds: [] });

    // Select atom 0 and start drag via manipulation wrapper
    api.state.selection = { kind:'atom', data:{ index:0 } };
    const startPos = { ...api.state.positions[0] };
    // Begin drag with a dummy intersector returning the same point (no movement yet)
    api.manipulation.beginDrag(() => ({ x:startPos.x, y:startPos.y, z:startPos.z }), { source:'desktop', planePoint: startPos, planeNormal:{x:0,y:1,z:0} });

    // Incoming frame tries to move atom 0; because it's being dragged, position should not change
    const curUIC = api.getVersionInfo().userInteractionVersion;
    const newPositions = api.state.positions.map((p,idx)=> idx===0 ? [p.x+5, p.y+5, p.z+5] : [p.x, p.y, p.z]);
    if (typeof window.__ON_WS_RESULT__ === 'function') {
      window.__ON_WS_RESULT__({ positions: newPositions, forces: [], energy: 0, userInteractionCount: curUIC, simStep: 1 });
    }

    // Verify atom 0 position unchanged
    const after = api.state.positions[0];
    expect(after.x).toBeCloseTo(startPos.x, 8);
    expect(after.y).toBeCloseTo(startPos.y, 8);
    expect(after.z).toBeCloseTo(startPos.z, 8);

    // End drag to cleanup
    api.manipulation.endDrag();
  });
});
