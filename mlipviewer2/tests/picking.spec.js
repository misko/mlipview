import { createPickingService } from '../public/core/pickingService.js';

// Some ESM test environments may not expose jest global early; provide a light fallback.
if (typeof jest === 'undefined') {
  globalThis.jest = { fn: (impl = ()=>{}) => { const f = (...a)=>impl(...a); f.mock = { calls:[] }; const wrap = (...a)=>{ f.mock.calls.push(a); return impl(...a); }; return wrap; } };
}

describe('picking service (logic only)', () => {
  test('returns null when no hit', () => {
    const scene = { pointerX:0, pointerY:0, pick(){ return { hit:false }; }, onPointerObservable:{ add(){} } };
    const view = { resolveAtomPick:()=>null, resolveBondPick:()=>null };
    const selection = { clickAtom:jest.fn(), clickBond:jest.fn(), clear:jest.fn() };
    const picking = createPickingService(scene, view, selection);
    expect(picking.pickAtPointer()).toBeNull();
  });
});
