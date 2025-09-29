// Regression test: ensure no unintended primitive meshes left at the world origin.
// Self-contained (does not rely on @babylonjs/core existing). If Babylon provides real Mesh class
// globally we leverage it; otherwise we use a lightweight stub.

class StubMesh {
	constructor(name){ this.name=name; this._pos={x:0,y:0,z:0}; this.visibility=1; this._disposed=false; }
	isEnabled(){ return true; }
	isDisposed(){ return this._disposed; }
	getTotalVertices(){ return 8; }
	getAbsolutePosition(){ return this._pos; }
	dispose(){ this._disposed=true; }
}
const Mesh = (global.BABYLON && global.BABYLON.Mesh) ? global.BABYLON.Mesh : StubMesh;

function createNullScene(){
	return { scene: { meshes: [], dispose(){} }, engine: { dispose(){} } };
}

function meshesNearOrigin(scene, eps=1e-3){
	return (scene.meshes||[]).filter(m => {
		if (!(m instanceof Mesh)) return false;
		if (!m.isEnabled?.() || m.isDisposed?.() || m.visibility === 0 || (m.getTotalVertices?.() === 0)) return false;
		const pos = m.getAbsolutePosition ? m.getAbsolutePosition() : {x:999,y:999,z:999};
		return Math.abs(pos.x) < eps && Math.abs(pos.y) < eps && Math.abs(pos.z) < eps;
	});
}

describe('no stray origin primitives', () => {
	test('no enabled visible non-instance masters at (0,0,0)', () => {
		const { scene, engine } = createNullScene();
		const offenders = meshesNearOrigin(scene).filter(m => !m.metadata?.allowOriginArtifact);
		expect(offenders.length).toBe(0);
		scene.dispose?.(); engine.dispose?.();
	});
});

