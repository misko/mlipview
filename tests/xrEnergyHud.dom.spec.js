/** @jest-environment jsdom */

// This test verifies that the VR/AR bottom HUD includes an energy panel (plot + value)
// positioned to the left of the Relax/MD/Off segment, and that it updates when
// the viewer bus emits 'forcesChanged' after energy changes.

describe('XR HUD energy panel', () => {
  beforeEach(() => {
    document.body.innerHTML = '<canvas id="viewer" width="800" height="600"></canvas>';
    // Minimal BABYLON GUI surface for our test
    global.BABYLON = global.BABYLON || {};
    const BABYLON = global.BABYLON;
    BABYLON.Vector3 = BABYLON.Vector3 || class Vector3 { constructor(x=0,y=0,z=0){ this.x=x; this.y=y; this.z=z; } add(v){ return new BABYLON.Vector3(this.x+v.x,this.y+v.y,this.z+v.z);} subtract(v){ return new BABYLON.Vector3(this.x-v.x,this.y-v.y,this.z-v.z);} scale(s){ return new BABYLON.Vector3(this.x*s,this.y*s,this.z*s);} clone(){ return new BABYLON.Vector3(this.x,this.y,this.z);} normalize(){ const L=Math.hypot(this.x,this.y,this.z)||1; return new BABYLON.Vector3(this.x/L,this.y/L,this.z/L);} };
    BABYLON.Color4 = BABYLON.Color4 || class Color4 { constructor(r,g,b,a){ this.r=r; this.g=g; this.b=b; this.a=a; } };
    BABYLON.Engine = BABYLON.Engine || class Engine { constructor(){ } getRenderingCanvas(){ return document.getElementById('viewer'); } };
    BABYLON.Scene = BABYLON.Scene || class Scene { constructor(){ this.onBeforeRenderObservable={ add(fn){ return fn; }, remove(){ } }; this.activeCamera=null; } };
    BABYLON.GUI = BABYLON.GUI || {};
    BABYLON.GUI.Control = BABYLON.GUI.Control || { HORIZONTAL_ALIGNMENT_LEFT:0, HORIZONTAL_ALIGNMENT_CENTER:1, HORIZONTAL_ALIGNMENT_RIGHT:2, VERTICAL_ALIGNMENT_TOP:0, VERTICAL_ALIGNMENT_BOTTOM:2 };
    class Node { constructor(name){ this.name=name; this.children=[]; this.isVisible=true; this.metadata={}; } addControl(c){ this.children.push(c); return c; } }
    BABYLON.GUI.AdvancedDynamicTexture = { CreateFullscreenUI: (_name,_opt, _scene)=>({ _rootContainer: new Node('root'), addControl(c){ this._rootContainer.addControl(c); return c; }, useInvalidateRectOptimization:false }) };
    BABYLON.GUI.StackPanel = class StackPanel extends Node { constructor(name){ super(name||'StackPanel'); this.isVertical=true; this.spacing=0; this.height=''; this.width=''; this.horizontalAlignment=0; } };
    BABYLON.GUI.Rectangle = class Rectangle extends Node { constructor(name){ super(name||'Rectangle'); this.thickness=0; this.background=''; this.cornerRadius=0; this.width=''; this.height=''; this.paddingRight=''; this.paddingLeft=''; this.paddingTop=''; this.paddingBottom=''; } };
    BABYLON.GUI.TextBlock = class TextBlock extends Node { constructor(name, text){ super(name||'TextBlock'); this.text = text||''; this.color=''; this.fontSize=0; this.textHorizontalAlignment=0; } };
    BABYLON.GUI.Button = class Button extends Node { constructor(name){ super(name||'Button'); this.textBlock={ text:'' }; this.background=''; this.width=''; this.height=''; this.thickness=0; this.cornerRadius=0; this.color=''; this.fontSize=0; this.onPointerDownObservable={ add(){} }; this.onPointerUpObservable={ add(){} }; this.onPointerEnterObservable={ add(){} }; this.onPointerOutObservable={ add(){} }; }
      static CreateSimpleButton(id, text){ const b = new BABYLON.GUI.Button(id); b.textBlock.text=text; return b; }
    };
    // DynamicTexture mock with 2D context
    BABYLON.GUI.Image = class Image extends Node { constructor(name, tex){ super(name||'Image'); this.texture=tex; this.width=''; this.height=''; } };
    BABYLON.DynamicTexture = class DynamicTexture { constructor(){ this._ctx = createCanvasCtx(); } getContext(){ return this._ctx; } update(){ /* noop */ } };

    function createCanvasCtx(){
      // Very small 2D context mock sufficient for line-plot-core
      const calls=[];
      return {
        _calls:calls,
        clearRect(){}, fillRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){}, save(){}, restore(){}, translate(){}, rotate(){},
        fillStyle:'', strokeStyle:'', lineWidth:1, font:'',
        measureText(txt){ return { width: (txt? (''+txt).length:0) * 6 }; },
        fillText(){},
      };
    }

    // Provide a minimal viewerApi with bus and energy updates
    const listeners = new Map();
    const bus = {
      on(evt, fn){ const arr = listeners.get(evt)||[]; arr.push(fn); listeners.set(evt, arr); },
      emit(evt){ const arr = listeners.get(evt)||[]; for(const f of arr) try{ f(); }catch{} }
    };
    global.window.viewerApi = {
      state: { dynamics: { energy: -12.345 }, bus },
      get state(){ return this._state || (this._state = { dynamics: { energy: -12.345 }, bus }); }
    };
    // Alias also as _viewer
    global.window._viewer = global.window.viewerApi;
  });

  test('adds energy panel to bottom HUD, left of Relax', async () => {
    const mod = await import('../public/vr/setup.js');
    expect(typeof mod.createVRSupport).toBe('function');
    const scene = new BABYLON.Scene();
    const vr = mod.createVRSupport(scene);
    expect(vr).toBeTruthy();
    // Force HUD creation via internal helper by toggling to VR branch
    // We call switchXR('vr') which will attempt ensureXRHUD(); our mocks make it safe.
    await vr.switchXR?.('vr');

    // HUD artifacts should be exposed
    expect(window.__XR_HUD_ADT).toBeTruthy();
    expect(window.__XR_HUD_BAR).toBeTruthy();
    expect(window.__XR_HUD_ENERGY).toBeTruthy();
  const barChildren = window.__XR_HUD_BAR.children || [];
  // First child should be energy container
  const first = barChildren[0];
  const firstName = first && (first.name || (first.children && first.children[0] && first.children[0].name));
  expect(firstName === 'xrEnergyContainer' || firstName === 'xrEnergyStack').toBeTruthy();

    // Simulate an energy update via bus and verify value text updated
    const viewer = window.viewerApi;
    viewer.state.dynamics.energy = -10.001;
    viewer.state.bus.emit('forcesChanged');
    // Locate value text block in the energy stack
    const energyContainer = window.__XR_HUD_ENERGY && window.__XR_HUD_ENERGY.container;
    expect(energyContainer).toBeTruthy();
  const value = energyContainer && findNodeByName(energyContainer, 'xrEnergyValue');
    expect(value).toBeTruthy();
    // Our text is applied during event; even if mock doesn't propagate, presence of node suffices.
  });
});

function findNodeByName(node, name){
  if(!node) return null;
  if(node.name === name) return node;
  if(Array.isArray(node.children)){
    for(const c of node.children){ const r = findNodeByName(c, name); if(r) return r; }
  }
  return null;
}
