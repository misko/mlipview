/** @jest-environment jsdom */

// Verify the Selection section renders and updates for atom and bond selections.

jest.mock('../public/render/scene.js', () => ({ createScene: async () => ({ engine:{ runRenderLoop:(fn)=>{} }, scene:{ meshes:[], render:()=>{}, onPointerObservable:{ _l:[], add(fn){this._l.push(fn);} } }, camera:{ attachControl:()=>{} } }) }));

// Minimal fetch stub to satisfy viewer init API calls
global.fetch = async () => ({ ok:true, status:200, json: async ()=> ({ results:{} }) });

function wait(ms){ return new Promise(r=>setTimeout(r, ms)); }

describe('Selection panel UI', () => {
  async function setup(){
    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements:['N','O'],
      positions:[ {x:0,y:0,z:0}, {x:1,y:0,z:0} ],
      bonds:[ { i:0, j:1 } ]
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    // tiny delay to allow initial syncs
    await wait(0);
    return viewer;
  }

  test('renders Selection section and initial empty state', async () => {
    await setup();
    const sec = document.getElementById('section-selection');
    expect(sec).toBeTruthy();
    const a = document.getElementById('selSphereA');
    const b = document.getElementById('selSphereB');
    const elName = document.getElementById('selElementName');
    const pos = document.getElementById('selPosition');
    const w = document.getElementById('selAtomicWeight');
    const vdw = document.getElementById('selVdw');
    expect(a).toBeTruthy(); expect(b).toBeTruthy(); expect(elName).toBeTruthy();
    expect(a.style.display).toBe('none');
    expect(b.style.display).toBe('none');
    expect(elName.textContent).toBe('—');
    expect(pos.textContent).toBe('(-,-,-)');
    expect(w.textContent).toBe('—');
    expect(vdw.textContent).toBe('—');
  });

  test('atom selection shows one sphere, highlights element and position', async () => {
    const v = await setup();
    v.selection.clickAtom(0); // N at (0,0,0)
    const a = document.getElementById('selSphereA');
    const b = document.getElementById('selSphereB');
    const elName = document.getElementById('selElementName');
    const pos = document.getElementById('selPosition');
    const w = document.getElementById('selAtomicWeight');
    const vdw = document.getElementById('selVdw');
    expect(a.style.display).toBe('block');
    expect(b.style.display).toBe('none');
  expect((elName.textContent||'').toLowerCase()).toContain('nitrogen');
  expect(pos.textContent).toContain('(0.0');
    expect(w.textContent).toBe('14.007');
    expect(vdw.textContent).toBe('1.55');
    // periodic table highlight
    const nCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="N"]');
    const oCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="O"]');
    expect(nCell && nCell.classList.contains('highlight')).toBe(true);
    expect(oCell && oCell.classList.contains('highlight')).toBe(false);
  });

  test('bond selection shows two spheres, highlights both and shows bond length', async () => {
    const v = await setup();
    v.selection.clickBond({ i:0, j:1, index:0, key:'0-1' });
    const a = document.getElementById('selSphereA');
    const b = document.getElementById('selSphereB');
    const elName = document.getElementById('selElementName');
    const pos = document.getElementById('selPosition');
    const w = document.getElementById('selAtomicWeight');
    const vdw = document.getElementById('selVdw');
    const bondLen = document.getElementById('bondLength');
    expect(a.style.display).toBe('block');
    expect(b.style.display).toBe('block');
    const name = (elName.textContent||'').toLowerCase();
    expect(name.includes('nitrogen') && name.includes('oxygen')).toBe(true);
    expect(pos.textContent).toContain(') – (');
    expect(w.textContent).toBe('14.007 – 15.999');
    expect(vdw.textContent).toBe('1.55 – 1.52');
    expect(bondLen.textContent).toContain('1.000');
    expect(bondLen.textContent).toContain('Å');
    const nCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="N"]');
    const oCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="O"]');
    expect(nCell && nCell.classList.contains('highlight')).toBe(true);
    expect(oCell && oCell.classList.contains('highlight')).toBe(true);
  });

  test('clearing selection hides spheres and resets highlights', async () => {
    const v = await setup();
    v.selection.clickAtom(0);
    v.selection.clear();
    const a = document.getElementById('selSphereA');
    const b = document.getElementById('selSphereB');
    const elName = document.getElementById('selElementName');
    const pos = document.getElementById('selPosition');
    const w = document.getElementById('selAtomicWeight');
    const vdw = document.getElementById('selVdw');
    const nCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="N"]');
    const oCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="O"]');
    expect(a.style.display).toBe('none');
    expect(b.style.display).toBe('none');
    expect(elName.textContent).toBe('—');
    expect(pos.textContent).toBe('(-,-,-)');
    expect(w.textContent).toBe('—');
    expect(vdw.textContent).toBe('—');
    expect(nCell && nCell.classList.contains('highlight')).toBe(false);
    expect(oCell && oCell.classList.contains('highlight')).toBe(false);
  });
});
