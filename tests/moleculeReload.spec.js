/**
 * Verifies that molecule selection triggers a full page reload via URL param,
 * and that loader honors ?mol=... on startup.
 */

describe('molecule reload selection flow', () => {
  test('buildReloadUrl encodes ?mol and preserves other params', async () => {
    const mod = await import('../public/ui/moleculeSelect.js');
    const href = 'http://localhost:4000/?debug=1#frag';
    const out = mod.buildReloadUrl(href, 'molecules/water.xyz');
    expect(out.startsWith('http://localhost:4000/?')).toBe(true);
    const url = new URL(out);
    expect(url.searchParams.get('debug')).toBe('1');
    expect(url.searchParams.get('mol')).toBe('molecules/water.xyz');
  });

  test('installMoleculeSelector sets location on change', async () => {
    const { installMoleculeSelector } = await import('../public/ui/moleculeSelect.js');
    const doc = { createElement: (tag)=> ({ tagName: tag, id:'', innerHTML:'', onchange:null, value:'', appendChild(){}, }) };
    const appended = [];
    const hudEl = { appendChild: (el)=> appended.push(el) };
    const loc = { href: 'http://localhost:4000/', assign: jest.fn() };
    const win = { location: loc };
    const { selectEl } = installMoleculeSelector({ hudEl, windowRef: win, documentRef: doc });
    expect(appended.length).toBe(1);
    selectEl.value = 'molecules/benzene.xyz';
    // simulate change
    selectEl.onchange();
    expect(loc.assign).toHaveBeenCalled();
    const dest = loc.assign.mock.calls[0][0];
    const url = new URL(dest);
    expect(url.searchParams.get('mol')).toBe('molecules/benzene.xyz');
  });

  test('loadDefault honors ?mol= from URL when valid', async () => {
    jest.resetModules();
    const { loadDefault, getRequestedMoleculeFromUrl } = await import('../public/util/moleculeLoader.js');
    // fake viewer API minimal surface with required state methods
    const viewerApi = { state: { markPositionsChanged(){}, markBondsChanged(){} }, recomputeBonds: jest.fn() };
    global.window = { location: { search: '?mol=molecules/water.xyz', origin:'http://localhost:4000' } };
    // mock fetch to return trivial xyz
    global.fetch = jest.fn(async ()=> ({ ok:true, text: async()=> '3\nwater\nH 0 0 0\nH 0 0 1\nO 1 0 0\n' }));
    expect(getRequestedMoleculeFromUrl()).toBe('molecules/water.xyz');
    const res = await loadDefault(viewerApi);
    expect(res.file).toBe('molecules/water.xyz');
    expect(Array.isArray(viewerApi.state.positions)).toBe(true);
    expect(viewerApi.recomputeBonds).toHaveBeenCalled();
  });
});
