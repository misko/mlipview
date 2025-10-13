/** @jest-environment jsdom */

// Tests for SMILES generation and XYZ upload flows using the System panel

describe('System panel SMILES + XYZ upload', () => {
  let origFetch;
  beforeAll(()=>{ origFetch = global.fetch; });
  afterAll(()=>{ global.fetch = origFetch; });

  beforeEach(() => {
    document.body.innerHTML = '<div id="app" style="position:relative"></div>';
    // Provide minimal window.location behavior for buildReloadUrlWithParam usage
    delete window.location; // jsdom allows reassignment
    window.location = { href: 'http://localhost:4000/?foo=bar', origin: 'http://localhost:4000', assign: function(h){ this.href = h; } };
  });

  async function buildPanel(){
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('app') });
  }

  test('SMILES generation success triggers navigation with smiles param', async () => {
    await buildPanel();
    const input = document.getElementById('smilesInput');
    const btn = document.getElementById('smilesGoBtn');
    input.value = 'CCO'; // ethanol
    input.dispatchEvent(new Event('input'));
    expect(btn.disabled).toBe(false);
    btn.click();
    expect(window.location.href.includes('smiles=CCO')).toBe(true);
  });

  test('SMILES generation failure warns and does not navigate', async () => {
    await buildPanel();
    const input = document.getElementById('smilesInput');
    const btn = document.getElementById('smilesGoBtn');
    input.value = 'bad smiles with spaces';
    input.dispatchEvent(new Event('input'));
    // Button still enabled because non-empty; click should validate fail
    btn.click();
    // No status element anymore; ensure navigation did not occur
    // href should not change to contain smiles=
    expect(window.location.href.includes('smiles=')).toBe(false);
  });

  test('XYZ upload success reloads with molxyz param', async () => {
    await buildPanel();
    const fileInput = document.getElementById('xyzFileInput');
    const uploadBtn = document.getElementById('uploadXyzBtn');
    expect(uploadBtn).toBeTruthy();
    // Spy assign to observe navigation
    window.location.assign = jest.fn(function(h){ this.href = h; });
    // Prepare a valid small XYZ (H2O)
    const xyz = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0';
  const file = new File([xyz], 'water.xyz', { type: 'text/plain' });
  // Some jsdom versions may not implement File.prototype.text reliably; force it
  Object.defineProperty(file, 'text', { configurable: true, value: () => Promise.resolve(xyz) });
  // Simulate file selection (jsdom: files is a readonly getter)
  Object.defineProperty(fileInput, 'files', { configurable: true, get: () => [file] });
    fileInput.dispatchEvent(new Event('change'));
    // Allow async handler to run (file.text() + navigation). Retry briefly.
    const start = Date.now();
    while (Date.now() - start < 50) {
      // If assign was used, assert and exit
      if (window.location.assign.mock.calls.length > 0) break;
      // If href already updated, we can also stop
      if (String(window.location.href).includes('molxyz=')) break;
      // wait a tick
      // eslint-disable-next-line no-await-in-loop
      await new Promise(r=>setTimeout(r, 5));
    }
    if (window.location.assign.mock.calls.length > 0) {
      const calledUrl = window.location.assign.mock.calls[0][0];
      expect(String(calledUrl)).toContain('molxyz=');
    } else {
      expect(String(window.location.href)).toContain('molxyz=');
    }
  });

  test('XYZ upload with blank comment line reloads with molxyz param', async () => {
    await buildPanel();
    const fileInput = document.getElementById('xyzFileInput');
    const uploadBtn = document.getElementById('uploadXyzBtn');
    expect(uploadBtn).toBeTruthy();
    window.location.assign = jest.fn(function(h){ this.href = h; });
    // XYZ with empty second line (valid XYZ format)
    const xyz = '3\n\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n';
    const file = new File([xyz], 'water_blank_comment.xyz', { type: 'text/plain' });
    Object.defineProperty(file, 'text', { configurable: true, value: () => Promise.resolve(xyz) });
    Object.defineProperty(fileInput, 'files', { configurable: true, get: () => [file] });
    fileInput.dispatchEvent(new Event('change'));
    // Wait for async to complete
    const start = Date.now();
    while (Date.now() - start < 50) {
      if (window.location.assign.mock.calls.length > 0) break;
      if (String(window.location.href).includes('molxyz=')) break;
      // eslint-disable-next-line no-await-in-loop
      await new Promise(r=>setTimeout(r, 5));
    }
    if (window.location.assign.mock.calls.length > 0) {
      const calledUrl = window.location.assign.mock.calls[0][0];
      expect(String(calledUrl)).toContain('molxyz=');
    } else {
      expect(String(window.location.href)).toContain('molxyz=');
    }
  });

  test('XYZ upload failure does not navigate', async () => {
    await buildPanel();
    const fileInput = document.getElementById('xyzFileInput');
    // Invalid XYZ (declares 2 but provides 1)
    const bad = '2\nbad\nH 0 0 0';
  const file = new File([bad], 'bad.xyz', { type: 'text/plain' });
  Object.defineProperty(fileInput, 'files', { configurable: true, get: () => [file] });
    const prevHref = window.location.href;
  fileInput.dispatchEvent(new Event('change'));
  await new Promise(r=>setTimeout(r,0));
    expect(window.location.href).toBe(prevHref);
  });
});
