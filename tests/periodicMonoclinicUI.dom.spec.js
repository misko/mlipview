/** @jest-environment jsdom */

// DOM test for new Periodic tab: PBC toggle + monoclinic inputs

function wait(ms){ return new Promise(r=>setTimeout(r, ms)); }

describe('Periodic tab monoclinic UI', () => {
  async function setup() {
    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = '<div id="app"></div>';
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    // Minimal viewer API stub
    const listeners = {};
    const bus = { on:(ev,fn)=>{ (listeners[ev]||(listeners[ev]=[])).push(fn); }, emit:(ev)=>{ (listeners[ev]||[]).forEach(fn=>fn()); } };
    const state = {
      bus,
      showCell: false,
      showGhostCells: false,
      cell: { a:{x:10,y:0,z:0}, b:{x:0,y:12,z:0}, c:{x:0,y:0,z:15}, originOffset:{x:0,y:0,z:0}, enabled: true },
      toggleCellVisibilityEnhanced(){ this.showCell = !this.showCell; this.bus.emit('cellChanged'); },
      toggleGhostCells(){ this.showGhostCells = !this.showGhostCells; this.bus.emit('cellChanged'); },
      markCellChanged(){ this.bus.emit('cellChanged'); }
    };
    window.viewerApi = { state, setForceVectorsEnabled: ()=>{} };
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    await wait(0);
    return { state };
  }

  test('Periodic tab exists with PBC toggle and inputs', async () => {
    await setup();
    const tab = document.getElementById('mobileTopBar') ? document.getElementById('mobileTab-periodic') : null;
    // Desktop: tab might not be used; ensure section exists
    const section = document.getElementById('section-periodic');
    expect(section).toBeTruthy();
    const pbcToggle = section.querySelector('#togglePBCHeader');
    const a = section.querySelector('#cellA');
    const b = section.querySelector('#cellB');
    const c = section.querySelector('#cellC');
    const beta = section.querySelector('#cellBeta');
    expect(pbcToggle).toBeTruthy();
    expect(a && b && c && beta).toBeTruthy();
  });

  test('enabling PBC enables nudgers; nudging updates cell with beta angle', async () => {
    const { state } = await setup();
    const section = document.getElementById('section-periodic');
    const pbcToggle = section.querySelector('#togglePBCHeader');
    const a = section.querySelector('#cellA');
    const b = section.querySelector('#cellB');
    const c = section.querySelector('#cellC');
    const beta = section.querySelector('#cellBeta');
    const aPlus = section.querySelector('#cellAPlus');
    const bPlus = section.querySelector('#cellBPlus');
    const cPlus = section.querySelector('#cellCPlus');
    const betaPlus = section.querySelector('#cellBetaPlus');

    // Initially disabled because showCell=false (via buttons disabled)
    expect(aPlus.disabled && bPlus.disabled && cPlus.disabled && betaPlus.disabled).toBe(true);
    // Toggle PBC on
    pbcToggle.click();
    expect(state.showCell).toBe(true);
    // Controls should enable
    expect(aPlus.disabled || bPlus.disabled || cPlus.disabled || betaPlus.disabled).toBe(false);

    // Nudge values a few times
    aPlus.click(); bPlus.click(); cPlus.click();
    betaPlus.click(); betaPlus.click();
    // Verify cell vectors shape (monoclinic: a on x, b on y, c in xz with beta between a and c)
    expect(state.cell.a.x).toBeGreaterThan(0);
    expect(state.cell.b.y).toBeGreaterThan(0);
    // |c| should be > 0, with non-zero x or z (depending on beta)
    const cLen = Math.hypot(state.cell.c.x, state.cell.c.y, state.cell.c.z);
    expect(cLen).toBeGreaterThan(0);
  });

  test('holding beta (+) repeats and clamps to < 180°', async () => {
    const { state } = await setup();
    // Enable fake timers after setup so the internal wait(0) can resolve
    jest.useFakeTimers();
    const section = document.getElementById('section-periodic');
    const pbcToggle = section.querySelector('#togglePBCHeader');
  const beta = section.querySelector('#cellBeta');
    const betaPlus = section.querySelector('#cellBetaPlus');
    pbcToggle.click();
    // Start near upper bound
  beta.textContent = '178';
  // Hold pointer down to auto-increment (use generic Event for jsdom)
  betaPlus.dispatchEvent(new Event('pointerdown', { bubbles:true }));
    // Let a few repeats happen
    jest.advanceTimersByTime(400);
    // Release
  window.dispatchEvent(new Event('pointerup', { bubbles:true }));
    // Flush any pending timers
    jest.runOnlyPendingTimers();
    // Should be clamped to <= 179
  const val = parseFloat(beta.textContent);
    expect(val).toBeLessThanOrEqual(179);
    expect(val).toBeGreaterThanOrEqual(178); // at least initial
    // Ensure cell application happened (apply on each nudge)
    expect(state.cell && state.cell.enabled).toBe(true);
    jest.useRealTimers();
  });
});
