/** @jest-environment jsdom */

import { buildDesktopPanel } from '../public/ui/desktopPanel.js';

function setViewport(width) {
  Object.defineProperty(window, 'innerWidth', { configurable: true, value: width });
  // jsdom doesn't implement matchMedia fully; provide a stub to satisfy our code
  if (!window.matchMedia) {
    window.matchMedia = (q) => ({
      matches: /max-width:\s*800px/.test(q) ? width <= 800 : false,
      media: q,
      addListener: () => {},
      removeListener: () => {},
      addEventListener: () => {},
      removeEventListener: () => {},
    });
  }
}

describe('responsive: desktop panel vs mobile top bar', () => {
  beforeEach(() => {
    document.body.innerHTML = '<div id="app"></div>';
    // default: desktop size
    setViewport(1200);
  });

  test('desktop: panel visible, mobile top bar hidden', () => {
    const { panelEl } = buildDesktopPanel({ attachTo: document.getElementById('app') });
    expect(panelEl).toBeTruthy();
    const topBar = document.getElementById('mobileTopBar');
    expect(topBar).toBeTruthy();
    // Desktop: top bar hidden
    expect(topBar.style.display).toBe('none');
    // Desktop: panel visible
    expect(panelEl.style.display).toBe('block');
  });

  test('mobile: panel hidden by default, top bar visible; tapping tabs toggles sheet', () => {
    setViewport(390);
    const { panelEl } = buildDesktopPanel({ attachTo: document.getElementById('app') });
    const topBar = document.getElementById('mobileTopBar');
    expect(topBar.style.display).toBe('flex');
    // Hidden until a tab is clicked
    expect(panelEl.style.display).toBe('none');

    // Tap Simulation tab
    const simBtn = document.getElementById('mobileTab-simulation');
    simBtn.click();
    expect(panelEl.getAttribute('data-mobile-open')).toBe('true');
    // Panel stays hidden on mobile; content renders inside the top bar
    expect(panelEl.style.display).toBe('none');
    expect(topBar.getAttribute('data-open')).toBe('true');
    const sheet = topBar.querySelector('.mobile-sheet');
    expect(sheet && sheet.style.display).toBe('block');
    // Verify that Simulation content is mounted into the sheet (has MD toggle)
    expect(sheet.querySelector('#toggleMD')).toBeTruthy();

    // Tap again to close
    simBtn.click();
    expect(panelEl.getAttribute('data-mobile-open')).toBeNull();
    expect(panelEl.style.display).toBe('none');
    expect(topBar.getAttribute('data-open')).toBeNull();
    const sheet2 = topBar.querySelector('.mobile-sheet');
    expect(sheet2 && sheet2.style.display).toBe('none');
  });
});
