/** @jest-environment jsdom */

// Ensures that when an XYZ with cell parameters is loaded, the periodic toggle reflects 'On'.

import { parseXYZ } from '../public/util/xyzLoader.js';
import { applyParsedToViewer } from '../public/util/moleculeLoader.js';
import { buildDesktopPanel } from '../public/ui/desktopPanel.js';

function makeViewer() {
  const listeners = {};
  const bus = {
    on: (ev, fn) => {
      (listeners[ev] || (listeners[ev] = [])).push(fn);
    },
    emit: (ev) => {
      (listeners[ev] || []).forEach((fn) => fn());
    },
  };
  const state = {
    bus,
    elements: [],
    positions: [],
    bonds: [],
    cell: {
      a: { x: 1, y: 0, z: 0 },
      b: { x: 0, y: 1, z: 0 },
      c: { x: 0, y: 0, z: 1 },
      enabled: false,
      originOffset: { x: 0, y: 0, z: 0 },
    },
    showCell: false,
    markCellChanged() {
      bus.emit('cellChanged');
    },
    markPositionsChanged() {
      bus.emit('positionsChanged');
    },
    markBondsChanged() {
      bus.emit('bondsChanged');
    },
    dynamics: {},
  };
  const viewerApi = { state, recomputeBonds: () => {} };
  return viewerApi;
}

// Minimal DOM hook that desktopPanel expects
beforeEach(() => {
  document.body.innerHTML = '<div id="app"></div>';
  // Expose a fake viewerApi to window so desktopPanel can read it
  window.viewerApi = makeViewer();
});

function aceticHeader() {
  return `32\ncell: a=13.31, b=4.09, c=5.77, alpha=90.00, beta=107.0, gamma=90.00, spacegroup=P21/c, temp=400K\nC 0 0 0\nH 0 0 1\nH 0 1 0\nH 1 0 0\nC 2 0 0\nH 2 0 1\nH 2 1 0\nH 3 0 0\nC 4 0 0\nO 4 0 1\nO 4 1 0\nC 5 0 0\nH 5 0 1\nH 5 1 0\nH 6 0 0\nH 6 0 1\nC 7 0 0\nO 7 0 1\nO 7 1 0\nC 8 0 0\nH 8 0 1\nH 8 1 0\nH 9 0 0\nH 9 0 1\nC 10 0 0\nO 10 0 1\nO 10 1 0\nC 11 0 0\nH 11 0 1\nH 11 1 0\nH 12 0 0\nH 12 0 1\n`;
}

describe('Periodic UI toggle reflects PBC On after XYZ load with cell', () => {
  test('desktop panel toggle shows On', () => {
    // Build panel first so it subscribes to cellChanged
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    // Parse and apply XYZ
    const parsed = parseXYZ(aceticHeader());
    applyParsedToViewer(window.viewerApi, parsed);

    // Find the PBC toggle and state label
    const pbcToggle = document.getElementById('togglePBCHeader');
    const pbcStateLabel = document.getElementById('pbcStateLabel');

    // Toggle should reflect On, and state label if present too
    expect(pbcToggle).toBeTruthy();
    expect(pbcToggle.getAttribute('data-on')).toBe('true');
    if (pbcStateLabel) {
      expect(pbcStateLabel.textContent).toMatch(/On/i);
    }
  });
});
