/** @jest-environment jsdom */

import { OMOL25_ELEMENTS } from '../public/data/periodicTable.js';

beforeEach(() => {
  document.body.innerHTML = '<div id="mount"></div>';
});

describe('periodic table OMol25 highlighting', () => {
  test('all OMol25 elements are decorated in the selection periodic table', async () => {
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    buildDesktopPanel({ attachTo: document.getElementById('mount') });

    const table = document.getElementById('miniPeriodic');
    expect(table).toBeTruthy();

    const lookup = (sym) => table.querySelector(`.pt-el[data-symbol="${sym}"]`);

    const missing = OMOL25_ELEMENTS.filter((sym) => !lookup(sym));
    expect(missing).toEqual([]);

    OMOL25_ELEMENTS.forEach((sym) => {
      const cell = lookup(sym);
      expect(cell.classList.contains('omol25')).toBe(true);
      expect(cell.getAttribute('title') || '').toContain(sym);
    });

    const nonOmol = ['Po', 'At', 'Rn'];
    nonOmol.forEach((sym) => {
      const cell = lookup(sym);
      expect(cell).toBeTruthy();
      expect(cell.classList.contains('omol25')).toBe(false);
    });
  });
});
