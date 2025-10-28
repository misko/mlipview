// Single Cell toggle control: toggles both cell box visibility and ghost cells.
// Usage: initCellToggle({ getViewer })

function labelFor(state) {
  const on = !!(state && (state.showCell || state.showGhostCells));
  return on ? 'Cell OFF' : 'Cell ON';
}

export function initCellToggle({ getViewer }) {
  const btn = document.getElementById('btnCell');
  if (!btn) return;
  function refresh() {
    try {
      const v = getViewer && getViewer();
      const st = v && v.state;
      btn.textContent = labelFor(st);
    } catch {}
  }
  btn.addEventListener('click', () => {
    const v = getViewer && getViewer();
    if (!v) return;
    try {
      const st = v.state;
      // Toggle cell visibility (supports legacy or enhanced)
      (st.toggleCellVisibilityEnhanced || st.toggleCellVisibility).call(st);
      // Ensure ghost cells match cell visibility (unified behavior)
      // If cell is now shown, enable ghosts; else disable ghosts.
      const wantGhosts = !!st.showCell;
      if (!!st.showGhostCells !== wantGhosts) st.toggleGhostCells();
      document.getElementById('status').textContent = st.showCell ? 'Cell ON' : 'Cell OFF';
    } finally {
      refresh();
    }
  });
  // Initial label
  refresh();
}
