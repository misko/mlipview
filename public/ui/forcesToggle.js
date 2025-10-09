// Wires a HUD button to toggle force vectors and keep label/state in sync
// Contract:
// - Expects a button element with id 'btnToggleForces'
// - Updates button.textContent to 'forces on' or 'forces off'
// - Uses viewer.state.showForces and viewer.state.toggleForceVectorsVisibility()
export function initForcesToggle({ getViewer }){
  const btn = document.getElementById('btnToggleForces');
  if(!btn) return;
  function refresh(){
    const viewer = typeof getViewer === 'function' ? getViewer() : null;
    const on = !!(viewer && viewer.state && viewer.state.showForces);
    btn.textContent = on ? 'forces off' : 'forces on';
    btn.style.background = on ? '#284' : '';
  }
  btn.addEventListener('click', ()=>{
    const viewer = typeof getViewer === 'function' ? getViewer() : null;
    if(!viewer || !viewer.state || typeof viewer.state.toggleForceVectorsVisibility !== 'function') return;
    viewer.state.toggleForceVectorsVisibility();
    refresh();
  });
  // Also update when forcesChanged events are emitted (e.g., external toggles)
  try {
    const viewer = typeof getViewer === 'function' ? getViewer() : null;
    viewer?.state?.bus?.on?.('forcesChanged', refresh);
  } catch {}
  // When viewer becomes ready (index.html emits this), update label to match real state
  try { window.addEventListener('mlipviewer-ready', refresh); } catch {}
  refresh();
}
