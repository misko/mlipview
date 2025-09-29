// mlipviewer2/public/ui/adapters/rotationHUDDesktop.js
// Desktop DOM implementation of bond rotation HUD.

export function createRotationHUDDesktop({ controller, selectionService }) {
  const root = document.createElement('div');
  root.style.position = 'absolute';
  root.style.bottom = '14px';
  root.style.left = '50%';
  root.style.transform = 'translateX(-50%)';
  root.style.display = 'flex';
  root.style.gap = '8px';
  root.style.padding = '8px 10px';
  root.style.background = 'rgba(15,18,24,0.72)';
  root.style.border = '1px solid rgba(255,255,255,0.08)';
  root.style.borderRadius = '10px';
  root.style.font = '13px system-ui,Segoe UI,Roboto,sans-serif';
  root.style.color = '#d7e6ff';
  root.style.backdropFilter = 'blur(4px)';

  const label = document.createElement('span');
  label.textContent = 'Select a bond…';
  label.style.minWidth = '180px';

  const btnMinus = document.createElement('button');
  btnMinus.textContent = '⟲ −';
  const btnPlus = document.createElement('button');
  btnPlus.textContent = '⟳ +';
  const btnCycle = document.createElement('button');
  btnCycle.textContent = '↺ ori';

  [btnMinus, btnPlus, btnCycle].forEach(b => {
    b.style.padding = '8px 12px';
    b.style.fontSize = '16px';
    b.style.background = 'rgba(255,255,255,0.08)';
    b.style.border = '1px solid rgba(255,255,255,0.18)';
    b.style.color = '#f2fbff';
    b.style.borderRadius = '6px';
    b.style.cursor = 'pointer';
  });

  btnMinus.onclick = () => controller.rotateMinus();
  btnPlus.onclick = () => controller.rotatePlus();
  btnCycle.onclick = () => controller.cycleOrientation();

  root.append(label, btnMinus, btnPlus, btnCycle);
  document.body.appendChild(root);

  function refresh() {
    const sel = selectionService.getSelection();
    if (!sel || sel.kind !== 'bond') {
      label.textContent = 'Select a bond…';
      btnMinus.disabled = btnPlus.disabled = btnCycle.disabled = true;
      root.style.opacity = '0.55';
    } else {
      const { i, j, orientation } = sel.data;
      label.textContent = `Bond ${i}-${j} (ori ${orientation})`; // step display optional
      btnMinus.disabled = btnPlus.disabled = btnCycle.disabled = false;
      root.style.opacity = '1';
    }
  }
  refresh();
  selectionService.on('selectionChanged', refresh);

  function dispose() { root.remove(); }
  return { dispose, root };
}
