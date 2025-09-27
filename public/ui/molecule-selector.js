// ui/molecule-selector.js - Simple molecule selection overlay
// Provides a function to create/show selector and reload page with chosen molecule

export async function showMoleculeSelector({ onClose } = {}) {
  // Reuse existing overlay if present
  if (document.getElementById('molSelectorOverlay')) {
    document.getElementById('molSelectorOverlay').style.display = 'flex';
    return;
  }
  const overlay = document.createElement('div');
  overlay.id = 'molSelectorOverlay';
  Object.assign(overlay.style, {
    position: 'fixed', top: 0, left: 0, right: 0, bottom: 0,
    background: 'rgba(0,0,0,0.65)', display: 'flex', flexDirection: 'column',
    alignItems: 'center', justifyContent: 'center', zIndex: 10000,
    fontFamily: 'sans-serif', color: '#fff'
  });
  const panel = document.createElement('div');
  Object.assign(panel.style, {
    background: '#1e1e1e', padding: '16px 20px', borderRadius: '8px',
    minWidth: '300px', maxWidth: '420px', maxHeight: '60%', overflowY: 'auto',
    boxShadow: '0 4px 24px rgba(0,0,0,0.4)'
  });
  const title = document.createElement('h2');
  title.textContent = 'Select Molecule';
  title.style.marginTop = '0';
  title.style.fontSize = '20px';
  panel.appendChild(title);

  const list = document.createElement('div');
  list.textContent = 'Loading...';
  panel.appendChild(list);

  const btnRow = document.createElement('div');
  Object.assign(btnRow.style, { display: 'flex', gap: '8px', marginTop: '16px', justifyContent: 'flex-end' });
  const closeBtn = document.createElement('button');
  closeBtn.textContent = 'Close';
  closeBtn.onclick = () => { overlay.style.display = 'none'; onClose && onClose(); };
  styleBtn(closeBtn);
  btnRow.appendChild(closeBtn);
  panel.appendChild(btnRow);

  overlay.appendChild(panel);
  document.body.appendChild(overlay);

  try {
    const res = await fetch('/api/molecules');
    if (!res.ok) throw new Error('Failed to fetch molecules');
    const data = await res.json();
    list.textContent = '';

    if (!data.molecules || !data.molecules.length) {
      list.textContent = 'No molecules found.';
    } else {
      data.molecules.forEach(m => {
        const item = document.createElement('button');
        item.textContent = m.name;
        item.title = m.file;
        styleBtn(item);
        item.style.width = '100%';
        item.style.textAlign = 'left';
        item.onclick = () => {
          const base = window.location.pathname;
          const params = new URLSearchParams(window.location.search);
          params.set('molecule', m.name); // use short name
          // Reload page with new parameter (retain vr.html if in VR)
          window.location.href = `${base}?${params.toString()}`;
        };
        list.appendChild(item);
      });
    }
  } catch (e) {
    list.textContent = 'Error loading molecules: ' + e.message;
  }
}

export function addMoleculeSelectButton(targetParent, { label = 'Molecules', small = false } = {}) {
  const btn = document.createElement('button');
  btn.textContent = label;
  styleBtn(btn);
  if (small) {
    btn.style.fontSize = '12px';
    btn.style.padding = '4px 8px';
  }
  btn.onclick = () => showMoleculeSelector();
  targetParent.appendChild(btn);
  return btn;
}

function styleBtn(btn) {
  Object.assign(btn.style, {
    background: '#333', color: '#fff', border: '1px solid #555',
    borderRadius: '4px', padding: '6px 12px', cursor: 'pointer',
    fontFamily: 'inherit', fontSize: '14px'
  });
  btn.onmouseenter = () => btn.style.background = '#444';
  btn.onmouseleave = () => btn.style.background = '#333';
  btn.onmousedown = () => btn.style.background = '#222';
  btn.onmouseup = () => btn.style.background = '#444';
}
