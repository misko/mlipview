// Molecule selection UI that triggers a full page reload with a URL param.
// Exports small helpers for easier unit testing.

/**
 * Build a URL with ?mol=<file> (preserving other query params) to force
 * a clean reload into the desired molecule.
 * @param {string} currentHref - current window.location.href
 * @param {string} file - e.g. 'molecules/water.xyz'
 * @returns {string} new href
 */
export function buildReloadUrl(currentHref, file){
  const url = new URL(currentHref, typeof window!== 'undefined' ? window.location.origin : 'http://localhost');
  url.searchParams.set('mol', file);
  return url.toString();
}

/**
 * Install a <select id="moleculeSelect"> into hudEl that reloads the page on change.
 * Returns the created select element for testing.
 * @param {{ hudEl: HTMLElement, windowRef?: Window, documentRef?: Document }} opts
 */
export function installMoleculeSelector(opts){
  const { hudEl } = opts;
  const windowRef = opts.windowRef || window;
  const documentRef = opts.documentRef || document;
  const select = documentRef.createElement('select');
  select.id = 'moleculeSelect';
  select.innerHTML = '<option value="molecules/roy.xyz">ROY</option>'+
                     '<option value="molecules/benzene.xyz">Benzene</option>'+
                     '<option value="molecules/water.xyz">Water</option>';
  select.onchange = () => {
    const file = select.value;
    const href = buildReloadUrl(windowRef.location.href, file);
    // Full navigation replaces current page; avoids any retained in-memory state
    if (typeof windowRef.location.assign === 'function') {
      windowRef.location.assign(href);
    } else {
      // Fallback for tests
      windowRef.location.href = href;
    }
  };
  hudEl.appendChild(select);
  return { selectEl: select };
}
