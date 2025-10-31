// Molecule selection UI that triggers a full page reload with a URL param.
// Exports small helpers for easier unit testing.

/**
 * Build a URL with ?mol=<file> (preserving other query params) to force
 * a clean reload into the desired molecule.
 * @param {string} currentHref - current window.location.href
 * @param {string} file - e.g. 'molecules/water.xyz'
 * @returns {string} new href
 */
export function buildReloadUrl(currentHref, file) {
  const url = new URL(
    currentHref,
    typeof window !== 'undefined' ? window.location.origin : 'http://localhost'
  );
  url.searchParams.set('mol', file);
  return url.toString();
}

export function buildReloadUrlWithParam(currentHref, key, value) {
  const url = new URL(
    currentHref,
    typeof window !== 'undefined' ? window.location.origin : 'http://localhost'
  );
  url.searchParams.set(key, value);
  return url.toString();
}

// UTF-8 safe base64 helpers (browser + Node)
export function base64EncodeUtf8(str) {
  if (typeof Buffer !== 'undefined') return Buffer.from(str, 'utf8').toString('base64');
  // Browser path: encodeURIComponent to UTF-8 bytes, then map percent bytes to chars for btoa
  const utf8Bytes = encodeURIComponent(str).replace(/%([0-9A-F]{2})/g, (_, p1) =>
    String.fromCharCode(parseInt(p1, 16))
  );
  return btoa(utf8Bytes);
}
export function base64DecodeUtf8(b64) {
  if (typeof Buffer !== 'undefined') return Buffer.from(b64, 'base64').toString('utf8');
  const binary = atob(b64);
  const esc = Array.prototype.map
    .call(binary, (c) => '%' + ('00' + c.charCodeAt(0).toString(16)).slice(-2))
    .join('');
  return decodeURIComponent(esc);
}

/**
 * Install a <select id="moleculeSelect"> into hudEl that reloads the page on change.
 * Returns the created select element for testing.
 * @param {{ hudEl: HTMLElement, windowRef?: Window, documentRef?: Document }} opts
 */
export function installMoleculeSelector(opts) {
  if (!opts?.windowRef && typeof window !== 'undefined' && window.__MLIPVIEW_REACT_UI_ACTIVE) {
    return { selectEl: null };
  }
  const { hudEl } = opts;
  const windowRef = opts.windowRef || window;
  const documentRef = opts.documentRef || document;
  const select = documentRef.createElement('select');
  select.id = 'moleculeSelect';
  // Static list (served directly from /public) — keep labels short for UI
  select.innerHTML =
    '' +
    '<option value="molecules/roy.xyz">ROY</option>' +
    '<option value="molecules/benzene.xyz">Benzene</option>' +
    '<option value="molecules/water.xyz">Water</option>' +
    '<option value="molecules/water10.xyz">Water x10</option>' +
    '<option value="molecules/azomethane.xyz">Azomethane</option>' +
    '<option value="molecules/hydronium_hydroxide.xyz">Hydronium + Hydroxide</option>' +
    '<option value="molecules/sn2.xyz">SN2: Cl- + CH3I</option>' +
    '<option value="molecules/twomethane.xyz">Methyl radicals</option>' +
    '<option value="molecules/rna.xyz">RNA fragment</option>' +
    '<option value="molecules/imidazole.xyz">Imidazole</option>' +
    '<option value="molecules/naphthalene.xyz">Naphthalene</option>' +
    '<option value="molecules/benzoyl_peroxide.xyz">Benzoyl Peroxide</option>' +
    '<option value="molecules/urea.xyz">Urea</option>' +
    '<option value="molecules/acetic_acid.xyz">Acetic Acid</option>';
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
