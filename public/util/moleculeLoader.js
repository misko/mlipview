import { parseXYZ, applyXYZToState } from './xyzLoader.js';
import { __count } from './funcCount.js';

export function getRequestedMoleculeFromUrl(win){
  try {
    const w = win || (typeof globalThis!=='undefined' && globalThis.window) || (typeof window!=='undefined' ? window : undefined);
    if(!w || !w.location) return null;
    const q = new URLSearchParams(w.location.search||'');
    const mol = q.get('mol');
    if(!mol) return null;
    // crude allowlist: only allow paths under molecules/
    if(/^molecules\/[\w.-]+\.xyz$/i.test(mol)) return mol;
  } catch {}
  return null;
}

function toAbsolute(url){
  try {
    // If already absolute, new URL will succeed and return as-is
    const origin = (typeof globalThis!=='undefined' && globalThis.window && globalThis.window.location && globalThis.window.location.origin)
      ? globalThis.window.location.origin
      : ((typeof window!=='undefined' && window.location && window.location.origin) ? window.location.origin : 'http://localhost');
    const u = new URL(url, origin);
    return u.toString();
  } catch { return url; }
}

export async function fetchXYZ(url) {
  __count('moleculeLoader#fetchXYZ');
  const href = toAbsolute(url);
  const res = await fetch(href);
  if (!res.ok) throw new Error(`Failed fetch ${url}: ${res.status}`);
  return await res.text();
}

export async function loadXYZIntoViewer(viewerApi, url) {
  __count('moleculeLoader#loadXYZIntoViewer');
  const txt = await fetchXYZ(url);
  const parsed = parseXYZ(txt);
  applyXYZToState(viewerApi.state, parsed);
  // Recompute bonds for new structure
  viewerApi.recomputeBonds();
  return parsed;
}

export async function loadDefault(viewerApi) {
  __count('moleculeLoader#loadDefault');
  // If URL requested a specific molecule, honor it first
  const requested = getRequestedMoleculeFromUrl();
  if (requested) {
    try {
      const parsed = await loadXYZIntoViewer(viewerApi, requested);
      return { file: requested, parsed };
    } catch (e) {
      console.warn('[moleculeLoader] requested mol failed', requested, e);
    }
  }
  // Default order: ROY then Benzene as fallback
  const candidates = ['molecules/roy.xyz','molecules/benzene.xyz'];
  for (const c of candidates) {
    try { const parsed = await loadXYZIntoViewer(viewerApi, c); return { file:c, parsed }; }
    catch (e) { console.warn('[moleculeLoader] failed', c, e); }
  }
  throw new Error('All default molecule loads failed');
}
