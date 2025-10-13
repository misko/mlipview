import { parseXYZ, applyXYZToState } from './xyzLoader.js';
import { __count } from './funcCount.js';
import { validateParsedXYZ } from './constraints.js';
import { smilesToXYZ, isLikelySmiles } from './smilesLoader.js';
import { base64DecodeUtf8 } from '../ui/moleculeSelect.js';

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

export function getRequestedMolXYZFromUrl(win){
  try {
    const w = win || (typeof globalThis!=='undefined' && globalThis.window) || (typeof window!=='undefined' ? window : undefined);
    if(!w || !w.location) return null;
    const q = new URLSearchParams(w.location.search||'');
    const b64 = q.get('molxyz');
    if(!b64) return null;
    return b64;
  } catch {}
  return null;
}

export function getRequestedSmilesFromUrl(win){
  try {
    const w = win || (typeof globalThis!=='undefined' && globalThis.window) || (typeof window!=='undefined' ? window : undefined);
    if(!w || !w.location) return null;
    const q = new URLSearchParams(w.location.search||'');
    const s = q.get('smiles');
    if(!s) return null;
    if(isLikelySmiles(s)) return s;
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
  // Reuse common path so initial positions get cached consistently
  return applyParsedToViewer(viewerApi, parsed);
}

export function applyParsedToViewer(viewerApi, parsed){
  const valid = validateParsedXYZ(parsed);
  if (!valid.ok) throw new Error(valid.error || 'Validation failed');
  applyXYZToState(viewerApi.state, parsed);
  // Cache the freshly loaded positions as the reset baseline for VR/AR
  try {
    const st = viewerApi.state;
    if (Array.isArray(st?.positions)) {
      st.__initialPositions = st.positions.map(p=>({ x:p.x, y:p.y, z:p.z }));
      try { console.log('[loader] cached initial positions', { count: st.__initialPositions.length }); } catch {}
    }
  } catch {}
  viewerApi.recomputeBonds();
  return parsed;
}

export async function loadDefault(viewerApi) {
  __count('moleculeLoader#loadDefault');
  // Priority: inline ?molxyz, then ?smiles, then static ?mol
  try {
    const b64 = getRequestedMolXYZFromUrl();
    if (b64) {
      const txt = base64DecodeUtf8(b64);
      const parsed = parseXYZ(txt);
      applyParsedToViewer(viewerApi, parsed);
      return { file: 'inline.xyz', parsed };
    }
  } catch (e) { console.warn('[moleculeLoader] inline molxyz failed', e); }
  try {
    const smi = getRequestedSmilesFromUrl();
    if (smi) {
      const xyzText = await smilesToXYZ(smi);
      const parsed = parseXYZ(xyzText);
      applyParsedToViewer(viewerApi, parsed);
      return { file: `smiles:${smi}`, parsed };
    }
  } catch (e) { console.warn('[moleculeLoader] smiles load failed', e); }
  // If URL requested a specific molecule file, honor it
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
