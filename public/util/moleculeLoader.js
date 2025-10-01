import { parseXYZ, applyXYZToState } from './xyzLoader.js';
import { __count } from './funcCount.js';

export async function fetchXYZ(url) {
  __count('moleculeLoader#fetchXYZ');
  const res = await fetch(url);
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
  // Attempt roy.xyz then fallback benzene.xyz
  const candidates = ['molecules/roy.xyz','molecules/benzene.xyz'];
  for (const c of candidates) {
    try { const parsed = await loadXYZIntoViewer(viewerApi, c); return { file:c, parsed }; }
    catch (e) { console.warn('[moleculeLoader] failed', c, e); }
  }
  throw new Error('All default molecule loads failed');
}
