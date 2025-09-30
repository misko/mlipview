// FairChem HTTP provider: fetch energy & forces from FastAPI server.
// Default behavior: use SAME origin via relative path '/simple_calculate' (served or reverse-proxied).
// Override precedence (highest -> lowest):
//   1. Explicit global window.__FAIRCHEM_URL (full base URL, no trailing slash required)
//   2. Fallback to relative '/simple_calculate'
// Prior hardcoded host: http://127.0.0.1:8002 (removed to avoid mixed-origin & CORS hassles in LAN / HTTPS dev).
const FC_BASE = (typeof window !== 'undefined' && window.__FAIRCHEM_URL) || null;

export async function fairchemCalculate(positions, elements, { properties = ['energy','forces'] } = {}) {
  // positions: array of [x,y,z]; elements: array of symbols
  const periodicTable = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'];
  const atomic_numbers = elements.map(el => {
    const Z = periodicTable.indexOf(el) + 1;
    if (Z <= 0) throw new Error('Unsupported element '+el);
    return Z;
  });
  const payload = { atomic_numbers, coordinates: positions, properties };
  const maxAttempts = 3;
  let attempt = 0; let lastErr = null;
  const endpoint = FC_BASE ? (FC_BASE.replace(/\/$/, '') + '/simple_calculate') : '/simple_calculate';
  if (typeof window !== 'undefined' && !window.__FAIRCHEM_FALLBACK_LOGGED) {
    window.__FAIRCHEM_FALLBACK_LOGGED = true;
    console.debug('[fairchem-provider:standalone]', 'endpoint =', endpoint, 'origin =', (typeof location!=='undefined'?location.origin:'(no window)'));
  }
  while (attempt < maxAttempts) {
    try {
      const started = performance && performance.now ? performance.now() : Date.now();
      console.debug('[fairchem-provider:standalone] request attempt', attempt+1, 'endpoint=', endpoint, 'atoms=', positions.length, 'props=', properties.join(','));
      const resp = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
      });
      if (!resp.ok) {
        // 4xx likely not transient; bail immediately except maybe 429, treat 500+ as retriable
        if (resp.status >= 500 && resp.status < 600 && attempt < maxAttempts - 1) {
          attempt++;
          await new Promise(r=>setTimeout(r, 100 * attempt));
          continue;
        }
        throw new Error('FairChem HTTP error '+resp.status);
      }
      const data = await resp.json();
      const elapsed = (performance && performance.now ? performance.now() : Date.now()) - started;
      console.debug('[fairchem-provider:standalone] success attempt', attempt+1, 'elapsed_ms=', elapsed.toFixed(1), 'energy=', (data?.results?.energy!=null?data.results.energy:'?'));
      const res = data.results || {};
      return { energy: res.energy, forces: res.forces };
    } catch (e) {
      lastErr = e;
      console.debug('[fairchem-provider:standalone] failure attempt', attempt+1, e?.message||e);
      // Network errors or thrown above: retry if attempts remain
      if (attempt < maxAttempts - 1) {
        attempt++;
        await new Promise(r=>setTimeout(r, 100 * attempt));
        continue;
      }
      break;
    }
  }
  throw lastErr || new Error('FairChem calculation failed');
}

export function createFairChemForcefield(state){
  return {
    /**
     * Compute energy & forces.
     * @param {Array<[number,number,number]>} overridePos Optional array of xyz used instead of state.positions.
     */
    async computeForces(overridePos){
      const pos = overridePos ? overridePos : state.positions.map(p=>[p.x,p.y,p.z]);
      const { energy, forces } = await fairchemCalculate(pos, state.elements);
      // If we used override positions (line-search trial), do not mutate state positions here; caller handles commit.
      if (!overridePos) {
        // Commit energy to state only when evaluating current accepted geometry
        state.dynamics = state.dynamics || {}; state.dynamics.energy = energy;
      }
      return { energy, forces };
    }
  };
}
