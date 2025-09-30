// FairChem HTTP provider: fetch energy & forces from local FastAPI server
// Uses /simple_calculate for efficiency (atomic numbers + coordinates)

const FC_URL = (typeof window !== 'undefined' && window.__FAIRCHEM_URL) || 'http://127.0.0.1:8000';

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
  while (attempt < maxAttempts) {
    try {
      const resp = await fetch(FC_URL + '/simple_calculate', {
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
      const res = data.results || {};
      return { energy: res.energy, forces: res.forces };
    } catch (e) {
      lastErr = e;
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
