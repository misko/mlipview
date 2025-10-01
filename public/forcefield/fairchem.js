// Consolidated FairChem shim: legacy createFairChemForceField now delegates to modern fairchem_provider.
// This preserves older test imports while avoiding duplicate HTTP logic.
export { fairchemCalculate as legacyFairchemCalculate } from '../fairchem_provider.js';
import { fairchemCalculate } from '../fairchem_provider.js';

// Legacy API signature wrapper
export function createFairChemForceField(opts = {}) {
  const baseUrl = (opts.baseUrl || null);
  // If a custom baseUrl is provided, temporarily override via direct call; otherwise defer to provider logic.
  return {
    async compute({ elements = [], positions = [] }) {
      const coords = positions.map(p => (Array.isArray(p) ? p : [p.x, p.y, p.z]));
      if (baseUrl) {
        try {
          const resp = await fetch(baseUrl.replace(/\/$/, '') + '/simple_calculate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ elements, positions: coords, properties: ['energy','forces'] })
          });
          if (!resp.ok) throw new Error('HTTP ' + resp.status);
          const data = await resp.json();
          return { energy: data?.results?.energy ?? data.energy ?? 0, forces: data?.results?.forces || data.forces || [] };
        } catch (e) {
          const mockEnergy = elements.length ? elements.length * 0.123 : 0;
          return { energy: mockEnergy, forces: positions.map(() => [0,0,0]) };
        }
      }
      // Fallback to canonical provider behavior (state-less form)
      const res = await fairchemCalculate(coords, elements, { properties:['energy','forces'] });
      return { energy: res.energy, forces: res.forces };
    }
  };
}
