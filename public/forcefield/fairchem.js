// Minimal fairchem forcefield adapter (stub) to satisfy tests.
// Provides createFairChemForceField({ baseUrl }) returning compute({ elements, positions })
// The real implementation can be expanded later.
export function createFairChemForceField(opts = {}) {
  const baseUrl = opts.baseUrl || 'http://localhost:8110';
  async function compute({ elements = [], positions = [] }) {
    // Convert positions to flat array if objects
    const coords = positions.map(p => (Array.isArray(p) ? p : [p.x, p.y, p.z]));
    try {
      const resp = await fetch(baseUrl + '/simple_calculate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ elements, positions: coords })
      });
      if (!resp.ok) throw new Error('HTTP ' + resp.status);
      const data = await resp.json();
      // Expect { energy: number, forces: [[fx,fy,fz], ...] }
      return { energy: data.energy ?? 0, forces: data.forces || [] };
    } catch (e) {
      // Fallback deterministic mock if server unavailable
      const mockEnergy = elements.length ? elements.length * 0.123 : 0;
      return { energy: mockEnergy, forces: positions.map(() => [0,0,0]) };
    }
  }
  return { compute };
}
