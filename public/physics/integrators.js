export function createDynamics(molState, { massPerElement = {}, timestep = 0.5 } = {}) {
  if (!molState.dynamics.velocities.length) {
    molState.dynamics.velocities = molState.positions.map(()=>({x:0,y:0,z:0}));
  }
  if (!molState.dynamics.mass.length) {
    molState.dynamics.mass = molState.elements.map(e => massPerElement[e] || 12);
  }
  function zeroForces() {
    molState.dynamics.forces = molState.positions.map(()=>({x:0,y:0,z:0}));
  }
  function applyForcesCallback(forceFn) {
    zeroForces();
    forceFn(molState);
  }
  function stepRelax({ gamma = 0.1, forceFn }) {
    applyForcesCallback(forceFn);
    for (let i=0;i<molState.positions.length;i++) {
      const f = molState.dynamics.forces[i];
      const p = molState.positions[i];
      p.x += gamma * f.x; p.y += gamma * f.y; p.z += gamma * f.z;
    }
    molState.markPositionsChanged();
    molState.markDynamicsChanged();
  }
  function stepMD({ targetTemp = 300, tau = 100, forceFn }) {
    applyForcesCallback(forceFn);
    const dt = timestep;
    for (let i=0;i<molState.positions.length;i++) {
      const v = molState.dynamics.velocities[i];
      const f = molState.dynamics.forces[i];
      const m = molState.dynamics.mass[i];
      const p = molState.positions[i];
      v.x += 0.5 * dt * f.x / m; v.y += 0.5 * dt * f.y / m; v.z += 0.5 * dt * f.z / m;
      p.x += dt * v.x; p.y += dt * v.y; p.z += dt * v.z;
    }
    applyForcesCallback(forceFn);
    let kinetic=0;
    for (let i=0;i<molState.positions.length;i++) {
      const v = molState.dynamics.velocities[i];
      const f = molState.dynamics.forces[i];
      const m = molState.dynamics.mass[i];
      v.x += 0.5 * dt * f.x / m; v.y += 0.5 * dt * f.y / m; v.z += 0.5 * dt * f.z / m;
      kinetic += 0.5 * m * (v.x*v.x + v.y*v.y + v.z*v.z);
    }
    const temp = (2*kinetic) / (3 * molState.positions.length);
    molState.dynamics.temperature = temp;
    const lambda = Math.sqrt(1 + dt/tau * ((targetTemp / temp) - 1));
    for (const v of molState.dynamics.velocities) { v.x*=lambda; v.y*=lambda; v.z*=lambda; }
    molState.markPositionsChanged();
    molState.markDynamicsChanged();
  }
  return { stepRelax, stepMD };
}
