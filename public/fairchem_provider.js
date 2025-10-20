// HTTP provider removed. All force/energy retrieval must use WebSocket protobuf API.
export async function fairchemCalculate() {
  throw new Error('[REST removed] fairchemCalculate via HTTP removed. Use getWS().userInteraction() + getWS().waitForEnergy() after initSystem().');
}

export function createFairChemForcefield() {
  throw new Error('[REST removed] createFairChemForcefield removed. Use WS client in index.js for energy/forces.');
}
