export function createFairChemForceField() {
  throw new Error(
    '[REST removed] createFairChemForceField shim removed. Migrate callers to WebSocket flows.'
  );
}

export const legacyFairchemCalculate = () => {
  throw new Error('[REST removed] legacyFairchemCalculate removed.');
};
