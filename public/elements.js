import { getElementInfo } from './data/periodicTable.js';

// Backwards compatible surface for existing code
export function elInfo(sym) {
  return getElementInfo(sym);
}

// Named export for completeness if some code still imports ELEMENTS
export const ELEMENTS = new Proxy(
  {},
  {
    get: (_, sym) => getElementInfo(sym),
  }
);
