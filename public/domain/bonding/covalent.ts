import { getElementInfo } from '../../data/periodicTable.js';

interface ElementInfo {
  covRad?: number;
}

export function covalentRadius(symbol: string): number {
  const info = getElementInfo(symbol) as ElementInfo | undefined;
  if (info && typeof info.covRad === 'number' && Number.isFinite(info.covRad)) {
    return info.covRad;
  }
  return 0.8;
}
