import { getElementInfo } from '../../data/periodicTable.js';

export function covalentRadius(symbol) {
  const info = getElementInfo(symbol);
  if (!info) return 0.8;
  if (typeof info.covRad === 'number' && Number.isFinite(info.covRad)) return info.covRad;
  return 0.8;
}
