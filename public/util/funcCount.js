// Lightweight function call counter instrumentation helper.
// Usage: import { __count } from './util/funcCount.js'; then call __count('file#fn');
// Adds global window.__FUNC_COUNTS and window.dumpFunctionCounts().

if (typeof window !== 'undefined') {
  window.__FUNC_COUNTS = window.__FUNC_COUNTS || Object.create(null);
  if (!window.dumpFunctionCounts) {
    window.dumpFunctionCounts = function(options={}) {
      const { min=0, sort='desc', limit=null } = options;
      const rows = Object.entries(window.__FUNC_COUNTS)
        .filter(([k,v]) => v >= min)
        .sort((a,b)=> sort==='asc' ? a[1]-b[1] : b[1]-a[1]);
      const out = rows.map(([k,v])=>({ fn:k, count:v }));
      if (limit != null) return out.slice(0, limit);
      return out;
    };
    if (!window.resetFunctionCounts) {
      window.resetFunctionCounts = function(){ for (const k of Object.keys(window.__FUNC_COUNTS)) delete window.__FUNC_COUNTS[k]; };
    }
  }
}
export function __count(key){
  if (typeof window === 'undefined') return;
  try {
    const store = window.__FUNC_COUNTS || (window.__FUNC_COUNTS = Object.create(null));
    store[key] = (store[key] || 0) + 1;
  } catch {}
}
