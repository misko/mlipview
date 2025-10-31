interface DumpOptions {
  min?: number;
  sort?: 'asc' | 'desc';
  limit?: number | null;
}

type CountStore = Record<string, number>;

declare global {
  interface Window {
    __FUNC_COUNTS?: CountStore;
    dumpFunctionCounts?: (options?: DumpOptions) => Array<{ fn: string; count: number }>;
    resetFunctionCounts?: () => void;
  }
}

if (typeof window !== 'undefined') {
  window.__FUNC_COUNTS = window.__FUNC_COUNTS || Object.create(null);
  if (!window.dumpFunctionCounts) {
    window.dumpFunctionCounts = function dumpFunctionCounts(options: DumpOptions = {}) {
      const { min = 0, sort = 'desc', limit = null } = options;
      const rows = Object.entries(window.__FUNC_COUNTS as CountStore)
        .filter(([, count]) => count >= min)
        .sort((a, b) => (sort === 'asc' ? a[1] - b[1] : b[1] - a[1]));
      const out = rows.map(([fn, count]) => ({ fn, count }));
      return limit != null ? out.slice(0, limit) : out;
    };
    if (!window.resetFunctionCounts) {
      window.resetFunctionCounts = function resetFunctionCounts() {
        const store = window.__FUNC_COUNTS as CountStore;
        for (const key of Object.keys(store)) delete store[key];
      };
    }
  }
}

export function __count(key: string): void {
  if (typeof window === 'undefined') return;
  try {
    const store = window.__FUNC_COUNTS || (window.__FUNC_COUNTS = Object.create(null));
    store[key] = (store[key] || 0) + 1;
  } catch {
    // ignore failures in restricted environments
  }
}
