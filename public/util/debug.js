// public/util/debug.js
// Centralized debug flag & logging helper.
// Enable by setting window.DEBUG = true before scripts load or adding ?debug=1 to URL.

const DEBUG = (() => {
  try {
    if (typeof window !== 'undefined') {
      if (typeof window.DEBUG === 'boolean') return window.DEBUG;
      const params = new URLSearchParams(window.location.search);
      return params.get('debug') === '1';
    }
  } catch {}
  return false;
})();

function dbg(...args) {
  if (DEBUG) {
    try { console.debug(...args); } catch {}
  }
}

export { DEBUG, dbg };
