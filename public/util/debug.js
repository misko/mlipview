export function dbg(...args){ if (typeof process !== 'undefined' && process.env && process.env.DEBUG_PLOT) console.log('[dbg]', ...args); }
