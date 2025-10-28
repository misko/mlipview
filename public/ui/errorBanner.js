// Lightweight top-of-screen error banner for transient messages
// Usage: showErrorBanner('Message', { timeoutMs: 4000 })

let __bannerEl = null;
let __hideTimer = null;

function ensureBanner() {
  if (__bannerEl && __bannerEl.parentNode) return __bannerEl;
  const host =
    typeof document !== 'undefined' ? document.getElementById('app') || document.body : null;
  if (!host) return null;
  const el = document.createElement('div');
  el.id = 'mlip-top-error-banner';
  el.style.position = 'absolute';
  el.style.top = '0';
  el.style.left = '0';
  el.style.right = '0';
  el.style.zIndex = '10000';
  el.style.padding = '8px 12px';
  el.style.textAlign = 'center';
  el.style.fontSize = '13px';
  el.style.background = '#fde2e2'; // light red
  el.style.color = '#5a0000';
  el.style.borderBottom = '1px solid #f5b5b5';
  el.style.boxShadow = '0 2px 6px rgba(0,0,0,0.1)';
  el.style.transform = 'translateY(-100%)';
  el.style.transition = 'transform 150ms ease-in-out, opacity 150ms ease-in-out';
  el.style.opacity = '0.95';
  host.appendChild(el);
  __bannerEl = el;
  return __bannerEl;
}

export function hideErrorBanner() {
  if (!__bannerEl) return;
  try {
    __bannerEl.style.transform = 'translateY(-100%)';
    __bannerEl.style.opacity = '0.0';
    if (__hideTimer) clearTimeout(__hideTimer);
  } catch {}
  __hideTimer = null;
}

export function showErrorBanner(message, opts = {}) {
  try {
    const el = ensureBanner();
    if (!el) return;
    if (opts.html) {
      el.innerHTML = opts.html;
    } else {
      el.textContent = String(message || 'Error');
    }
    if (opts.className) {
      el.className = opts.className;
    } else {
      el.className = '';
    }
    // Show (slide down)
    el.style.transform = 'translateY(0)';
    el.style.opacity = '0.98';
    if (typeof opts.onRender === 'function') {
      try { opts.onRender(el); } catch {}
    }
    const persist = !!opts.persist;
    const timeoutMs = Number.isFinite(opts.timeoutMs) ? opts.timeoutMs : 4000;
    if (__hideTimer) { clearTimeout(__hideTimer); __hideTimer = null; }
    if (!persist) {
      __hideTimer = setTimeout(
        () => {
          try {
            el.style.transform = 'translateY(-100%)';
            el.style.opacity = '0.0';
          } catch {}
        },
        Math.max(500, timeoutMs)
      );
    } else {
      __hideTimer = null;
    }
  } catch {}
}

// Convenience for other severities if needed later
export function showInfoBanner(message, opts = {}) {
  showErrorBanner(message, opts);
}

// Expose globally for ad-hoc use/debug
try {
  if (typeof window !== 'undefined') window.showErrorBanner = showErrorBanner;
} catch {}
