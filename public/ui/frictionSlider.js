// Friction slider initialization for MD integrator.
// Controls window.__MLIP_CONFIG.mdFriction; default from constants. Range [0, 2.0] with 0.01 step.
import { DEFAULT_MD_FRICTION, DEFAULT_MIN_STEP_INTERVAL_MS } from '../util/constants.js';

export function initFrictionSlider({ hudEl }) {
  if (!hudEl) return;
  if (hudEl.querySelector('#frictionSliderWrapper')) return;

  if (typeof window !== 'undefined') {
    window.__MLIP_CONFIG = window.__MLIP_CONFIG || {
      minStepIntervalMs: DEFAULT_MIN_STEP_INTERVAL_MS,
      mdFriction: DEFAULT_MD_FRICTION,
    };
    if (window.__MLIP_CONFIG.minStepIntervalMs == null)
      window.__MLIP_CONFIG.minStepIntervalMs = DEFAULT_MIN_STEP_INTERVAL_MS;
    if (window.__MLIP_CONFIG.mdFriction == null)
      window.__MLIP_CONFIG.mdFriction = DEFAULT_MD_FRICTION;
  }

  const wrapper = document.createElement('div');
  wrapper.id = 'frictionSliderWrapper';
  wrapper.style.display = 'inline-flex';
  wrapper.style.alignItems = 'center';
  wrapper.style.marginLeft = '8px';
  wrapper.style.gap = '4px';
  wrapper.title = 'MD friction (dimensionless, UMA)';

  const label = document.createElement('span');
  label.id = 'frictionLabel';
  label.style.fontSize = '11px';
  label.style.opacity = '0.85';

  const slider = document.createElement('input');
  slider.type = 'range';
  slider.min = '0';
  slider.max = '2';
  slider.step = '0.01';
  slider.id = 'mdFrictionSlider';
  slider.style.width = '120px';

  function updateLabel(v) {
    label.textContent = `ζ=${Number(v).toFixed(2)}`;
  }
  function setFriction(v) {
    const num = Math.max(0, Math.min(5, Number(v)));
    if (typeof window !== 'undefined') {
      window.__MLIP_CONFIG.mdFriction = num;
    }
    updateLabel(num);
  }

  // Initialize from config
  const initial =
    typeof window !== 'undefined' && window.__MLIP_CONFIG?.mdFriction != null
      ? window.__MLIP_CONFIG.mdFriction
      : DEFAULT_MD_FRICTION;
  slider.value = String(initial);
  setFriction(initial);

  slider.addEventListener('input', () => setFriction(slider.value));

  wrapper.appendChild(label);
  wrapper.appendChild(slider);
  hudEl.appendChild(wrapper);

  return {
    getFriction: () => (typeof window !== 'undefined' ? window.__MLIP_CONFIG.mdFriction : initial),
    setFriction: (v) => {
      slider.value = String(v);
      setFriction(v);
    },
  };
}
