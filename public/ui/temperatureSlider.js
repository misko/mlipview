// Temperature slider initialization.
// Provides a discrete 30-step slider from 0K to 2000K (inclusive) with emphasized ticks at 0, 298, 2000.
// Maintains target temperature in window.__MLIP_TARGET_TEMPERATURE (default 298) and mirrors to
// viewerApi.state.dynamics.targetTemperature (non-intrusive additional field) when viewer is available.

export function initTemperatureSlider({ hudEl, getViewer }) {
  if (!hudEl) return;
  // Avoid duplicate init
  if (hudEl.querySelector('#tempSliderWrapper')) return;

  // Build a lookup of 30 temperatures across 0-2000K; ensure 298K is exactly represented by replacing the closest value.
  const steps = 30;
  const MAX_K = 2000;
  let temps = Array.from({ length: steps }, (_, i) => Math.round(i * MAX_K / (steps - 1)));
  // Force include 298 by finding closest index and setting it; adjust neighbor if collision.
  let closestIdx = 0; let closestDiff = Infinity;
  for (let i = 0; i < temps.length; i++) {
    const d = Math.abs(temps[i] - 298);
    if (d < closestDiff) { closestDiff = d; closestIdx = i; }
  }
  temps[closestIdx] = 298;
  // De-duplicate sequence if forcing 298 created a duplicate with neighbor.
  for (let i = 1; i < temps.length; i++) {
    if (temps[i] === temps[i - 1]) temps[i] += 1; // small nudge; within 0-2000 still safe (rare occurrence)
  }

  // Wrapper UI
  const wrapper = document.createElement('div');
  wrapper.id = 'tempSliderWrapper';
  wrapper.style.display = 'inline-flex';
  wrapper.style.alignItems = 'center';
  wrapper.style.marginLeft = '8px';
  wrapper.style.gap = '4px';
  wrapper.title = 'Target MD temperature (K)';

  const label = document.createElement('span');
  label.id = 'tempLabel';
  label.style.fontSize = '11px';
  label.style.opacity = '0.85';
  label.textContent = 'T=298K';

  const slider = document.createElement('input');
  slider.type = 'range';
  slider.min = '0';
  slider.max = String(steps - 1);
  slider.step = '1';
  slider.value = String(closestIdx);
  slider.id = 'mdTempSlider';
  slider.style.width = '120px';

  // Tick labels (0, 298, 400) positioned under slider using a simple flex container.
  const ticks = document.createElement('div');
  ticks.style.position = 'relative';
  ticks.style.width = slider.style.width;
  ticks.style.height = '10px';
  ticks.style.fontSize = '9px';
  ticks.style.display = 'flex';
  ticks.style.justifyContent = 'space-between';
  ticks.style.pointerEvents = 'none';
  const tick0 = document.createElement('span'); tick0.textContent = '0K';
  const tickMid = document.createElement('span'); tickMid.textContent = '298K';
  const tickMax = document.createElement('span'); tickMax.textContent = '2000K';
  ticks.appendChild(tick0); ticks.appendChild(tickMid); ticks.appendChild(tickMax);

  function getTempForIndex(idx) { return temps[idx] ?? 298; }
  function updateTarget(idx) {
    const T = getTempForIndex(idx);
    window.__MLIP_TARGET_TEMPERATURE = T;
    label.textContent = `T=${T}K`;
    try {
      const v = getViewer && getViewer();
      if (v && v.state && v.state.dynamics) v.state.dynamics.targetTemperature = T;
    } catch {}
  }

  // Initialize global target
  if (typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE == null) {
    window.__MLIP_TARGET_TEMPERATURE = 298;
  }
  updateTarget(Number(slider.value));

  slider.addEventListener('input', () => updateTarget(Number(slider.value)));

  const col = document.createElement('div');
  col.style.display = 'flex';
  col.style.flexDirection = 'column';
  col.style.alignItems = 'stretch';
  col.appendChild(slider); col.appendChild(ticks);

  wrapper.appendChild(label);
  wrapper.appendChild(col);
  hudEl.appendChild(wrapper);

  return { getTemperature: () => window.__MLIP_TARGET_TEMPERATURE, setIndex: (i) => { slider.value = String(i); updateTarget(i); } };
}
