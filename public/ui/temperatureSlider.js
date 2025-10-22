// Temperature slider initialization.
// Provides a discrete 30-step slider from 0K to 3000K (inclusive) with emphasized ticks at 298, 400, 3000 (0K tick removed to avoid label overlap).
// Maintains target temperature in window.__MLIP_TARGET_TEMPERATURE (default 298) and mirrors to
// viewerApi.state.dynamics.targetTemperature (non-intrusive additional field) when viewer is available.

function buildTemps() {
  const steps = 30;
  const MAX_K = 3000;
  let temps = Array.from({ length: steps }, (_, i) => Math.round((i * MAX_K) / (steps - 1)));
  // Force include 298 by finding closest index and setting it; adjust neighbor if collision.
  let idx298 = 0;
  let closestDiff = Infinity;
  for (let i = 0; i < temps.length; i++) {
    const d = Math.abs(temps[i] - 298);
    if (d < closestDiff) {
      closestDiff = d;
      idx298 = i;
    }
  }
  temps[idx298] = 298;
  for (let i = 1; i < temps.length; i++) {
    if (temps[i] === temps[i - 1]) temps[i] += 1;
  }
  // Force include 1500
  let idx1500 = 0;
  let diff1500 = Infinity;
  for (let i = 0; i < temps.length; i++) {
    const d = Math.abs(temps[i] - 1500);
    if (d < diff1500) {
      diff1500 = d;
      idx1500 = i;
    }
  }
  temps[idx1500] = 1500;
  for (let i = 1; i < temps.length; i++) {
    if (temps[i] === temps[i - 1]) temps[i] += 1;
  }
  return { temps, idx1500, MAX_K };
}

export function initTemperatureSlider({ hudEl, getViewer }) {
  if (!hudEl) return;
  // Avoid duplicate init
  if (hudEl.querySelector('#tempSliderWrapper')) return;
  const { temps, idx1500, MAX_K } = buildTemps();

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
  // Keep width stable when value length changes
  label.style.whiteSpace = 'pre';
  try { label.style.fontVariantNumeric = 'tabular-nums'; } catch {}
  label.textContent = 'T=    1500K';

  const slider = document.createElement('input');
  slider.type = 'range';
  slider.min = '0';
  slider.max = String(temps.length - 1);
  slider.step = '1';
  slider.value = String(idx1500);
  slider.id = 'mdTempSlider';
  slider.style.width = '120px';

  // Tick labels (298 K, 400 K, 3000 K) positioned under slider using absolute positioning for accurate value mapping.
  const ticks = document.createElement('div');
  ticks.style.position = 'relative';
  ticks.style.width = slider.style.width;
  ticks.style.height = '10px';
  ticks.style.fontSize = '9px';
  ticks.style.display = 'block';
  ticks.style.pointerEvents = 'none';
  function addTick(val, labelText) {
    const span = document.createElement('span');
    span.textContent = labelText;
    span.style.position = 'absolute';
    span.style.transform = 'translateX(-50%)';
    const pct = (val / MAX_K) * 100;
    span.style.left = pct + '%';
    ticks.appendChild(span);
  }
  addTick(298, '298K');
  //addTick(400, '400K');
  addTick(MAX_K, MAX_K + 'K');

  function getTempForIndex(idx) {
    return temps[idx] ?? 1500;
  }
  function updateTarget(idx) {
    const T = getTempForIndex(idx);
    window.__MLIP_TARGET_TEMPERATURE = T;
    // Pad numeric portion to fixed width (8 digits) so layout doesn't shift
    const tStr = String(Math.round(T)).padStart(8, ' ');
    label.textContent = `T=${tStr}K`;
    try {
      const v = getViewer && getViewer();
      if (v && v.state && v.state.dynamics) v.state.dynamics.targetTemperature = T;
    } catch { }
    // Notify listeners (viewer can live-update running MD params)
    try {
      const evt = new Event('mlip:temperature-changed');
      window.dispatchEvent(evt);
    } catch { }
  }

  // Initialize from existing global if provided by loader (XYZ temperature), else default to 1500
  if (typeof window !== 'undefined') {
    const preset =
      window.__MLIP_TARGET_TEMPERATURE != null ? Number(window.__MLIP_TARGET_TEMPERATURE) : null;
    if (preset == null || !Number.isFinite(preset)) {
      window.__MLIP_TARGET_TEMPERATURE = 1500;
      slider.value = String(idx1500);
    } else {
      // Find nearest slider index for pre-set temperature
      let bestI = 0,
        bestD = Infinity;
      const want = preset;
      for (let i = 0; i < temps.length; i++) {
        const d = Math.abs(temps[i] - want);
        if (d < bestD) {
          bestD = d;
          bestI = i;
        }
      }
      slider.value = String(bestI);
    }
  }
  updateTarget(Number(slider.value));

  slider.addEventListener('input', () => updateTarget(Number(slider.value)));

  const col = document.createElement('div');
  col.style.display = 'flex';
  col.style.flexDirection = 'column';
  col.style.alignItems = 'stretch';
  col.appendChild(slider);
  col.appendChild(ticks);

  wrapper.appendChild(label);
  wrapper.appendChild(col);
  hudEl.appendChild(wrapper);

  // Mark widget so a global sync can find and update it later
  wrapper.setAttribute('data-mlip-temp-widget', 'true');

  return {
    getTemperature: () => window.__MLIP_TARGET_TEMPERATURE,
    setIndex: (i) => {
      slider.value = String(i);
      updateTarget(i);
    },
  };
}

// Synchronize an existing temperature slider widget to the current global temperature
export function syncTemperatureSliderToGlobal(root = document) {
  try {
    const wrapper = root.querySelector('#tempSliderWrapper');
    if (!wrapper) return;
    const slider = wrapper.querySelector('#mdTempSlider');
    const label = wrapper.querySelector('#tempLabel');
    if (!slider || !label) return;
    const { temps } = buildTemps();
    const want =
      typeof window !== 'undefined' && window.__MLIP_TARGET_TEMPERATURE != null
        ? Number(window.__MLIP_TARGET_TEMPERATURE)
        : null;
    if (want == null || !Number.isFinite(want)) return;
    let bestI = 0,
      bestD = Infinity;
    for (let i = 0; i < temps.length; i++) {
      const d = Math.abs(temps[i] - want);
      if (d < bestD) {
        bestD = d;
        bestI = i;
      }
    }
    slider.value = String(bestI);
    label.textContent = `T=${want}K`;
  } catch { }
}

// Install a single global listener the first time this module is evaluated
try {
  if (typeof window !== 'undefined' && !window.__MLIP_TEMP_SLIDER_EVT) {
    window.__MLIP_TEMP_SLIDER_EVT = true;
    document.addEventListener('mlip:temperature-changed', () => {
      try {
        syncTemperatureSliderToGlobal(document);
      } catch { }
    });
  }
} catch { }
