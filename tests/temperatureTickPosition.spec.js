/** @jest-environment jsdom */
import { initTemperatureSlider } from '../public/ui/temperatureSlider.js';

test('298K tick positioned proportionally for 0-3000K range', () => {
  const hud = document.createElement('div');
  document.body.appendChild(hud);
  initTemperatureSlider({ hudEl: hud, getViewer: () => null });
  const wrapper = hud.querySelector('#tempSliderWrapper');
  expect(wrapper).toBeTruthy();
  const tickContainer = wrapper.querySelector('div > div > div'); // wrapper -> col -> ticks
  // Fallback: search by text
  const span298 = Array.from(tickContainer.querySelectorAll('span')).find(
    (s) => s.textContent === '298K'
  );
  expect(span298).toBeTruthy();
  const left = parseFloat(span298.style.left.replace('%', ''));
  const expected = (298 / 3000) * 100; // ~9.9333%
  // Allow small rounding differences (<0.5%)
  expect(Math.abs(left - expected)).toBeLessThan(0.5);
});
