/** @jest-environment jsdom */

import { initTemperatureSlider } from '../public/ui/temperatureSlider.js';

describe('x-temperature-tick-position', () => {
  test('298K tick is positioned proportionally for 0-3000K range', () => {
    const hud = document.createElement('div');
    document.body.appendChild(hud);
    initTemperatureSlider({ hudEl: hud, getViewer: () => null });
    const wrapper = hud.querySelector('#tempSliderWrapper');
    expect(wrapper).toBeTruthy();
    const tickContainer = wrapper.querySelector('div > div > div');
    const tick = Array.from(tickContainer.querySelectorAll('span')).find(
      (el) => el.textContent === '298K'
    );
    expect(tick).toBeTruthy();
    const leftPercent = parseFloat(tick.style.left.replace('%', ''));
    const expectedPercent = (298 / 3000) * 100;
    expect(Math.abs(leftPercent - expectedPercent)).toBeLessThan(0.5);
  });
});
