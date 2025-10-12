/** @jest-environment jsdom */

import { buildDesktopPanel } from '../public/ui/desktopPanel.js';

describe('mobile top bar tabs', () => {
	function setup() {
		document.body.innerHTML = '<div id="app"></div>';
		const host = document.getElementById('app');
		buildDesktopPanel({ attachTo: host });
		const panel = document.getElementById('controlPanel');
		const bar = document.getElementById('mobileTopBar');
		return { host, panel, bar };
	}

	test('renders three tabs at top', () => {
		const { bar } = setup();
		expect(bar).toBeTruthy();
		const buttons = Array.from(bar.querySelectorAll('button.tab')).map(b=>b.textContent.trim());
		expect(buttons).toEqual(['Live Metrics','Simulation','System']);
	});

	test('only one section expanded at a time', () => {
		const { panel, bar } = setup();
		const btnLive = bar.querySelector('#mobileTab-live');
		const btnSim = bar.querySelector('#mobileTab-simulation');

		// Open Live
		btnLive.click();
		expect(panel.getAttribute('data-mobile-open')).toBe('true');
		const active1 = panel.querySelector('.panel-section[data-mobile-active="true"]');
		expect(active1 && active1.id).toBe('section-live-stats');

		// Switch to Simulation
		btnSim.click();
		const active2 = panel.querySelector('.panel-section[data-mobile-active="true"]');
		expect(active2 && active2.id).toBe('section-simulation');
	});

	test('tapping active tab collapses the sheet', () => {
		const { panel, bar } = setup();
		const btnSys = bar.querySelector('#mobileTab-system');

		btnSys.click();
		expect(panel.getAttribute('data-mobile-open')).toBe('true');
		let active = panel.querySelector('.panel-section[data-mobile-active="true"]');
		expect(active && active.id).toBe('section-system');

		// Tap again to collapse
		btnSys.click();
		expect(panel.getAttribute('data-mobile-open')).toBe(null);
		active = panel.querySelector('.panel-section[data-mobile-active="true"]');
		expect(active).toBe(null);
	});
});
