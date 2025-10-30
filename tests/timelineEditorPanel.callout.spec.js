/** @jest-environment jsdom */

/**
 * Regression coverage for the timeline editor panel.
 * These tests ensure callout edits persist into the serialized control message payload.
 */

import { createTimelineEditorPanel } from '../public/ui/timelineEditorPanel.js';

function setupPanel(initialMessages, hooks = {}) {
  document.body.innerHTML = '';
  const host = document.createElement('div');
  host.id = 'app';
  document.body.appendChild(host);

  const recordedMessages = [];
  const recordedPlayback = [];

  const panel = createTimelineEditorPanel({
    attachTo: host,
    getMessages: () => JSON.parse(JSON.stringify(initialMessages)),
    setMessages: (msgs) => {
      recordedMessages.push(JSON.parse(JSON.stringify(msgs)));
    },
    refreshControlEngine: jest.fn(),
    getPlaybackConfig: () => (hooks.playbackSnapshot ? hooks.playbackSnapshot() : { defaultFps: 20 }),
    setPlaybackConfig: (cfg) => recordedPlayback.push(cfg),
    offsetToFrameId: () => 'frame-00000',
    offsetToIndex: () => 0,
    getCurrentOffset: () => -1,
  });

  return { host, panel, recordedMessages, recordedPlayback };
}

function clickButton(container, label) {
  const btn = Array.from(container.querySelectorAll('button')).find((node) => node.textContent === label);
  if (!btn) throw new Error(`Button "${label}" not found`);
  btn.click();
}

describe('timeline editor panel – callout persistence', () => {
  const baseMessage = [
    {
      id: 'control-1',
      priority: 0,
      label: 'Control 1',
      range: {
        start: { frameId: 'frame-00010' },
        end: { frameId: 'frame-00040', inclusive: true },
      },
      actions: [
        {
          type: 'overlay.callout',
          text: 'Initial text',
          anchor: { mode: 'atom', atoms: [2] },
          offset: [0, 1, 0],
          panelSize: { width: 2, height: 1 },
          textSize: 1.2,
        },
      ],
    },
  ];

  test('saving callout edits persists anchor atoms', () => {
    const { host, recordedMessages } = setupPanel(baseMessage);

    const atomInput = host.querySelector('input[placeholder="e.g. 5"]');
    expect(atomInput).toBeTruthy();
    atomInput.value = '7';
    atomInput.dispatchEvent(new Event('input', { bubbles: true }));

    clickButton(host, 'Save Changes');
    expect(recordedMessages).toHaveLength(1);
    const saved = recordedMessages[0][0];
    const callout = saved.actions.find((action) => action.type === 'overlay.callout');
    expect(callout).toBeTruthy();
    expect(callout.anchor).toBeDefined();
    expect(callout.anchor.mode).toBe('atom');
    expect(callout.anchor.atoms).toEqual([7]);
  });

  test('offset and panel fields persist after save', () => {
    const { host, recordedMessages } = setupPanel(baseMessage);

    const [offsetX, offsetY, offsetZ] = ['x', 'y', 'z'].map((placeholder) =>
      host.querySelector(`input[placeholder="${placeholder}"]`)
    );
    expect(offsetX).toBeTruthy();
    expect(offsetY).toBeTruthy();
    expect(offsetZ).toBeTruthy();
    offsetX.value = '0.5';
    offsetY.value = '1.5';
    offsetZ.value = '-0.25';
    offsetX.dispatchEvent(new Event('input', { bubbles: true }));
    offsetY.dispatchEvent(new Event('input', { bubbles: true }));
    offsetZ.dispatchEvent(new Event('input', { bubbles: true }));

    const textSizeField = Array.from(host.querySelectorAll('label'))
      .find((label) => label.textContent.includes('Text size'))
      ?.querySelector('input');
    expect(textSizeField).toBeTruthy();
    if (textSizeField) {
      textSizeField.value = '2.5';
      textSizeField.dispatchEvent(new Event('input', { bubbles: true }));
    }

    clickButton(host, 'Save Changes');
    expect(recordedMessages).toHaveLength(1);
    const saved = recordedMessages[0][0];
    const callout = saved.actions.find((action) => action.type === 'overlay.callout');
    expect(callout).toBeTruthy();
    expect(callout.offset).toEqual([0.5, 1.5, -0.25]);
    expect(callout.textSize).toBe(2.5);
  });
});
