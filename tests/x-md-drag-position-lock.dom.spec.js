/**
 * @jest-environment jsdom
 */

import { getWS } from '../public/fairchem_ws_client.js';
import { stubWebSocketAndHook } from './utils/wsTestStub.js';

describe('x-md-drag-position-lock', () => {
  test('dragging atom emits USER_INTERACTION while md frames arrive', async () => {
    const { sent, emit } = stubWebSocketAndHook();
    const ws = getWS();
    await ws.ensureConnected();

    ws.initSystem({
      atomic_numbers: [1, 1],
      positions: [
        [0, 0, 0],
        [0.74, 0, 0],
      ],
    });

    ws.userInteraction({
      positions: [
        [0.3, 0, 0],
        [0.74, 0, 0],
      ],
    });

    emit({
      seq: 2,
      positions: [
        [0.25, 0, 0],
        [0.74, 0, 0],
      ],
      simStep: 2,
    });

    const ui = sent.find((m) => m && (m.type === 'USER_INTERACTION' || m.type === 6));
    expect(ui).toBeTruthy();
  });
});

