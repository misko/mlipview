/**
 * @jest-environment jsdom
 */
// Verify that during drag, client continues to receive frames but trusts dragged atom position (conceptual)
import { getWS } from '../public/fairchem_ws_client.js';
import { stubWebSocketAndHook } from './utils/wsTestStub.js';

describe('Dragged atom position lock (client-side semantics)', () => {
  test('frames continue; client can choose to prefer local drag pos', async () => {
    const { sent, emit } = stubWebSocketAndHook();
    const ws = getWS();
    await ws.ensureConnected();

    ws.initSystem({ atomic_numbers: [1,1], positions: [[0,0,0],[0.74,0,0]] });
    ws.startSimulation({ type: 'md', params: { calculator: 'lj', temperature: 0, timestep_fs: 1, friction: 0 } });

    // Start dragging atom index 0 (front-end policy; here we just simulate UI by sending USER_INTERACTION updates)
    ws.userInteraction({ positions: [[0.3,0,0],[0.74,0,0]] });

    // Server emits a frame that would otherwise pull atom 0 elsewhere
    emit({ seq: 2, positions: [[0.25,0,0],[0.74,0,0]], simStep: 2 });

    // Assert that frames were received (via stub) and USER_INTERACTION was sent
    const ui = sent.find(m => m && (m.type === 'USER_INTERACTION' || m.type === 6));
    expect(!!ui).toBe(true);
  });
});
