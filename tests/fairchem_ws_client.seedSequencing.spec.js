/** @jest-environment node */

import { createFairchemWS } from '../public/fairchem_ws_client.js';

describe('fairchem_ws_client seedSequencing', () => {
  const originalGlobalWs = globalThis.__fairchem_ws__;
  let api;

  afterEach(() => {
    try { api?.close?.(); } catch { }
    api = null;
    if (originalGlobalWs === undefined) {
      delete globalThis.__fairchem_ws__;
    } else {
      globalThis.__fairchem_ws__ = originalGlobalWs;
    }
  });

  it('honours stored next sequence when replaying snapshots', () => {
    jest.useFakeTimers();
    try {
      api = createFairchemWS();
      api.seedSequencing({
        nextSeq: 21,
        lastSeq: 20,
        ack: 18,
        userInteractionCount: 3,
        simStep: 2,
      });

      expect(api.getState().seq).toBe(20);
      expect(api.getState().clientAck).toBe(18);

      const firstSeq = api.userInteraction({
        natoms: 1,
        positions: [[0, 0, 0]],
        full_update: true,
      });
      expect(firstSeq).toBe(21);

      const secondSeq = api.userInteraction({
        natoms: 1,
        positions: [[0, 0, 0]],
        full_update: true,
      });
      expect(secondSeq).toBe(22);

      jest.runOnlyPendingTimers();
    } finally {
      jest.useRealTimers();
    }
  });
});
