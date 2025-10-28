import { jest } from '@jest/globals';

export function createWsStub(overrides = {}) {
  const listeners = new Set();
  const stub = {
    connect: jest.fn(async () => {}),
    ensureConnected: jest.fn(async () => true),
    setCounters: jest.fn(),
    userInteraction: jest.fn(() => 0),
    waitForClientSeq: jest.fn(async () => {}),
    requestSingleStep: jest.fn(async () => ({
      positions: [],
      forces: [],
      energy: -1,
    })),
    waitForEnergy: jest.fn(async () => ({
      energy: -1,
      forces: [],
    })),
    onResult: jest.fn((fn) => {
      listeners.add(fn);
      return () => listeners.delete(fn);
    }),
    onFrame: jest.fn(() => () => {}),
    startSimulation: jest.fn(),
    stopSimulation: jest.fn(),
    setAck: jest.fn(),
    ack: jest.fn(),
    setTestHook: jest.fn(),
    emit(result) {
      listeners.forEach((fn) => fn(result));
    },
    ...overrides,
  };
  return stub;
}

