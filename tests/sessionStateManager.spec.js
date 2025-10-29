import { createSessionStateManager } from '../public/core/sessionStateManager.js';

describe('SessionStateManager', () => {
  it('rehydrates viewer and seeds websocket state', async () => {
    const Mode = { Idle: 'idle', MD: 'md', Relax: 'relax', Timeline: 'timeline' };
    let currentMode = Mode.Idle;
    const state = {
      elements: ['H', 'O'],
      positions: [{ x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }],
      showCell: false,
      showGhostCells: false,
      dynamics: {},
      markCellChanged: jest.fn(),
    };
    const applyFullSnapshot = jest.fn();
    const energyPlot = {
      exportSeries: () => ({ series: [{ energy: -10, kind: 'idle' }], markerIndex: 0 }),
      importSeries: jest.fn(),
      setMarker: jest.fn(),
    };
    const frameStore = [];
    const frameBuffer = {
      exportFrames: () => frameStore.slice(),
      importFrames: jest.fn((frames) => {
        frameStore.splice(0, frameStore.length, ...(frames || []));
      }),
      size: () => frameStore.length,
      getByOffset: (offset) => {
        if (!frameStore.length || offset !== -1) return null;
        return {
          ...frameStore[frameStore.length - 1],
          positions: frameStore[frameStore.length - 1].positions,
        };
      },
    };
    const installBaseline = jest.fn();
    const setInteractionCounters = jest.fn();
    const resetWsInitState = jest.fn();
    const seedSequencing = jest.fn();
    const userInteraction = jest.fn();
    const ensureWsInit = jest.fn(() => Promise.resolve());
    const applyTimelineFrame = jest.fn();
    const rememberResume = jest.fn();

    const manager = createSessionStateManager({
      getViewerState: () => state,
      applyFullSnapshot,
      normalizeElement: (el) => el,
      energyPlot,
      frameBuffer,
      captureBaselineFromState: () => ({ elements: state.elements, positions: state.positions }),
      installBaseline,
      getInteractionCounters: () => ({ user: 0, total: 0, lastApplied: 0 }),
      setInteractionCounters,
      resetWsInitState,
      getWsClient: () => ({ seedSequencing, userInteraction }),
      ensureWsInit,
      posToTriples: (s) => s.positions.map((p) => [p.x, p.y, p.z]),
      zOf: (el) => (el === 'H' ? 1 : 8),
      stateCellToArray: () => null,
      getVelocitiesForSnapshot: () => null,
      setMode: (next) => { currentMode = next; },
      Mode,
      timelineState: { active: false, playing: false, suppressEnergy: false, offset: -1, resumeMode: null },
      timelineUiRef: { value: null },
      clearTimelinePlayback: jest.fn(),
      setTimelineOverlayVisible: jest.fn(),
      setTimelineInteractionLock: jest.fn(),
      applyTimelineFrame,
      refreshTimelineUi: jest.fn(),
      setTimelineUiMode: jest.fn(),
      lastContinuousOpts: { md: null, relax: null },
      rememberResume,
      stopSimulation: jest.fn(),
      getMode: () => currentMode,
    });

    const snapshot = manager.captureSnapshot({ kind: 'xyz', label: 'test.xyz' });
    expect(snapshot.viewer.elements).toEqual(['H', 'O']);
    expect(snapshot.energyPlot.series).toHaveLength(1);

    const loadSnapshot = {
      schemaVersion: 1,
      savedAt: new Date().toISOString(),
      source: { kind: 'json', label: 'session.json' },
      viewer: {
        elements: ['H', 'O'],
        positions: [[0, 0, 0], [1.1, 0.1, 0]],
        showCell: false,
        showGhostCells: false,
      },
      energyPlot: { series: [{ energy: -9.8, kind: 'md' }], markerIndex: 0 },
      timeline: {
        capacity: 500,
        frames: [{
          kind: 'md',
          seq: 15,
          simStep: 2,
          userInteractionCount: 3,
          energy: -9.8,
          positions: [[0, 0, 0], [1.1, 0.1, 0]],
        }],
        lastLiveMode: 'md',
        wasRunning: true,
        pendingSimParams: { temperature: 300 },
      },
      websocket: { seq: 20, clientAck: 18, userInteractionCount: 3, simStep: 2 },
    };

    await manager.loadSnapshot(loadSnapshot);

    expect(applyFullSnapshot).toHaveBeenCalled();
    expect(energyPlot.importSeries).toHaveBeenCalledWith(loadSnapshot.energyPlot);
    expect(frameBuffer.importFrames).toHaveBeenCalledWith(loadSnapshot.timeline.frames);
    expect(seedSequencing).toHaveBeenCalledWith(expect.objectContaining({ nextSeq: 20 }));
    expect(userInteraction).toHaveBeenCalledWith(expect.objectContaining({ full_update: true }));
    expect(ensureWsInit).toHaveBeenCalled();
    expect(currentMode).toBe(Mode.Timeline);
    expect(rememberResume).toHaveBeenCalledWith('md', expect.any(Object));
    expect(setInteractionCounters).toHaveBeenCalledWith(expect.objectContaining({ user: 3 }));
  });
});
