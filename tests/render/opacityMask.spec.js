import { combineOpacity, applyMaterialTransparency } from '../../public/render/moleculeView.js';

const MATERIAL_CONSTANTS = {
  MATERIAL_ALPHABLEND: 2,
  MATERIAL_OPAQUE: 0,
};

beforeEach(() => {
  global.BABYLON = {
    Material: MATERIAL_CONSTANTS,
  };
});

afterEach(() => {
  delete global.BABYLON;
});

describe('combineOpacity', () => {
  it('multiplies base and mask when using multiply mode', () => {
    expect(combineOpacity(0.8, 0.5, 'multiply')).toBeCloseTo(0.4);
  });

  it('overrides base when using override mode', () => {
    expect(combineOpacity(0.2, 0.6, 'override')).toBeCloseTo(0.6);
  });

  it('clamps to the minimum when using clamp mode', () => {
    expect(combineOpacity(0.7, 0.3, 'clamp')).toBeCloseTo(0.3);
  });

  it('normalises invalid blend modes to multiply', () => {
    expect(combineOpacity(0.5, 0.5, 'unknown')).toBeCloseTo(0.25);
  });

  it('clamps inputs outside the 0-1 range', () => {
    expect(combineOpacity(1.6, 2.0, 'multiply')).toBeCloseTo(1);
    expect(combineOpacity(-1, -5, 'override')).toBeCloseTo(0);
  });
});

describe('applyMaterialTransparency', () => {
  it('enables depth pre-pass without depth writes when transparent and requested', () => {
    const mat = {};
    applyMaterialTransparency(mat, true, { usePrePass: true });
    expect(mat.transparencyMode).toBe(MATERIAL_CONSTANTS.MATERIAL_ALPHABLEND);
    expect(mat.forceDepthWrite).toBe(false);
    expect(mat.needDepthPrePass).toBe(true);
    expect(mat.separateCullingPass).toBe(true);
  });

  it('disables pre-pass and restores opaque defaults', () => {
    const mat = {
      transparencyMode: MATERIAL_CONSTANTS.MATERIAL_ALPHABLEND,
      forceDepthWrite: true,
      needDepthPrePass: true,
      separateCullingPass: true,
    };
    applyMaterialTransparency(mat, false);
    expect(mat.transparencyMode).toBe(MATERIAL_CONSTANTS.MATERIAL_OPAQUE);
    expect(mat.forceDepthWrite).toBe(false);
    expect(mat.needDepthPrePass).toBe(false);
    expect(mat.separateCullingPass).toBe(false);
  });

  it('leaves pre-pass disabled if not requested', () => {
    const mat = {};
    applyMaterialTransparency(mat, true, { usePrePass: false });
    expect(mat.needDepthPrePass).toBe(false);
    expect(mat.forceDepthWrite).toBe(false);
  });
});
