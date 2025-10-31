export type Vec3 = [number, number, number];

export interface Cartesian {
  x: number;
  y: number;
  z: number;
}

export interface Fractional {
  u: number;
  v: number;
  w: number;
}

export interface NormalizedAtom {
  element: string;
  pos: Vec3;
}

export interface PeriodicCell {
  enabled: boolean;
  a: Cartesian;
  b: Cartesian;
  c: Cartesian;
  originOffset?: Cartesian;
}

export interface BaseBondEntry {
  i: number;
  j: number;
  length: number | null;
  weight: number | null;
  opacity: number;
  inRing?: boolean;
}

export interface BondRecord extends BaseBondEntry {
  crossing: boolean;
  imageDelta: Vec3;
  cellOffsetA: Vec3;
  cellOffsetB: Vec3;
}

export interface GhostAtom {
  atomIndex: number;
  shift: Vec3;
  position: Vec3;
}

export interface GhostBondMeta {
  base: { i: number; j: number };
  shiftA: Vec3;
  shiftB: Vec3;
  imageDelta: Vec3;
  length?: number | null;
  weight?: number | null;
  opacity?: number;
  crossing?: boolean;
}

export interface BondDiagnostics {
  [key: string]: unknown;
}

export interface BaseBondResult {
  bonds: BaseBondEntry[];
  diagnostics?: BondDiagnostics;
}

export interface PeriodicAugmentResult {
  bonds: BondRecord[];
  ghostAtoms: GhostAtom[];
  ghostBondMeta: GhostBondMeta[];
  diagnostics?: BondDiagnostics;
}

export interface ComputeBondsResult {
  bonds: BondRecord[];
  ghostAtoms: GhostAtom[];
  ghostBondMeta: GhostBondMeta[];
  diagnostics: {
    base?: BondDiagnostics;
    periodic?: BondDiagnostics;
  };
}

export type BondComputationOptions = {
  maxNeighbors?: number;
  debugMultiplier?: number;
  [key: string]: unknown;
};
