import { computeBaseBonds } from './baseBonds.ts';
import { MAX_BOND_DISTANCE } from './constants.ts';
import type {
  BaseBondResult,
  BondComputationOptions,
  BondRecord,
  GhostAtom,
  GhostBondMeta,
  NormalizedAtom,
  PeriodicAugmentResult,
  PeriodicCell,
  Vec3,
} from './types.ts';

interface Cartesian {
  x: number;
  y: number;
  z: number;
}

type Matrix3 = number[];

interface AugmentedAtom {
  element: string;
  pos: Cartesian;
}

interface AtomMeta {
  baseIndex: number;
  shiftIndex: number;
  shift: Vec3;
}

const AXIS_SHIFTS: readonly Vec3[] = [
  [0, 0, 0],
  [1, 0, 0],
  [-1, 0, 0],
  [0, 1, 0],
  [0, -1, 0],
  [0, 0, 1],
  [0, 0, -1],
];

const PRIMARY_GHOST_OFFSETS: readonly Vec3[] = [
  [1, 0, 0],
  [-1, 0, 0],
  [0, 1, 0],
  [0, -1, 0],
  [0, 0, 1],
  [0, 0, -1],
];

function safeMatrixInverse(a: Cartesian, b: Cartesian, c: Cartesian): Matrix3 | null {
  const M: Matrix3 = [a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z];
  const det =
    M[0] * (M[4] * M[8] - M[5] * M[7]) -
    M[1] * (M[3] * M[8] - M[5] * M[6]) +
    M[2] * (M[3] * M[7] - M[4] * M[6]);
  if (Math.abs(det) <= 1e-14) return null;
  const invDet = 1 / det;
  return [
    (M[4] * M[8] - M[5] * M[7]) * invDet,
    (M[2] * M[7] - M[1] * M[8]) * invDet,
    (M[1] * M[5] - M[2] * M[4]) * invDet,
    (M[5] * M[6] - M[3] * M[8]) * invDet,
    (M[0] * M[8] - M[2] * M[6]) * invDet,
    (M[2] * M[3] - M[0] * M[5]) * invDet,
    (M[3] * M[7] - M[4] * M[6]) * invDet,
    (M[1] * M[6] - M[0] * M[7]) * invDet,
    (M[0] * M[4] - M[1] * M[3]) * invDet,
  ];
}

function wrapFractional(f: { u: number; v: number; w: number }) {
  const { u, v, w } = f;
  return {
    u: u - Math.floor(u),
    v: v - Math.floor(v),
    w: w - Math.floor(w),
  };
}

function cartesianFromFractional(coords: { u: number; v: number; w: number }, cell: PeriodicCell, origin: Cartesian): Cartesian {
  const { u, v, w } = coords;
  return {
    x: origin.x + cell.a.x * u + cell.b.x * v + cell.c.x * w,
    y: origin.y + cell.a.y * u + cell.b.y * v + cell.c.y * w,
    z: origin.z + cell.a.z * u + cell.b.z * v + cell.c.z * w,
  };
}

function fractionalFromCartesian(p: Cartesian, origin: Cartesian, inv: Matrix3 | null) {
  const dx = p.x - origin.x;
  const dy = p.y - origin.y;
  const dz = p.z - origin.z;
  if (!inv) {
    return { u: dx, v: dy, w: dz };
  }
  return {
    u: inv[0] * dx + inv[1] * dy + inv[2] * dz,
    v: inv[3] * dx + inv[4] * dy + inv[5] * dz,
    w: inv[6] * dx + inv[7] * dy + inv[8] * dz,
  };
}

function longThreshold(cell: PeriodicCell, options: BondComputationOptions = {}): number {
  try {
    if (!cell?.enabled) return 6;
    const vec = {
      x: (cell.a?.x || 0) + (cell.b?.x || 0) + (cell.c?.x || 0),
      y: (cell.a?.y || 0) + (cell.b?.y || 0) + (cell.c?.y || 0),
      z: (cell.a?.z || 0) + (cell.b?.z || 0) + (cell.c?.z || 0),
    };
    const diag = Math.hypot(vec.x, vec.y, vec.z) || 0;
    const mul = Number.isFinite(options.debugMultiplier) ? (options.debugMultiplier as number) : 0.5;
    return diag * mul || 6;
  } catch {
    return 6;
  }
}

function shiftVector(shift: Vec3, cell: PeriodicCell): Cartesian {
  const sx = shift[0] || 0;
  const sy = shift[1] || 0;
  const sz = shift[2] || 0;
  return {
    x: sx * cell.a.x + sy * cell.b.x + sz * cell.c.x,
    y: sx * cell.a.y + sy * cell.b.y + sz * cell.c.y,
    z: sx * cell.a.z + sy * cell.b.z + sz * cell.c.z,
  };
}

function isZeroOffset(offset: Vec3 | number[] | undefined | null): boolean {
  return !offset || (offset[0] === 0 && offset[1] === 0 && offset[2] === 0);
}

function shiftsEqual(a: Vec3, b: Vec3): boolean {
  return a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
}

function normalizeShift(arr: unknown): Vec3 {
  if (!Array.isArray(arr)) return [0, 0, 0];
  const out: Vec3 = [
    Number.isFinite(arr[0]) ? Math.round(arr[0] as number) : 0,
    Number.isFinite(arr[1]) ? Math.round(arr[1] as number) : 0,
    Number.isFinite(arr[2]) ? Math.round(arr[2] as number) : 0,
  ];
  return out;
}

function compareShift(a: Vec3, b: Vec3): number {
  for (let idx = 0; idx < 3; idx++) {
    if (a[idx] === b[idx]) continue;
    return a[idx] < b[idx] ? -1 : 1;
  }
  return 0;
}

export function augmentWithPeriodicImages({
  atoms,
  cell,
  options = {},
}: {
  atoms: NormalizedAtom[];
  cell: PeriodicCell;
  options?: BondComputationOptions;
}): PeriodicAugmentResult {
  const origin: Cartesian = cell.originOffset || { x: 0, y: 0, z: 0 };
  const inv = safeMatrixInverse(cell.a, cell.b, cell.c);
  const debug = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;

  const augmentedAtoms: AugmentedAtom[] = [];
  const metadata: AtomMeta[] = [];

  for (let si = 0; si < AXIS_SHIFTS.length; si++) {
    const shift = AXIS_SHIFTS[si];
    const vec = shiftVector(shift, cell);
    for (let i = 0; i < atoms.length; i++) {
      const base = atoms[i];
      augmentedAtoms.push({
        element: base.element,
        pos: {
          x: base.pos[0] + vec.x,
          y: base.pos[1] + vec.y,
          z: base.pos[2] + vec.z,
        },
      });
      metadata.push({ baseIndex: i, shiftIndex: si, shift: [...shift] as Vec3 });
    }
  }

  const augmented = computeBaseBonds(
    augmentedAtoms.map((a) => ({
      element: a.element,
      pos: [a.pos.x, a.pos.y, a.pos.z] as Vec3,
    })),
    options,
  ) as BaseBondResult;

  const primaryBondMap = new Map<string, BondRecord>();
  const crossingBondMap = new Map<
    string,
    { bond: BondRecord; deltaMagSq: number; shiftScore: number }
  >();
  const ghostBondMeta: GhostBondMeta[] = [];
  const ghostMetaSeen = new Set<string>();
  const LONG_THR = longThreshold(cell, options);
  const ghostShiftSet = new Set<string>();

  for (const bond of augmented.bonds) {
    const metaA = metadata[bond.i];
    const metaB = metadata[bond.j];
    if (!metaA || !metaB) continue;
    if (metaA.baseIndex === metaB.baseIndex) continue;

    const shiftVecA = shiftVector(metaA.shift, cell);
    const shiftVecB = shiftVector(metaB.shift, cell);
    const posA =
      metadata[bond.i].shiftIndex === 0
        ? atoms[metaA.baseIndex].pos
        : ([
            atoms[metaA.baseIndex].pos[0] + shiftVecA.x,
            atoms[metaA.baseIndex].pos[1] + shiftVecA.y,
            atoms[metaA.baseIndex].pos[2] + shiftVecA.z,
          ] as Vec3);
    const posB =
      metadata[bond.j].shiftIndex === 0
        ? atoms[metaB.baseIndex].pos
        : ([
            atoms[metaB.baseIndex].pos[0] + shiftVecB.x,
            atoms[metaB.baseIndex].pos[1] + shiftVecB.y,
            atoms[metaB.baseIndex].pos[2] + shiftVecB.z,
          ] as Vec3);

    const fracRawA = fractionalFromCartesian({ x: posA[0], y: posA[1], z: posA[2] }, origin, inv);
    const fracRawB = fractionalFromCartesian({ x: posB[0], y: posB[1], z: posB[2] }, origin, inv);
    const fracA = wrapFractional(fracRawA);
    const fracB = wrapFractional(fracRawB);

    const duRaw = fracB.u - fracA.u;
    const dvRaw = fracB.v - fracA.v;
    const dwRaw = fracB.w - fracA.w;
    const wrapU = Math.round(duRaw);
    const wrapV = Math.round(dvRaw);
    const wrapW = Math.round(dwRaw);
    const du = duRaw - wrapU;
    const dv = dvRaw - wrapV;
    const dw = dwRaw - wrapW;

    const minVec: Cartesian = {
      x: cell.a.x * du + cell.b.x * dv + cell.c.x * dw,
      y: cell.a.y * du + cell.b.y * dv + cell.c.y * dw,
      z: cell.a.z * du + cell.b.z * dv + cell.c.z * dw,
    };
    const minDist = Math.sqrt(minVec.x * minVec.x + minVec.y * minVec.y + minVec.z * minVec.z);
    if (minDist < 1e-4 || minDist > MAX_BOND_DISTANCE) {
      if (debug) {
        console.log('[Bonding][periodicSkipDistance]', {
          baseI: metaA.baseIndex,
          baseJ: metaB.baseIndex,
          shiftA: metaA.shift,
          shiftB: metaB.shift,
          minDist,
        });
      }
      continue;
    }

    let idxA = metaA.baseIndex;
    let idxB = metaB.baseIndex;
    if (!Number.isInteger(idxA) || !Number.isInteger(idxB) || idxA === idxB) continue;

    const crossImage = metaA.shiftIndex !== 0 || metaB.shiftIndex !== 0;
    const wrapCross = wrapU !== 0 || wrapV !== 0 || wrapW !== 0;

    let shiftA = normalizeShift(metaA.shift);
    let shiftB = normalizeShift(metaB.shift);

    if (idxB < idxA) {
      [idxA, idxB] = [idxB, idxA];
      [shiftA, shiftB] = [shiftB, shiftA];
    } else if (idxA === idxB && compareShift(shiftB, shiftA) < 0) {
      [shiftA, shiftB] = [shiftB, shiftA];
    }

    const isPrimary = isZeroOffset(shiftA) && isZeroOffset(shiftB);
    const crossing = crossImage || wrapCross || !isPrimary;

    if (debug) {
      console.log('[Bonding][periodicPair]', {
        i: idxA,
        j: idxB,
        shiftIndexA: metaA.shiftIndex,
        shiftIndexB: metaB.shiftIndex,
        shiftA: metaA.shift,
        shiftB: metaB.shift,
        shiftAResolved: shiftA,
        shiftBResolved: shiftB,
        minDist: Number(minDist.toFixed(6)),
        isPrimary,
        crossing,
      });
    }

    if (isPrimary) {
      for (const offset of PRIMARY_GHOST_OFFSETS) {
        const ghostKey = `${idxA}_${idxB}_${offset.join(',')}`;
        if (!ghostMetaSeen.has(ghostKey)) {
          ghostMetaSeen.add(ghostKey);
          ghostBondMeta.push({
            base: { i: idxA, j: idxB },
            shiftA: [...offset] as Vec3,
            shiftB: [...offset] as Vec3,
            imageDelta: [0, 0, 0],
            length: minDist,
            weight: bond.weight ?? null,
            opacity: 0,
            crossing: true,
          });
        }
        if (!isZeroOffset(offset)) ghostShiftSet.add(offset.join(','));
      }

      const mapKey = `${idxA}_${idxB}`;
      const existing = primaryBondMap.get(mapKey);
      if (!existing || (bond.opacity ?? 0) > existing.opacity) {
        primaryBondMap.set(mapKey, {
          i: idxA,
          j: idxB,
          length: minDist,
          weight: bond.weight ?? null,
          opacity: bond.opacity ?? 0,
          inRing: bond.inRing,
          crossing: false,
          imageDelta: [0, 0, 0],
          cellOffsetA: [...shiftA] as Vec3,
          cellOffsetB: [...shiftB] as Vec3,
        });
      }
      continue;
    }

    const pushGhostOffset = (offset: Vec3) => {
      if (isZeroOffset(offset)) return;
      const key = `${idxA}_${idxB}_${offset.join(',')}`;
      if (ghostMetaSeen.has(key)) return;
      ghostMetaSeen.add(key);
      ghostBondMeta.push({
        base: { i: idxA, j: idxB },
        shiftA: [...offset] as Vec3,
        shiftB: [...offset] as Vec3,
        imageDelta: [0, 0, 0],
        length: minDist,
        weight: bond.weight ?? null,
        opacity: 0,
        crossing: true,
      });
      ghostShiftSet.add(offset.join(','));
    };

    pushGhostOffset(shiftA);
    pushGhostOffset(shiftB);

    if (shiftsEqual(shiftA, shiftB)) {
      continue;
    }

    const delta: Vec3 = [
      (shiftB[0] || 0) - (shiftA[0] || 0),
      (shiftB[1] || 0) - (shiftA[1] || 0),
      (shiftB[2] || 0) - (shiftA[2] || 0),
    ];
    const deltaMagSq = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
    const shiftScore =
      Math.abs(shiftA[0] || 0) +
      Math.abs(shiftA[1] || 0) +
      Math.abs(shiftA[2] || 0) +
      Math.abs(shiftB[0] || 0) +
      Math.abs(shiftB[1] || 0) +
      Math.abs(shiftB[2] || 0);

    const crossingKey = `${idxA}_${idxB}`;
    const existingCross = crossingBondMap.get(crossingKey);
    if (
      !existingCross ||
      deltaMagSq < existingCross.deltaMagSq ||
      (deltaMagSq === existingCross.deltaMagSq && shiftScore < existingCross.shiftScore)
    ) {
      crossingBondMap.set(crossingKey, {
        bond: {
          i: idxA,
          j: idxB,
          length: minDist,
          weight: bond.weight ?? null,
          opacity: 0,
          inRing: bond.inRing,
          crossing: true,
          imageDelta: delta,
          cellOffsetA: [...shiftA] as Vec3,
          cellOffsetB: [...shiftB] as Vec3,
        },
        deltaMagSq,
        shiftScore,
      });
    }
  }

  const ghostAtoms: GhostAtom[] = [];
  for (const key of ghostShiftSet) {
    const shift = key.split(',').map((v) => Number(v) || 0) as Vec3;
    const vec = shiftVector(shift, cell);
    for (let i = 0; i < atoms.length; i++) {
      const base = atoms[i];
      const position: Vec3 = [
        base.pos[0] + vec.x,
        base.pos[1] + vec.y,
        base.pos[2] + vec.z,
      ];
      ghostAtoms.push({ atomIndex: i, shift: [...shift] as Vec3, position });
    }
  }

  const bonds: BondRecord[] = [
    ...primaryBondMap.values(),
    ...Array.from(crossingBondMap.values()).map((entry) => entry.bond),
  ];

  return {
    bonds,
    ghostAtoms,
    ghostBondMeta,
    diagnostics: {
      augmentedAtoms: augmentedAtoms.length,
      candidateBonds: augmented.bonds.length,
      acceptedBonds: primaryBondMap.size,
      ghostShifts: ghostShiftSet.size,
    },
  };
}

declare global {
  interface Window {
    __MLIPVIEW_DEBUG_BONDS?: boolean;
  }
}
