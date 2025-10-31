import { computeBaseBonds } from './baseBonds.js';
import { MIN_BOND_DISTANCE, MAX_BOND_DISTANCE } from './constants.js';

function safeMatrixInverse(a, b, c) {
  const M = [
    a.x,
    b.x,
    c.x,
    a.y,
    b.y,
    c.y,
    a.z,
    b.z,
    c.z,
  ];
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

function wrapFractional(f) {
  const { u, v, w } = f;
  return {
    u: u - Math.floor(u),
    v: v - Math.floor(v),
    w: w - Math.floor(w),
  };
}

function cartesianFromFractional({ u, v, w }, cell, origin) {
  return {
    x: origin.x + cell.a.x * u + cell.b.x * v + cell.c.x * w,
    y: origin.y + cell.a.y * u + cell.b.y * v + cell.c.y * w,
    z: origin.z + cell.a.z * u + cell.b.z * v + cell.c.z * w,
  };
}

function fractionalFromCartesian(p, origin, inv) {
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

function longThreshold(cell, options = {}) {
  try {
    if (!cell?.enabled) return 6;
    const vec = {
      x: (cell.a?.x || 0) + (cell.b?.x || 0) + (cell.c?.x || 0),
      y: (cell.a?.y || 0) + (cell.b?.y || 0) + (cell.c?.y || 0),
      z: (cell.a?.z || 0) + (cell.b?.z || 0) + (cell.c?.z || 0),
    };
    const diag = Math.hypot(vec.x, vec.y, vec.z) || 0;
    const mul = Number.isFinite(options.debugMultiplier) ? options.debugMultiplier : 0.5;
    return diag * mul || 6;
  } catch {
    return 6;
  }
}

function shiftVector(shift, cell) {
  const sx = shift[0] || 0;
  const sy = shift[1] || 0;
  const sz = shift[2] || 0;
  return {
    x: sx * cell.a.x + sy * cell.b.x + sz * cell.c.x,
    y: sx * cell.a.y + sy * cell.b.y + sz * cell.c.y,
    z: sx * cell.a.z + sy * cell.b.z + sz * cell.c.z,
  };
}

export function augmentWithPeriodicImages({ atoms, cell, options = {} }) {
  const origin = cell.originOffset || { x: 0, y: 0, z: 0 };
  const inv = safeMatrixInverse(cell.a, cell.b, cell.c);
  const debug = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;
  const shifts = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
  ];

  const augmentedAtoms = [];
  const metadata = [];
  for (let si = 0; si < shifts.length; si++) {
    const shift = shifts[si];
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
      metadata.push({ baseIndex: i, shiftIndex: si, shift });
    }
  }

  const augmented = computeBaseBonds(
    augmentedAtoms.map((a) => ({
      element: a.element,
      pos: [a.pos.x, a.pos.y, a.pos.z],
    })),
    options
  );

  const primaryBondMap = new Map();
  const crossingBondMap = new Map();
  const ghostBondMeta = [];
  const ghostMetaSeen = new Set();
  const LONG_THR = longThreshold(cell, options);
  const ghostShiftSet = new Set();
  const DEBUG = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;

  function isZeroOffset(offset) {
    return !offset || (offset[0] === 0 && offset[1] === 0 && offset[2] === 0);
  }

  function shiftsEqual(a, b) {
    return a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
  }

  function decomposeShiftPath(from, to) {
    if (shiftsEqual(from, to)) return [];
    const segments = [];
    const current = from.slice();
    const target = to.slice();
    while (!shiftsEqual(current, target)) {
      const next = current.slice();
      for (let axis = 0; axis < 3; axis++) {
        if (current[axis] < target[axis]) next[axis] += 1;
        else if (current[axis] > target[axis]) next[axis] -= 1;
      }
      segments.push({ from: current.slice(), to: next.slice() });
      current[0] = next[0];
      current[1] = next[1];
      current[2] = next[2];
    }
    return segments;
  }

  for (const bond of augmented.bonds) {
    const metaA = metadata[bond.i];
    const metaB = metadata[bond.j];
    if (!metaA || !metaB) continue;
    if (metaA.baseIndex === metaB.baseIndex) continue;

    const baseI = Math.min(metaA.baseIndex, metaB.baseIndex);
    const baseJ = Math.max(metaA.baseIndex, metaB.baseIndex);
    const shiftVecA = shiftVector(metaA.shift, cell);
    const shiftVecB = shiftVector(metaB.shift, cell);
    const posA =
      metadata[bond.i].shiftIndex === 0
        ? atoms[metaA.baseIndex].pos
        : [
            atoms[metaA.baseIndex].pos[0] + shiftVecA.x,
            atoms[metaA.baseIndex].pos[1] + shiftVecA.y,
            atoms[metaA.baseIndex].pos[2] + shiftVecA.z,
          ];
    const posB =
      metadata[bond.j].shiftIndex === 0
        ? atoms[metaB.baseIndex].pos
        : [
            atoms[metaB.baseIndex].pos[0] + shiftVecB.x,
            atoms[metaB.baseIndex].pos[1] + shiftVecB.y,
            atoms[metaB.baseIndex].pos[2] + shiftVecB.z,
          ];

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

    const minVec = {
      x: cell.a.x * du + cell.b.x * dv + cell.c.x * dw,
      y: cell.a.y * du + cell.b.y * dv + cell.c.y * dw,
      z: cell.a.z * du + cell.b.z * dv + cell.c.z * dw,
    };
    const minDist = Math.sqrt(minVec.x * minVec.x + minVec.y * minVec.y + minVec.z * minVec.z);
    if (minDist < 1e-4 || minDist > MAX_BOND_DISTANCE) {
      if (DEBUG) {
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

    const normalizeShift = (arr) =>
      Array.isArray(arr)
        ? arr.slice(0, 3).map((v) => {
            const num = Number(v);
            return Number.isFinite(num) ? Math.round(num) : 0;
          })
        : [0, 0, 0];

    let shiftA = normalizeShift(metaA.shift);
    let shiftB = normalizeShift(metaB.shift);

    const compareShift = (a, b) => {
      for (let idx = 0; idx < 3; idx++) {
        if (a[idx] === b[idx]) continue;
        return a[idx] < b[idx] ? -1 : 1;
      }
      return 0;
    };

    if (idxB < idxA) {
      // Canonical orientation: smaller atom index first.
      const tmpIdx = idxA;
      idxA = idxB;
      idxB = tmpIdx;
      const tmpShift = shiftA;
      shiftA = shiftB;
      shiftB = tmpShift;
    } else if (idxA === idxB && compareShift(shiftB, shiftA) < 0) {
      const tmpShift = shiftA;
      shiftA = shiftB;
      shiftB = tmpShift;
    }

    const isPrimary = isZeroOffset(shiftA) && isZeroOffset(shiftB);
    const crossing = crossImage || wrapCross || !isPrimary;

    if (DEBUG) {
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
      const axisOffsets = [
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
      ];
      for (const offset of axisOffsets) {
        const ghostKey = `${idxA}_${idxB}_${offset.join(',')}`;
        if (!ghostMetaSeen.has(ghostKey)) {
          ghostMetaSeen.add(ghostKey);
          ghostBondMeta.push({
            base: { i: idxA, j: idxB },
            shiftA: offset.slice(),
            shiftB: offset.slice(),
            imageDelta: [0, 0, 0],
            length: minDist,
            weight: bond.weight,
            opacity: 0,
            crossing: true,
          });
        }
        if (!isZeroOffset(offset)) ghostShiftSet.add(offset.join(','));
      }

      const mapKey = `${idxA}_${idxB}`;
      const existing = primaryBondMap.get(mapKey);
      if (!existing || existing.opacity < bond.opacity) {
        primaryBondMap.set(mapKey, {
          i: idxA,
          j: idxB,
          length: minDist,
          weight: bond.weight,
          opacity: bond.opacity,
          inRing: bond.inRing,
          crossing: false,
          imageDelta: [0, 0, 0],
          cellOffsetA: shiftA.slice(),
          cellOffsetB: shiftB.slice(),
        });
      }
      continue;
    }

    const pushGhostOffset = (offset) => {
      if (!offset || isZeroOffset(offset)) return;
      const key = `${idxA}_${idxB}_${offset.join(',')}`;
      if (ghostMetaSeen.has(key)) return;
      ghostMetaSeen.add(key);
      ghostBondMeta.push({
        base: { i: idxA, j: idxB },
        shiftA: offset.slice(),
        shiftB: offset.slice(),
        imageDelta: [0, 0, 0],
        length: minDist,
        weight: bond.weight,
        opacity: 0,
        crossing: true,
      });
      ghostShiftSet.add(offset.join(','));
    };

    pushGhostOffset(shiftA);
    pushGhostOffset(shiftB);

    if (shiftsEqual(shiftA, shiftB)) {
      // Both endpoints live in the same translated cell; this duplicates a primary bond, skip.
      continue;
    }

    const delta = [
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
        weight: bond.weight,
        opacity: 0,
        inRing: bond.inRing,
        crossing: true,
        imageDelta: delta,
        cellOffsetA: shiftA.slice(),
        cellOffsetB: shiftB.slice(),
        },
        deltaMagSq,
        shiftScore,
      });
    }
  }

  const ghostAtomList = [];
  for (const key of ghostShiftSet) {
    const shift = key.split(',').map((v) => Number(v) || 0);
    const vec = shiftVector(shift, cell);
    for (let i = 0; i < atoms.length; i++) {
      const base = atoms[i];
      const position = [
        base.pos[0] + vec.x,
        base.pos[1] + vec.y,
        base.pos[2] + vec.z,
      ];
      ghostAtomList.push({ atomIndex: i, shift: shift.slice(), position });
    }
  }

  return {
    bonds: [
      ...primaryBondMap.values(),
      ...Array.from(crossingBondMap.values()).map((entry) => entry.bond),
    ],
    ghostAtoms: ghostAtomList,
    ghostBondMeta,
    diagnostics: {
      augmentedAtoms: augmentedAtoms.length,
      candidateBonds: augmented.bonds.length,
      acceptedBonds: primaryBondMap.size,
      ghostShifts: ghostShiftSet.size,
    },
  };
}
