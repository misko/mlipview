// Unified manipulation service for atom dragging & bond rotation.
// Works headless; caller supplies ray -> plane intersection utilities from rendering layer.
// API:
//   beginDrag(pointerRay) if atom selected
//   updateDrag(pointerRay)
//   endDrag()
//   rotateBond(deltaAngleRadians) (uses current bond orientation side from selection)
// Emits position change via molState.markPositionsChanged(); does not recompute bonds.

import { orientationToSide } from '../selection-model.js';
import { computeBondRotationGroup } from './bondRotationUtils.js';
import { __count } from '../util/funcCount.ts';
import type {
  SelectionState,
  Selection,
  BondSelection,
  AtomSelection,
  BondRotationGroup,
} from './selectionService.ts';

type Vector3 = {
  x: number;
  y: number;
  z: number;
};

interface BondEntry {
  i: number;
  j: number;
  opacity?: number;
}

interface ManipulationState extends SelectionState {
  positions: Vector3[];
  bonds: BondEntry[];
  markPositionsChanged: () => void;
}

interface BondService {
  recomputeAndStore?: () => void;
}

type Intersector = (planePoint: Vector3, planeNormal: Vector3) => Vector3 | null;

type DragSource = 'vr' | 'desktop' | 'unknown';

type DragState = {
  atomIndex: number;
  startPos: Vector3;
  grabOffset: Vector3;
  planeNormal: Vector3;
  planePoint: Vector3;
  source: DragSource;
};

type RotateSide = 'i' | 'j';

type Options = {
  planePoint?: Vector3;
  planeNormal?: Vector3;
  source?: DragSource;
};

type ManipulationDebug = {
  getDragState: () => (DragState | null);
  getLastRotation: () => (BondRotationGroup | null);
};

export interface ManipulationService {
  beginDrag: (intersector: Intersector, opts?: Options) => boolean;
  updateDrag: (intersector: Intersector) => boolean;
  endDrag: () => void;
  setDragPlane: (point: Vector3, normal: Vector3) => void;
  rotateBond: (deltaAngle: number, sideOverride?: RotateSide) => boolean;
  setInteractionsEnabled: (on?: boolean) => boolean;
  _debug: ManipulationDebug;
}

export function createManipulationService(
  molState: ManipulationState,
  { bondService }: { bondService?: BondService } = {},
): ManipulationService {
  __count('manipulationService#createManipulationService');
  let dragState: DragState | null = null;
  const MAX_DRAG_RADIUS = 50;
  const DRAG_LOG = false;
  let interactionsEnabled = true;
  let lastRotationDebug: BondRotationGroup | null = null;

  const getAtomPosition = (i: number): Vector3 => molState.positions[i];
  const setAtomPosition = (i: number, p: Vector3): void => {
    molState.positions[i].x = p.x;
    molState.positions[i].y = p.y;
    molState.positions[i].z = p.z;
  };

  function beginDrag(intersector: Intersector, opts: Options = {}): boolean {
    __count('manipulationService#beginDrag');
    if (!interactionsEnabled) return false;
    if (dragState) return false;
    const selection = molState.selection as Selection | undefined;
    if (!selection || selection.kind !== 'atom') return false;

    const atomIndex = (selection as AtomSelection).data.index;
    const pos = getAtomPosition(atomIndex);

    const planePoint = opts.planePoint ? { ...opts.planePoint } : { ...pos };
    const planeNormal = opts.planeNormal ? { ...opts.planeNormal } : { x: 0, y: 1, z: 0 };

    const hit = intersector(planePoint, planeNormal);
    let grabOffset: Vector3 = hit
      ? { x: pos.x - hit.x, y: pos.y - hit.y, z: pos.z - hit.z }
      : { x: 0, y: 0, z: 0 };
    if (Math.hypot(grabOffset.x, grabOffset.y, grabOffset.z) < 1e-6) {
      grabOffset = { x: 0, y: 0, z: 0 };
    }

    const source: DragSource = opts.source || 'unknown';
    dragState = { atomIndex, startPos: { ...pos }, grabOffset, planeNormal, planePoint, source };
    if (DRAG_LOG) {
      console.log('[drag][start]', { atomIndex, startPos: { ...pos }, planePoint, planeNormal });
    }
    return true;
  }

  function setDragPlane(point: Vector3, normal: Vector3): void {
    __count('manipulationService#setDragPlane');
    if (!interactionsEnabled || !dragState) return;
    dragState.planePoint = { ...point };
    dragState.planeNormal = { ...normal };
  }

  function updateDrag(intersector: Intersector): boolean {
    __count('manipulationService#updateDrag');
    if (!interactionsEnabled || !dragState) return false;
    const hit = intersector(dragState.planePoint, dragState.planeNormal);
    if (!hit) return false;

    const prev = { ...molState.positions[dragState.atomIndex] };
    let newPos: Vector3 = {
      x: hit.x + dragState.grabOffset.x,
      y: hit.y + dragState.grabOffset.y,
      z: hit.z + dragState.grabOffset.z,
    };

    const r = Math.hypot(newPos.x, newPos.y, newPos.z);
    if (r > MAX_DRAG_RADIUS) {
      const s = MAX_DRAG_RADIUS / r;
      newPos = { x: newPos.x * s, y: newPos.y * s, z: newPos.z * s };
    }
    setAtomPosition(dragState.atomIndex, newPos);
    molState.markPositionsChanged();

    if (DRAG_LOG) {
      const dx = +(newPos.x - prev.x).toFixed(3);
      const dy = +(newPos.y - prev.y).toFixed(3);
      const dz = +(newPos.z - prev.z).toFixed(3);
      console.log('[drag][move]', {
        atomIndex: dragState.atomIndex,
        from: prev,
        to: { ...newPos },
        d: { dx, dy, dz },
      });
    }
    return true;
  }

  function endDrag(): void {
    __count('manipulationService#endDrag');
    if (!interactionsEnabled || !dragState) return;
    const { atomIndex, startPos } = dragState;
    const cur = molState.positions[atomIndex];
    const moved =
      Math.abs(cur.x - startPos.x) > 1e-9 ||
      Math.abs(cur.y - startPos.y) > 1e-9 ||
      Math.abs(cur.z - startPos.z) > 1e-9;

    if (DRAG_LOG) {
      const dx = +(cur.x - startPos.x).toFixed(3);
      const dy = +(cur.y - startPos.y).toFixed(3);
      const dz = +(cur.z - startPos.z).toFixed(3);
      console.log('[drag][drop]', { atomIndex, moved, delta: { dx, dy, dz }, startPos, finalPos: { ...cur } });
    }

    dragState = null;
    if (moved && bondService?.recomputeAndStore) {
      bondService.recomputeAndStore();
    }
  }

  function rotateBond(deltaAngle: number, sideOverride?: RotateSide): boolean {
    __count('manipulationService#rotateBond');
    if (!interactionsEnabled) return false;
    const selection = molState.selection as Selection | undefined;
    if (!selection || selection.kind !== 'bond') return false;

    const { i, j, orientation } = (selection as BondSelection).data;
    const side: RotateSide = sideOverride || orientationToSide(orientation) || 'j';
    let group = (selection as BondSelection).data.rotationGroup;
    if (!group || group.orientation !== orientation || group.i !== i || group.j !== j) {
      group = computeBondRotationGroup(molState as unknown as Parameters<typeof computeBondRotationGroup>[0], {
        i,
        j,
        orientation,
      });
      (selection as BondSelection).data.rotationGroup = group;
    }

    const { anchor, movingRoot, sideAtoms } = group;
    const pAnchor = molState.positions[anchor];
    const pMove = molState.positions[movingRoot];
    const ax = pMove.x - pAnchor.x;
    const ay = pMove.y - pAnchor.y;
    const az = pMove.z - pAnchor.z;
    const len = Math.hypot(ax, ay, az) || 1;
    const ux = ax / len;
    const uy = ay / len;
    const uz = az / len;
    const s = Math.sin(deltaAngle);
    const c = Math.cos(deltaAngle);

    const rotatePoint = (p: Vector3): void => {
      const dx = p.x - pAnchor.x;
      const dy = p.y - pAnchor.y;
      const dz = p.z - pAnchor.z;
      const dot = dx * ux + dy * uy + dz * uz;
      const rx = ux * dot + (dx - ux * dot) * c + (-uz * dy + uy * dz) * s;
      const ry = uy * dot + (dy - uy * dot) * c + (uz * dx - ux * dz) * s;
      const rz = uz * dot + (dz - uz * dot) * c + (-uy * dx + ux * dy) * s;
      p.x = pAnchor.x + rx;
      p.y = pAnchor.y + ry;
      p.z = pAnchor.z + rz;
    };

    for (const idx of sideAtoms) {
      rotatePoint(molState.positions[idx]);
    }

    lastRotationDebug = {
      ...group,
      side,
      sideAtoms: [...sideAtoms],
    };

    molState.markPositionsChanged();
    if (bondService?.recomputeAndStore) {
      bondService.recomputeAndStore();
    }
    return true;
  }

  const _debug: ManipulationDebug = {
    getDragState: () =>
      dragState
        ? {
            ...dragState,
            planeNormal: { ...dragState.planeNormal },
            planePoint: { ...dragState.planePoint },
            grabOffset: { ...dragState.grabOffset },
          }
        : null,
    getLastRotation: () => (lastRotationDebug ? { ...lastRotationDebug } : null),
  };

  return {
    beginDrag,
    updateDrag,
    endDrag,
    setDragPlane,
    rotateBond,
    setInteractionsEnabled(on: boolean = true): boolean {
      interactionsEnabled = !!on;
      if (!interactionsEnabled) {
        try {
          endDrag();
        } catch {}
      }
      return interactionsEnabled;
    },
    _debug,
  };
}
