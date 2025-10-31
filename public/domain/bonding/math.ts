import type { Vec3 } from './types.ts';

export function vsub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]] as Vec3;
}

export function vadd(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]] as Vec3;
}

export function vscale(v: Vec3, s: number): Vec3 {
  return [v[0] * s, v[1] * s, v[2] * s] as Vec3;
}

export function vlen(v: Vec3): number {
  return Math.hypot(v[0], v[1], v[2]);
}

export function vcross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ] as Vec3;
}

export function vnorm(v: Vec3): Vec3 {
  const length = vlen(v) || 1;
  return vscale(v, 1 / length);
}

export function keyForPair(i: number, j: number): string {
  return i < j ? `${i}_${j}` : `${j}_${i}`;
}

export type { Vec3 };
