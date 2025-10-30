export function vsub(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function vadd(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function vscale(v, s) {
  return [v[0] * s, v[1] * s, v[2] * s];
}

export function vlen(v) {
  return Math.hypot(v[0], v[1], v[2]);
}

export function vcross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function vnorm(v) {
  const L = vlen(v) || 1;
  return vscale(v, 1 / L);
}

export function keyForPair(i, j) {
  return i < j ? `${i}_${j}` : `${j}_${i}`;
}
