/**
 * Centralised application constants used across the viewer runtime.
 * Keeping these as typed exports helps TypeScript infer literal types
 * without sprinkling `as const` in consuming modules yet.
 */

export const DEFAULT_MD_FRICTION: number = 0.5;
export const DEFAULT_MIN_STEP_INTERVAL_MS: number = 30;
export const MAX_STEP: number = 0.01;
