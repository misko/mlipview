// Centralized application constants
// Default MD friction (dimensionless UMA parameter)
export const DEFAULT_MD_FRICTION = 0.5;
// Default minimum step interval for MD/relax loops (ms)
export const DEFAULT_MIN_STEP_INTERVAL_MS = 30;

// Maximum atomic position change (in simulation length units) allowed per relax step
// when requesting relaxation from the backend. This is injected into all /relax
// API payloads so it can be centrally tuned. Smaller values -> more conservative
// line search / trust region like behavior (if backend honors the hint).
export const MAX_STEP = 0.01;
