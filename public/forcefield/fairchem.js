// public/forcefield/fairchem.js
// HuggingFace FAIR-Chem remote inference adapter.
// Lightweight JSON POST wrapper modeled after hf_calculator.py behavior.
// NOTE: Endpoint & auth token must be configured by user.

import { atomsToZ, atomsToXYZ, rawForcesToVectors, validateZXYZ } from './interface.js';

// Local UMA HTTP server bridge (see fairchem_local_server). Environment / runtime can set
//   window.FAIRCHEM_SERVER = 'http://localhost:8000'
// We target its /simple_calculate endpoint.
let DEFAULT_SERVER;
let __fairchemWarnedMixed = false;
let __fairchemLockedEndpoint = null; // once we successfully use an upgraded endpoint
let __fairchemTriedHttps8444 = false; // track alternate port attempt
if (typeof window !== 'undefined') {
  if (window.FAIRCHEM_SERVER === 'AUTO_HOST_8000') {
    // explicit opt-in to legacy heuristic
    const host = window.location.hostname || 'localhost';
    DEFAULT_SERVER = `http://${host}:8000`;
    console.log('[fairchem] legacy AUTO_HOST_8000 server', DEFAULT_SERVER);
  } else if (window.FAIRCHEM_SERVER && window.FAIRCHEM_SERVER !== '' && window.FAIRCHEM_SERVER !== 'SAME_ORIGIN') {
    DEFAULT_SERVER = window.FAIRCHEM_SERVER.replace(/\/$/, '');
  } else {
    // Default to same-origin relative path
    DEFAULT_SERVER = '';
  }
} else {
  DEFAULT_SERVER = '';
}
let DEFAULT_ENDPOINT = DEFAULT_SERVER ? `${DEFAULT_SERVER}/simple_calculate` : '/simple_calculate';

// If page is HTTPS and default endpoint is HTTP, attempt auto-upgrade to HTTPS same host:port
try {
  if (typeof window !== 'undefined') {
    const loc = window.location;
    if (loc.protocol === 'https:' && DEFAULT_ENDPOINT.startsWith('http://')) {
      // If original was http://host:8000 try https://host:8444 first (script starts TLS there)
      const portMatch = DEFAULT_ENDPOINT.match(/^http:\/\/([^/]+):(\d+)(\/.*)$/);
      if (portMatch && portMatch[2] === '8000') {
        const host = portMatch[1];
        const rest = portMatch[3];
        const https8444 = `https://${host}:8444${rest}`;
        console.log('[fairchem] attempting HTTPS auto-upgrade (port 8444 first)', https8444);
        DEFAULT_ENDPOINT = https8444; // first candidate
      } else {
        const httpsCandidate = DEFAULT_ENDPOINT.replace('http://', 'https://');
        console.log('[fairchem] attempting HTTPS auto-upgrade', httpsCandidate);
        DEFAULT_ENDPOINT = httpsCandidate;
      }
    }
  }
} catch {}

function buildHeaders() { return { 'Content-Type': 'application/json' }; }

async function postJSON(endpoint, payload) {
  const t0 = performance.now();
  console.log('[fairchem] → POST', endpoint, '\n  payload:', {
    natoms: Array.isArray(payload.atomic_numbers) ? payload.atomic_numbers.length : 'n/a',
    keys: Object.keys(payload)
  });
  let res;
  try {
    res = await fetch(endpoint, {
      method: 'POST',
      headers: buildHeaders(),
      body: JSON.stringify(payload)
    });
  } catch (netErr) {
    const dt = (performance.now() - t0).toFixed(1);
    // Mixed content / TLS failure handling
    if (typeof window !== 'undefined') {
      const isHttpsPage = window.location.protocol === 'https:';
      const isHttpEndpoint = endpoint.startsWith('http://');
      if (isHttpsPage && isHttpEndpoint && !__fairchemWarnedMixed) {
        __fairchemWarnedMixed = true;
        console.warn('[fairchem] Mixed content blocked or network error accessing HTTP backend from HTTPS page. Options: 1) Run backend with TLS and use https://HOST:8000 2) Serve viewer over http:// 3) Put HTTPS reverse proxy in front. Endpoint:', endpoint);
      }
      // If we auto-upgraded to HTTPS on :8000 and that failed, attempt a single alternate HTTPS port (:8444) before falling back to HTTP.
      const isHttpsEndpoint = endpoint.startsWith('https://');
      if (isHttpsPage && isHttpsEndpoint && !__fairchemLockedEndpoint) {
        // If current endpoint is https://host:8000/... and we haven't tried 8444 yet, do so now.
        if (/:8000\//.test(endpoint) && !__fairchemTriedHttps8444) {
          __fairchemTriedHttps8444 = true;
            const alt = endpoint.replace(':8000/', ':8444/');
            console.warn('[fairchem] HTTPS :8000 failed; attempting alternate TLS port 8444', alt);
            try {
              return await postJSON(alt, payload);
            } catch (altErr) {
              // fall through to HTTP downgrade next
            }
        }
        // Candidate original HTTP (downgrade only for dev guidance; not repeated)
        const httpFallback = endpoint.replace('https://', 'http://');
        if (!__fairchemWarnedMixed) {
          console.warn('[fairchem] HTTPS auto-upgrade failed; trying HTTP fallback (may be blocked). Consider enabling TLS on backend.', endpoint);
        }
        __fairchemLockedEndpoint = httpFallback; // mark so we don't loop
        try {
          return await postJSON(httpFallback, payload); // recursion single attempt
        } catch (e2) {
          // After fallback failure, we surface original error
        }
      }
    }
    console.warn('[fairchem] ✗ network error after', dt, 'ms', netErr);
    throw netErr;
  }
  const dt = (performance.now() - t0).toFixed(1);
  let bodyText = '';
  let json;
  try {
    bodyText = await res.text();
  } catch (e) {
    console.warn('[fairchem] could not read body text', e);
  }
  if (!res.ok) {
  console.warn('[fairchem] ✗ HTTP', res.status, res.statusText, 'in', dt, 'ms body=', bodyText.slice(0, 400));
    throw new Error(`FAIR-Chem request failed: ${res.status} ${res.statusText}`);
  }
  try {
    json = bodyText ? JSON.parse(bodyText) : {};
  } catch (e) {
  console.warn('[fairchem] ✗ JSON parse error after', dt, 'ms body=', bodyText.slice(0,200));
    throw e;
  }
  if (!__fairchemLockedEndpoint && endpoint.startsWith('https://')) {
    __fairchemLockedEndpoint = endpoint; // lock successful secure endpoint
  }
  console.log('[fairchem] ← response OK in', dt, 'ms', {
    hasResults: !!json.results,
    energy: json?.results?.energy,
    forcesType: typeof json?.results?.forces,
    forcesLen: Array.isArray(json?.results?.forces) ? json.results.forces.length : 'n/a'
  });
  return json;
}

// Token handling: An optional file `public/hf-token.js` (ignored by git) can set window.HF_TOKEN = "hf_xxx".
// This avoids committing secrets. Example file content:
//   window.HF_TOKEN = 'hf_XXXXXXXXXXXXXXXXXXXXXXXXXXXX';
// Provide a token via that file or via an environment-driven injection mechanism during build/deploy.

export function createFairChemForceField({
  molecule = null,
  endpoint = (__fairchemLockedEndpoint || DEFAULT_ENDPOINT),
  units = { energy: 'eV', length: 'Angstrom' },
  defaultCharge = 0,
  defaultSpin = 1
} = {}) {
  const boundAtoms = molecule ? molecule.atoms : null;
  // Helper: extract cell if enabled & visible: returns [[ax,ay,az],[bx,by,bz],[cx,cy,cz]] or null
  function currentCellMatrix() {
    try {
      if (!molecule || !molecule.__cellState) return null;
      const cs = molecule.__cellState;
      if (!cs.visible || !cs.vectors) return null; // only pass when visible per requirement
      const { a, b, c } = cs.vectors;
      if (!a || !b || !c) return null;
      // Accept zero vectors? Skip if near-zero length
      const eps = 1e-8;
      if (a.length() < eps || b.length() < eps || c.length() < eps) return null;
      return [ [a.x,a.y,a.z], [b.x,b.y,b.z], [c.x,c.y,c.z] ];
    } catch {
      return null;
    }
  }

  async function computeRaw({ Z, xyz }) {
    validateZXYZ(Z, xyz);
    const payload = {
      atomic_numbers: Z,
      coordinates: xyz,
      charge: defaultCharge,
      spin_multiplicity: defaultSpin,
      properties: ['energy','forces']
    };
    const cell = currentCellMatrix();
    if (cell) {
      payload.cell = cell; // pass only when cell visualization enabled
    }
  console.log('[fairchem] computeRaw issuing request');
  const effEndpoint = __fairchemLockedEndpoint || endpoint;
  const data = await postJSON(effEndpoint, payload);
    const results = data.results || {};
    const energy = results.energy;
    const forces = results.forces;
    if (!Array.isArray(forces) || forces.length !== Z.length) {
  console.warn('[fairchem] invalid forces length', { got: Array.isArray(forces)? forces.length : typeof forces, expected: Z.length });
      throw new Error('Invalid forces array length from FAIR-Chem response');
    }
  console.log('[fairchem] computeRaw parsed energy/forces', { energy, forcesLen: forces.length });
    return { energy, forces };
  }

  async function compute() {
    if (!boundAtoms) throw new Error('No bound molecule provided to FAIR-Chem force field');
    const Z = atomsToZ(boundAtoms);
    const xyz = atomsToXYZ(boundAtoms);
    const t0 = performance.now();
  console.log('[fairchem] compute() start', { natoms: Z.length });
    try {
      const { energy, forces } = await computeRaw({ Z, xyz });
      const dt = (performance.now() - t0).toFixed(1);
  console.log('[fairchem] compute() done in', dt, 'ms', { energy });
      return { energy, forces: rawForcesToVectors(forces) };
    } catch (err) {
      const dt = (performance.now() - t0).toFixed(1);
      console.warn('[fairchem] compute() failed after', dt, 'ms', err);
      throw err;
    }
  }

  return {
    kind: 'fairchem',
    units,
  meta: { endpoint: (__fairchemLockedEndpoint || endpoint), provider: 'local-uma-server' },
    compute,
    computeRaw
  };
}
