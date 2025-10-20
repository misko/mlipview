#!/usr/bin/env node
/**
 * Fails if any mlipviewer2 source file imports legacy interaction modules from root public.
 * Converted to ESM to align with package type: module.
 */
import fs from 'fs';
import path from 'path';

const ROOT = process.cwd();
const SRC_DIR = path.join(ROOT, 'mlipviewer2'); // legacy folder (may be absent post-promotion)
const FORBIDDEN = ['public/interaction.js', 'public/bond_pick.js', 'public/selection_state.js'];

let failures = [];

function walk(dir) {
  for (const entry of fs.readdirSync(dir)) {
    const full = path.join(dir, entry);
    const stat = fs.statSync(full);
    if (stat.isDirectory()) {
      walk(full);
      continue;
    }
    if (!/\.(js|mjs|cjs|ts)$/i.test(entry)) continue;
    const rel = path.relative(ROOT, full);
    const text = fs.readFileSync(full, 'utf8');
    FORBIDDEN.forEach((fb) => {
      // naive pattern match is fine here; ensure not inside a comment? Keep simple, dev can refine.
      if (text.includes(fb)) {
        failures.push(`${rel} -> ${fb}`);
      }
    });
  }
}

if (fs.existsSync(SRC_DIR)) {
  walk(SRC_DIR);
} else {
  // If directory no longer exists, treat as success (no legacy code remains to scan).
}

if (failures.length) {
  console.error('[legacy-import-scan] Forbidden legacy import references found:');
  failures.forEach((f) => console.error('  -', f));
  process.exit(1);
} else {
  console.log('[legacy-import-scan] OK: no forbidden legacy imports in mlipviewer2');
}
