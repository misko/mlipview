#!/usr/bin/env node
/**
 * Fails (exit 1) if any source file (outside the deprecated folder itself) references
 * the deprecated directory name `public_to_be_deleted`.
 */
const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, '..');
const DEPRECATED = 'public_to_be_deleted';

function walk(dir, acc = []) {
  for (const entry of fs.readdirSync(dir)) {
    if (entry.startsWith('.')) continue;
    const full = path.join(dir, entry);
    const rel = path.relative(ROOT, full);
    const stat = fs.statSync(full);
    if (stat.isDirectory()) {
      // Skip node_modules and the deprecated directory itself
      if (entry === 'node_modules' || entry === DEPRECATED || entry === '.git') continue;
      walk(full, acc);
    } else if (stat.isFile()) {
      // Only check text-like files
      if (!/\.(js|jsx|ts|tsx|json|html|css|md)$/i.test(entry)) continue;
      acc.push(rel);
    }
  }
  return acc;
}

let failed = false;
for (const rel of walk(ROOT)) {
  const content = fs.readFileSync(path.join(ROOT, rel), 'utf8');
  if (content.includes(DEPRECATED)) {
    console.error(`[DEPRECATED GUARD] Reference to '${DEPRECATED}' found in ${rel}`);
    failed = true;
  }
}

if (failed) {
  console.error(
    '\n[DEPRECATED GUARD] One or more references to the deprecated folder were found. Please remove them.'
  );
  process.exit(1);
} else {
  console.log('[DEPRECATED GUARD] OK: no references to', DEPRECATED);
}
