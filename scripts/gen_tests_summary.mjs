#!/usr/bin/env node
/**
 * Generate tests_summary.md for UI tests.
 * - Scans tests/ and tests-browser/ for *.spec.js files
 * - Extracts first describe() title as intent summary (fallback: filename)
 * - Heuristically tags area from filename
 * - Detects REST vs WS usage by grepping for legacy endpoints or WS hooks
 * - Merges Jest pass/fail status from test-results/jest-ui.json when available
 */
import fs from 'fs';
import path from 'path';

const root = process.cwd();
const uiDirs = [
  path.join(root, 'tests'),
  path.join(root, 'tests-browser'),
  path.join(root, 'tests-e2e'),
];

function listSpecFiles(dir) {
  const out = [];
  function walk(p) {
    if (!fs.existsSync(p)) return;
    const st = fs.statSync(p);
    if (st.isDirectory()) {
      for (const name of fs.readdirSync(p)) walk(path.join(p, name));
    } else if (p.endsWith('.spec.js')) {
      out.push(p);
    }
  }
  walk(dir);
  return out;
}

function readText(p) {
  try {
    return fs.readFileSync(p, 'utf8');
  } catch {
    return '';
  }
}

function firstDescribeTitle(src) {
  const m = src.match(/describe\(\s*['"]([^'"]+)['"]/);
  return m ? m[1] : null;
}

function detectArea(file) {
  const f = path.basename(file).toLowerCase();
  if (f.includes('ws')) return 'websocket/protocol';
  if (f.startsWith('md') || f.includes('md')) return 'md/simulation';
  if (f.includes('relax')) return 'relaxation';
  if (f.includes('energy')) return 'energy/plot';
  if (f.includes('force')) return 'forces/render/cache';
  if (f.includes('pbc') || f.includes('cell')) return 'pbc/cell';
  if (f.includes('vr') || f.includes('xr')) return 'vr/xr';
  if (f.includes('mobile')) return 'mobile';
  if (f.includes('select')) return 'selection';
  if (f.includes('bond')) return 'bond/rotate';
  if (f.includes('molecule')) return 'molecule/load/state';
  if (f.includes('xyz')) return 'xyz/parser';
  return 'ui/core';
}

function designAlignment(src) {
  const hasREST = /\/serve\/(simple|relax)|INIT_SYSTEM|SIMPLE_CALCULATE|requestSimpleCalculate/.test(src);
  const hasWS = /getWS\(|userInteraction\(|startSimulation\(|injectTestResult\(|setTestHook\(/.test(src);
  if (hasREST && !hasWS) return 'Legacy REST usage (needs migration)';
  if (hasREST && hasWS) return 'Mixed (migrate to WS-only)';
  if (hasWS) return 'WS aligned (protobuf)';
  return 'Neutral (non-protocol UI)';
}

function loadJestStatus() {
  const p = path.join(root, 'test-results', 'jest-ui.json');
  if (!fs.existsSync(p)) return {};
  let json;
  try {
    const raw = fs.readFileSync(p, 'utf8');
    json = JSON.parse(raw);
  } catch (e) {
    return {};
  }
  const map = {};
  for (const tr of json.testResults || []) {
    const rel = tr.name.replace(root + '/', '');
    map[rel] = tr.status === 'passed' ? 'Pass' : tr.status === 'failed' ? 'Fail' : tr.status;
  }
  return map;
}

const jestStatus = loadJestStatus();

const files = uiDirs.flatMap(listSpecFiles).map((p) => path.relative(root, p)).sort();

const rows = [];
for (const rel of files) {
  const abs = path.join(root, rel);
  const src = readText(abs);
  const area = detectArea(rel);
  const title = firstDescribeTitle(src) || path.basename(rel);
  const align = designAlignment(src);
  const status = jestStatus[rel] || 'Unknown';
  rows.push({ file: rel, area, title, align, status });
}

function toMdTable(rows) {
  const header = ['File', 'Area', 'Test summary', 'Design alignment', 'Currently passing'];
  const lines = [
    `| ${header.join(' | ')} |`,
    `| ${header.map(() => '---').join(' | ')} |`,
  ];
  for (const r of rows) {
    lines.push(`| ${r.file} | ${r.area} | ${r.title.replace(/\|/g, '\\|')} | ${r.align} | ${r.status} |`);
  }
  return lines.join('\n');
}

const outPath = path.join(root, 'tests_summary.md');
const content = `# UI Tests Summary\n\nGenerated: ${new Date().toISOString()}\n\n${toMdTable(rows)}\n`;
fs.writeFileSync(outPath, content);
console.log(`Wrote ${outPath} with ${rows.length} rows.`);
