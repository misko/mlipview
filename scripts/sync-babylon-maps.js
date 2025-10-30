#!/usr/bin/env node

/**
 * Copies Babylon.js source maps from node_modules into the vendor directory so
 * Vite can serve them during development (avoids noisy 404 warnings).
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const vendorRoot = path.resolve(__dirname, '../public/vendor/babylon');
const modulesRoot = path.resolve(__dirname, '../node_modules');

const artifacts = [
  { from: 'babylonjs/babylon.js.map', to: 'babylon.js.map' },
  { from: 'babylonjs-loaders/babylonjs.loaders.js.map', to: 'babylonjs.loaders.js.map' },
  { from: 'babylonjs-gui/babylon.gui.js.map', to: 'babylon.gui.js.map' },
  { from: '@babylonjs/gui/babylon.gui.min.js.map', to: 'babylon.gui.min.js.map' },
  { from: '@babylonjs/loaders/babylonjs.loaders.min.js.map', to: 'babylonjs.loaders.min.js.map' },
];

function ensureDir(dir) {
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir, { recursive: true });
  }
}

function copyIfExists(source, target) {
  try {
    ensureDir(path.dirname(target));
    fs.copyFileSync(source, target);
    console.log(`[sync-babylon-maps] synced ${path.basename(target)}`);
  } catch (err) {
    if (err.code === 'ENOENT') {
      console.warn(`[sync-babylon-maps] skipped ${path.basename(target)} (missing in node_modules)`);
    } else {
      console.warn(`[sync-babylon-maps] failed ${path.basename(target)}: ${err.message}`);
    }
  }
}

artifacts.forEach(({ from, to }) => {
  const source = path.join(modulesRoot, from);
  const target = path.join(vendorRoot, to);
  copyIfExists(source, target);
});
