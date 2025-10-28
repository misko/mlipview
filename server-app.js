// server-app.js - exported Express app (no side-effect listeners) for tests & server.js
import express from 'express';
import path from 'path';
import fs from 'fs';
// Simpler: rely on process.cwd() as project root for static hosting (stable in tests & runtime).
const ROOT_DIR = process.cwd();

export function createApp() {
  const app = express();
  // Prefer serving bundled assets from dist/ if it exists (ensures bare imports are resolved)
  const distDir = path.join(ROOT_DIR, 'dist');
  const pubDir = path.join(ROOT_DIR, 'public');
  if (fs.existsSync(distDir)) {
    app.use(express.static(distDir));
  }
  // Always expose the fallback public/ tree so legacy assets like /vr/main-vr.js stay reachable.
  app.use(express.static(pubDir));
  // Expose google-protobuf browser runtime directly from node_modules for classic <script> usage
  // This avoids ESM wrappers and ensures window.goog/jspb globals are defined as expected by generated stubs.
  app.use(
    '/vendor/google-protobuf',
    express.static(path.join(ROOT_DIR, 'node_modules', 'google-protobuf'))
  );
  app.get('/health', (_req, res) => res.json({ ok: true }));
  return app;
}

export default createApp;
