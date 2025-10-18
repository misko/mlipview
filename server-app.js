// server-app.js - exported Express app (no side-effect listeners) for tests & server.js
import express from 'express';
import path from 'path';
// Simpler: rely on process.cwd() as project root for static hosting (stable in tests & runtime).
const ROOT_DIR = process.cwd();

export function createApp() {
  const app = express();
  app.use(express.static(path.join(ROOT_DIR, 'public')));
  // Expose google-protobuf browser runtime directly from node_modules for classic <script> usage
  // This avoids ESM wrappers and ensures window.goog/jspb globals are defined as expected by generated stubs.
  app.use('/vendor/google-protobuf', express.static(path.join(ROOT_DIR, 'node_modules', 'google-protobuf')));
  app.get('/health', (_req, res) => res.json({ ok: true }));
  return app;
}

export default createApp;
