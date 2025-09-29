// server-app.js - exported Express app (no side-effect listeners) for tests & server.js
import express from 'express';
import path from 'path';
// Simpler: rely on process.cwd() as project root for static hosting (stable in tests & runtime).
const ROOT_DIR = process.cwd();

export function createApp() {
  const app = express();
  app.use(express.static(path.join(ROOT_DIR, 'public')));
  app.get('/health', (_req, res) => res.json({ ok: true }));
  return app;
}

export default createApp;
