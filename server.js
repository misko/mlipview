import fs from 'fs';
import https from 'https';
import http from 'http';
import path from 'path';
import { fileURLToPath } from 'url';
import { createApp } from './server-app.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const PORT = process.env.PORT || 4000; // HTTP
const SSL_PORT = process.env.SSL_PORT || 4443; // HTTPS
const app = createApp();

function startServers() {
  // Always start HTTP (useful for local dev & redirect target)
  const httpServer = http.createServer(app);
  httpServer.listen(PORT, () => {
    console.log(`[mlipviewer2] HTTP  listening on  http://localhost:${PORT}`);
  });

  // Attempt HTTPS if certs exist
  const certDirCandidates = [
    path.join(__dirname, 'certs'),
    path.join(__dirname) // fallback if user placed certs at project root
  ];

  let keyPath = null;
  let certPath = null;
  for (const dir of certDirCandidates) {
    const k = path.join(dir, 'localhost-key.pem');
    const c = path.join(dir, 'localhost-cert.pem');
    if (fs.existsSync(k) && fs.existsSync(c)) { keyPath = k; certPath = c; break; }
  }

  if (keyPath && certPath) {
    try {
      const key = fs.readFileSync(keyPath);
      const cert = fs.readFileSync(certPath);
      const httpsServer = https.createServer({ key, cert }, app);
      httpsServer.listen(SSL_PORT, () => {
        console.log(`[mlipviewer2] HTTPS listening on https://localhost:${SSL_PORT}`);
      });
    } catch (e) {
      console.warn('[mlipviewer2] Failed to start HTTPS server:', e.message);
    }
  } else {
    console.warn('[mlipviewer2] HTTPS certs not found (expected localhost-key.pem & localhost-cert.pem in ./certs). Skipping HTTPS.');
  }
}

// Only auto-start when executed directly (node server.js), not when imported (tests).
if (process.argv[1] === __filename) {
  startServers();
}

export { app, startServers };
