const express = require("express");
const path = require("path");
const https = require("https");
const http = require("http");
const fs = require("fs");

const app = express();
// Parse JSON bodies for proxy forwarding
app.use(express.json({ limit: '1mb' }));
const PORT = process.env.PORT || 3000;
const HTTPS_PORT = process.env.HTTPS_PORT || 8443;

// Add security headers for WebXR (must come before static so they apply to assets too)
app.use((req, res, next) => {
  // Required for WebXR
  res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
  res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
  
  // Enable SharedArrayBuffer if needed
  res.setHeader('Cross-Origin-Resource-Policy', 'cross-origin');
  
  // Security headers
  res.setHeader('X-Frame-Options', 'SAMEORIGIN');
  res.setHeader('X-Content-Type-Options', 'nosniff');
  
  next();
});

// Middleware to serve static files
app.use(express.static(path.join(__dirname, "public")));
// New architecture preview served under /new
app.use('/new', express.static(path.join(__dirname, 'public_new')));

// API endpoint to list available molecule XYZ files
app.get('/api/molecules', (req, res) => {
  const molDir = path.join(__dirname, 'public', 'molecules');
  fs.readdir(molDir, (err, files) => {
    if (err) {
      console.error('Error reading molecules directory', err);
      return res.status(500).json({ error: 'Failed to read molecules directory' });
    }
    const xyzFiles = files.filter(f => f.toLowerCase().endsWith('.xyz'));
    // Return both filename and a short name (without extension)
    const molecules = xyzFiles.map(f => ({ file: f, name: f.replace(/\.xyz$/i, '') }));
    res.json({ molecules });
  });
});

// --- UMA FAIR-Chem local backend proxy (same-origin to avoid mixed content) ---
// Set UMA_BACKEND (default http://127.0.0.1:8000) to forward compute requests.
const UMA_BACKEND = process.env.UMA_BACKEND || 'http://127.0.0.1:8000';
const { URL } = require('url');
const backendURL = new URL(UMA_BACKEND);
const backendProto = backendURL.protocol === 'https:' ? https : http;

function forwardUMA(pathname) {
  return (req, res) => {
    const t0 = Date.now();
    const targetPath = pathname;
    const payload = req.body || {};
    const opts = {
      hostname: backendURL.hostname,
      port: backendURL.port,
      path: targetPath,
      method: 'POST',
      headers: { 'Content-Type': 'application/json' }
    };
    const prox = backendProto.request(opts, (pres) => {
      let buf = '';
      pres.on('data', (d) => buf += d);
      pres.on('end', () => {
        const ms = Date.now() - t0;
        res.status(pres.statusCode || 500);
        res.setHeader('Content-Type', pres.headers['content-type'] || 'application/json');
        res.send(buf);
        console.log(`[uma-proxy] ${pathname} ${pres.statusCode} in ${ms} ms`);
      });
    });
    prox.on('error', (e) => {
      const ms = Date.now() - t0;
      console.error(`[uma-proxy] error ${pathname} after ${ms} ms`, e.message);
      res.status(502).json({ error: 'backend_unreachable', detail: e.message });
    });
    try {
      prox.end(JSON.stringify(payload));
    } catch (e) {
      console.error('[uma-proxy] serialize error', e);
      res.status(500).json({ error: 'proxy_serialize_failed' });
    }
  };
}

app.post('/simple_calculate', forwardUMA('/simple_calculate'));
app.post('/calculate', forwardUMA('/calculate'));

// Start HTTP server
app.listen(PORT, () => {
  console.log(`HTTP Server running at http://localhost:${PORT}`);
});

// Try to start HTTPS server for WebXR (self-signed cert)
try {
  // Create self-signed certificate if it doesn't exist
  let key, cert;
  try {
    key = fs.readFileSync(path.join(__dirname, 'localhost-key.pem'));
    cert = fs.readFileSync(path.join(__dirname, 'localhost-cert.pem'));
  } catch (e) {
    console.log('HTTPS certificates not found. WebXR requires HTTPS.');
    console.log('To enable VR on Meta Quest, you need to:');
    console.log('1. Generate SSL certificates or use a service like ngrok');
    console.log('2. Access the site via HTTPS');
    console.log('3. Allow "insecure content" in browser if using self-signed certs');
    key = cert = null;
  }
  
  if (key && cert) {
    const httpsServer = https.createServer({ key, cert }, app);
    httpsServer.listen(HTTPS_PORT, () => {
      console.log(`HTTPS Server running at https://localhost:${HTTPS_PORT}`);
      console.log('Use the HTTPS URL for WebXR/VR functionality');
    });
  }
} catch (error) {
  console.log('Could not start HTTPS server:', error.message);
}
