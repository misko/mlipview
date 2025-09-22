const express = require("express");
const path = require("path");
const https = require("https");
const http = require("http");
const fs = require("fs");

const app = express();
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
