import { defineConfig } from 'vite';
import fs from 'fs';

// Attempt to load localhost TLS certs for HTTPS dev (optional)
function loadHttpsConfig(){
  const keyPath = './localhost-key.pem';
  const certPath = './localhost-cert.pem';
  try {
    if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {
      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    }
  } catch {}
  return null;
}
const httpsCfg = loadHttpsConfig();
import path from 'path';

// We serve static assets out of public/. Index.html already lives there and uses module scripts.
// Output will go to dist/ which can be hosted statically (no Node server required for basic usage).
// Simple plugin to redirect bare '/' requests explicitly to '/index.html' (helps when hostname restrictions / caching caused prior mismatch)
function rootRedirectPlugin(){
  return {
    name: 'root-redirect',
    configureServer(server){
      server.middlewares.use((req,res,next)=>{
        if(req.url === '/' || req.url === ''){
          res.statusCode = 302;
          res.setHeader('Location','/index.html');
          return res.end();
        }
        next();
      });
    }
  };
}

// Backend configuration: prefer explicit FASTAPI_URL, else compose from FASTAPI_PROTOCOL/HOST/PORT or default to http://127.0.0.1:8000
const backendProtocol = process.env.FASTAPI_PROTOCOL || 'http';
const backendHost = process.env.FASTAPI_HOST || '127.0.0.1';
const backendPort = process.env.FASTAPI_PORT || '8000';
const backendTarget = process.env.FASTAPI_URL || `${backendProtocol}://${backendHost}:${backendPort}`;

export default defineConfig({
  root: path.resolve(__dirname, 'public'),
  publicDir: path.resolve(__dirname, 'public'),
  build: {
    outDir: path.resolve(__dirname, 'dist'),
    emptyOutDir: true,
    sourcemap: true,
    rollupOptions: {
      // Entry is index.html in root; Vite auto-detects. We keep explicit config minimal.
    }
  },
  server: {
    port: 5173,
    host: '0.0.0.0', // listen on all interfaces for LAN / device access
    open: false,
    // Restrict or extend allowed host headers. Needed when accessing via hostname (e.g. http://kalman:5173).
    // Add more hostnames or patterns as needed; you can also expose an env var if desired.
    allowedHosts: ['kalman'],
  https: !!httpsCfg && !process.env.NO_VITE_HTTPS ? httpsCfg : false,
  proxy: {
    '/serve/simple': { target: backendTarget, changeOrigin: true, secure: false },
    '/serve/relax': { target: backendTarget, changeOrigin: true, secure: false },
    '/serve/md': { target: backendTarget, changeOrigin: true, secure: false },
    '/serve/health': { target: backendTarget, changeOrigin: true, secure: false }
    }
  },
  preview: {
    port: 5174,
    host: '0.0.0.0'
  },
  plugins: [rootRedirectPlugin()]
});
