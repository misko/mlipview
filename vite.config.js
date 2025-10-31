import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import fs from 'fs';

// Attempt to load localhost TLS certs for HTTPS dev (optional)
function loadHttpsConfig() {
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
function rootRedirectPlugin() {
  return {
    name: 'root-redirect',
    configureServer(server) {
      server.middlewares.use((req, res, next) => {
        if (req.url === '/' || req.url === '') {
          res.statusCode = 302;
          res.setHeader('Location', '/index.html');
          return res.end();
        }
        next();
      });
    },
  };
}

// Backend configuration: prefer explicit FASTAPI_URL, else compose from FASTAPI_PROTOCOL/HOST/PORT or default to http://127.0.0.1:8000
const backendProtocol = process.env.FASTAPI_PROTOCOL || 'http';
const backendHost = process.env.FASTAPI_HOST || '127.0.0.1';
const backendPort = process.env.FASTAPI_PORT || '8000';
const backendTarget =
  process.env.FASTAPI_URL || `${backendProtocol}://${backendHost}:${backendPort}`;

export default defineConfig({
  appType: 'mpa',
  root: path.resolve(__dirname, 'public'),
  publicDir: false,
  resolve: {
    alias: {
      '@src': path.resolve(__dirname, 'src'),
    },
  },
  fs: { allow: [path.resolve(__dirname, 'src'), path.resolve(__dirname, 'public')] },
  build: {
    outDir: path.resolve(__dirname, 'dist'),
    emptyOutDir: true,
    sourcemap: true,
    rollupOptions: {
      // Multi-page build: include ws-test.html alongside index.html
      input: {
        index: path.resolve(__dirname, 'public', 'index.html'),
        wsTest: path.resolve(__dirname, 'public', 'ws-test.html'),
      },
    },
  },
  server: {
    port: 5173,
    host: '0.0.0.0', // listen on all interfaces for LAN / device access
    open: false,
    // Restrict or extend allowed host headers. Needed when accessing via hostname (e.g. http://kalman:5173).
    // Add more hostnames or patterns as needed; you can also expose an env var if desired.
    allowedHosts: ['kalman'],
    fs: {
      // allow Vite to serve modules from project root (so /src/* works)
      allow: [__dirname],
    },
    https: !!httpsCfg && !process.env.NO_VITE_HTTPS ? httpsCfg : false,
    proxy: {
      '/serve/simple': { target: backendTarget, changeOrigin: true, secure: false },
      '/serve/relax': { target: backendTarget, changeOrigin: true, secure: false },
      '/serve/md': { target: backendTarget, changeOrigin: true, secure: false },
      '/serve/health': { target: backendTarget, changeOrigin: true, secure: false },
      '^/ws($|/)': { target: backendTarget, changeOrigin: true, secure: false, ws: true },
    },
  },
  preview: {
    port: 5174,
    host: '0.0.0.0',
  },
  plugins: [
    react(),
    rootRedirectPlugin(),
    // Ensure runtime-loaded assets like molecules/*.xyz exist in dist for preview server
    {
      name: 'copy-molecules-on-build',
      closeBundle() {
        try {
          const src = path.resolve(__dirname, 'public', 'molecules');
          const dst = path.resolve(__dirname, 'dist', 'molecules');
          if (fs.existsSync(src)) {
            fs.mkdirSync(dst, { recursive: true });
            // Node >=16 supports cpSync
            if (fs.cpSync) {
              fs.cpSync(src, dst, { recursive: true });
            } else {
              // Fallback: shallow copy files only (unlikely on our Node engine)
              for (const f of fs.readdirSync(src)) {
                const s = path.join(src, f);
                const d = path.join(dst, f);
                if (fs.statSync(s).isFile()) fs.copyFileSync(s, d);
              }
            }
          }
        } catch (e) {
          console.warn('[vite][copy-molecules] failed:', e?.message || e);
        }
      },
    },
    {
      name: 'copy-examples-on-build',
      closeBundle() {
        try {
          const src = path.resolve(__dirname, 'public', 'examples');
          const dst = path.resolve(__dirname, 'dist', 'examples');
          if (fs.existsSync(src)) {
            fs.mkdirSync(dst, { recursive: true });
            if (fs.cpSync) {
              fs.cpSync(src, dst, { recursive: true });
            } else {
              for (const f of fs.readdirSync(src)) {
                const s = path.join(src, f);
                const d = path.join(dst, f);
                if (fs.statSync(s).isDirectory()) {
                  fs.mkdirSync(d, { recursive: true });
                  for (const inner of fs.readdirSync(s)) {
                    const si = path.join(s, inner);
                    const di = path.join(d, inner);
                    if (fs.statSync(si).isFile()) fs.copyFileSync(si, di);
                  }
                } else if (fs.statSync(s).isFile()) {
                  fs.copyFileSync(s, d);
                }
              }
            }
          }
        } catch (e) {
          console.warn('[vite][copy-examples] failed:', e?.message || e);
        }
      },
    },
    {
      name: 'copy-vendor-on-build',
      closeBundle() {
        try {
          const src = path.resolve(__dirname, 'public', 'vendor');
          const dst = path.resolve(__dirname, 'dist', 'vendor');
          if (!fs.existsSync(src)) return;
          fs.mkdirSync(dst, { recursive: true });
          if (fs.cpSync) {
            fs.cpSync(src, dst, { recursive: true });
          } else {
            const entries = fs.readdirSync(src);
            for (const entry of entries) {
              const s = path.join(src, entry);
              const d = path.join(dst, entry);
              const stat = fs.statSync(s);
              if (stat.isDirectory()) {
                fs.mkdirSync(d, { recursive: true });
                const inner = fs.readdirSync(s);
                for (const file of inner) {
                  const sf = path.join(s, file);
                  const df = path.join(d, file);
                  if (fs.statSync(sf).isFile()) fs.copyFileSync(sf, df);
                }
              } else if (stat.isFile()) {
                fs.copyFileSync(s, d);
              }
            }
          }
        } catch (e) {
          console.warn('[vite][copy-vendor] failed:', e?.message || e);
        }
      },
    },
  ],
  optimizeDeps: {},
});
