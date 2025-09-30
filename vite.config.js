import { defineConfig } from 'vite';
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
    proxy: {
      // Forward API calls to local FastAPI (FairChem) backend so frontend can just fetch('/simple_calculate').
      '/simple_calculate': {
        target: 'http://127.0.0.1:8000',
        changeOrigin: true,
        // If backend expects the exact path, leave rewrite out; kept for clarity if extended later.
        // rewrite: path => path
      }
    }
  },
  preview: {
    port: 5174,
    host: '0.0.0.0'
  },
  plugins: [rootRedirectPlugin()]
});
