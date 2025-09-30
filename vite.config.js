import { defineConfig } from 'vite';
import path from 'path';

// We serve static assets out of public/. Index.html already lives there and uses module scripts.
// Output will go to dist/ which can be hosted statically (no Node server required for basic usage).
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
    open: false
  },
  preview: {
    port: 5174
  }
});
