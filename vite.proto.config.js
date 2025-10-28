import { defineConfig } from 'vite';

export default defineConfig({
  build: {
    lib: {
      entry: 'src/proto/index.js',
      name: 'ProtoBundle',
      formats: ['iife'],
      fileName: () => 'proto-bundle.js',
    },
    outDir: 'public/vendor',
    emptyOutDir: false,
    rollupOptions: {
      output: {
        extend: true,
      },
    },
  },
});
