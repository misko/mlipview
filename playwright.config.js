// tests run against Vite preview (dist/) unless BASE_URL is provided
import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './tests-e2e',
  timeout: 10_000,
  testIgnore: ['**/*.integration.spec.js'],
  globalSetup: './tests-e2e/global-setup.js',
  globalTeardown: './tests-e2e/global-teardown.js',
  use: {
    baseURL: process.env.BASE_URL || 'http://127.0.0.1:5174',
    headless: true,
    viewport: { width: 1280, height: 800 },
    ignoreHTTPSErrors: true,
  },
  workers: 1,
  retries: 0,
  reporter: [['list']],
});
