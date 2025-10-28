// tests run against Vite preview (dist/) unless BASE_URL is provided
import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './tests-e2e',
  testIgnore: '**/*.integration.spec.js',
  timeout: 10_000,
  fullyParallel: false,
  globalSetup: './tests-e2e/global-setup.js',
  globalTeardown: './tests-e2e/global-teardown.js',
  workers: 1,
  use: {
    baseURL: process.env.BASE_URL || 'http://127.0.0.1:5174',
    headless: true,
    viewport: { width: 1280, height: 800 },
    ignoreHTTPSErrors: true,
  },
  retries: 0,
  reporter: [['list']],
});
