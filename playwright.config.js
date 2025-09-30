// Basic Playwright config for e2e tests
// We run against the dev server (node server.js) assumed running on localhost:3000
// Tests can be invoked with: npx playwright test

import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './tests-e2e',
  timeout: 10_000,
  globalSetup: './tests-e2e/global-setup.js',
  use: {
    baseURL: 'http://localhost:4000',
    headless: true,
    viewport: { width: 1280, height: 800 },
    ignoreHTTPSErrors: true,
  },
  retries: 0,
  reporter: [['list']],
});
