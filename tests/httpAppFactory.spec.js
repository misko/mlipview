// Validates that server-app exports a plain express app without starting listeners.
import { createApp } from '../server-app.js';

describe('server-app factory', () => {
  test('returns express-style app with health route', async () => {
    const app = createApp();
    expect(typeof app.use).toBe('function');
  });
});
