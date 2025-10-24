// Integration test: HTTP GET /vr/main-vr.js returns 200 and contains initVRApp
import request from 'supertest';
import { createApp } from '../server-app.js';

describe('HTTP VR asset', () => {
  test('GET /vr/main-vr.js returns 200 & initVRApp token', async () => {
    const app = createApp();
    const res = await request(app).get('/vr/main-vr.js');
    expect(res.status).toBe(200);
    expect(res.text.includes('export async function initVRApp')).toBe(true);
    expect(res.text.includes('public_to_be_deleted')).toBe(false);
  });
});
