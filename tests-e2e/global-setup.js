import { startServers } from '../server.js';

export default async function globalSetup() {
  // Start only HTTP server; ignore HTTPS errors already handled internally
  process.env.NO_HTTPS = '1';
  startServers();
  // Give a brief delay for server to bind
  await new Promise(r=>setTimeout(r, 500));
}
