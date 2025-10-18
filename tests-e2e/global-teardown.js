export default async function globalTeardown() {
  for (const key of ['__WS_PID__', '__VITE_PREVIEW_PID__']) {
    const pid = Number(process.env[key]);
    if (pid) {
      try { process.kill(pid, 'SIGTERM'); } catch {}
    }
  }
}
