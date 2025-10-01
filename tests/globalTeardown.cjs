module.exports = async () => {
  const procs = [global.__MLIP_NODE_SERVER, global.__MLIP_PY_SERVER];
  for (const p of procs) {
    if (p && !p.killed) {
      try { p.kill('SIGTERM'); } catch {}
    }
  }
};
