// Test wrapper invoking the legacy import scan script so it runs under jest.
const { spawnSync } = require('node:child_process');
const path = require('node:path');

describe('no legacy import references inside mlipviewer2', () => {
  it('passes the legacy import guard script', () => {
    const script = path.join(process.cwd(), 'scripts', 'check-no-legacy-imports.js');
    const r = spawnSync('node', [script], { encoding: 'utf8' });
    if (r.status !== 0) {
      console.error(r.stdout, r.stderr);
    }
    expect(r.status).toBe(0);
    expect(r.stdout).toContain('OK: no forbidden legacy imports');
  });
});
