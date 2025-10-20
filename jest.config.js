export default {
  testEnvironment: 'node',
  // Include tests root only; public files imported via relative paths. (If we need module resolution, add moduleDirectories.)
  roots: ['<rootDir>/tests', '<rootDir>/tests-browser'],
  transform: {
    '^.+\\.js$': 'babel-jest',
  },
  moduleNameMapper: {
    '^/proto/(.*)$': '<rootDir>/public/proto/$1',
  },
  setupFilesAfterEnv: ['<rootDir>/tests/jest.setup.js', '<rootDir>/tests/jest.pertest.setup.js'],
  globalSetup: '<rootDir>/tests/globalSetup.cjs',
  globalTeardown: '<rootDir>/tests/globalTeardown.cjs',
  // Per-test timeout (ms); individual long-running tests should override explicitly if needed.
  testTimeout: 20000,
  // Force exit so a leaking handle never leaves CI hanging; keep detectOpenHandles on temporarily to surface offenders.
  // TODO: After leaks fixed and verified, consider removing forceExit + detectOpenHandles for stricter behavior.
  forceExit: true,
  detectOpenHandles: true,
};
