export default {
  testEnvironment: 'node',
  // Include tests root only; public files imported via relative paths. (If we need module resolution, add moduleDirectories.)
  roots: ['<rootDir>/tests', '<rootDir>/tests-browser', '<rootDir>/src/react'],
  transform: {
    '^.+\\.[tj]sx?$': 'babel-jest',
  },
  moduleNameMapper: {
    '^/proto/(.*)$': '<rootDir>/public/proto/$1',
  },
  moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
  setupFilesAfterEnv: ['<rootDir>/tests/jest.setup.js', '<rootDir>/tests/jest.pertest.setup.js'],
  globalSetup: '<rootDir>/tests/globalSetup.cjs',
  globalTeardown: '<rootDir>/tests/globalTeardown.cjs',
  // Per-test timeout (ms); individual long-running tests should override explicitly if needed.
  testTimeout: 20000,
  // Keep forceExit/detectOpenHandles on by default until lingering async leaks are resolved.
  // Set MLIPVIEW_JEST_FORCEEXIT=0 to opt out locally.
  forceExit: process.env.MLIPVIEW_JEST_FORCEEXIT !== '0',
  detectOpenHandles: process.env.MLIPVIEW_JEST_FORCEEXIT !== '0',
};
