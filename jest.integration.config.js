export default {
  testEnvironment: 'node',
  testMatch: ['<rootDir>/tests-e2e/**/*.integration.spec.js'],
  transform: {
    '^.+\\.js$': 'babel-jest'
  },
  setupFilesAfterEnv: ['<rootDir>/tests/jest.setup.js','<rootDir>/tests/jest.pertest.setup.js'],
  globalSetup: '<rootDir>/tests/globalSetup.cjs',
  globalTeardown: '<rootDir>/tests/globalTeardown.cjs',
  testTimeout: 300000,
  moduleNameMapper: {
    '^/proto/(.*)$': '<rootDir>/public/proto/$1'
  },
  forceExit: true,
  detectOpenHandles: true
};
