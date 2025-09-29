export default {
  testEnvironment: 'node',
  // Include tests root only; public files imported via relative paths. (If we need module resolution, add moduleDirectories.)
  roots: ['<rootDir>/tests'],
  transform: {
    '^.+\\.js$': 'babel-jest'
  },
  setupFilesAfterEnv: ['<rootDir>/tests/jest.setup.js']
};
