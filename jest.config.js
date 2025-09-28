module.exports = {
  testEnvironment: 'node',
  // Updated roots: include original legacy public/tests plus new mlipviewer2 project.
  roots: [
    '<rootDir>/public',
    '<rootDir>/tests',
    '<rootDir>/mlipviewer2/public',
    '<rootDir>/mlipviewer2/tests'
  ],
  moduleFileExtensions: ['js','mjs','cjs','jsx','ts','tsx','json'],
  setupFilesAfterEnv: ['<rootDir>/tests/jest.setup.js'],
  transform: {
    '^.+\\\.(js|mjs|cjs)$': 'babel-jest'
  },
  transformIgnorePatterns: ['/node_modules/'],
};
