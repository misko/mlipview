module.exports = {
  testEnvironment: 'node',
  roots: ['<rootDir>/public', '<rootDir>/tests'],
  moduleFileExtensions: ['js','mjs','cjs','jsx','ts','tsx','json'],
  transform: {
    '^.+\\.(js|mjs|cjs)$': 'babel-jest'
  },
  transformIgnorePatterns: ['/node_modules/'],
};
