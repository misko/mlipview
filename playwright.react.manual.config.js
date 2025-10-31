import baseConfig from './playwright.config.js';

const reactProject = baseConfig.projects?.find((proj) => proj.name === 'react');

if (!reactProject) {
  throw new Error('React project not defined in base Playwright config');
}

export default {
  ...baseConfig,
  globalSetup: undefined,
  globalTeardown: undefined,
  projects: [reactProject],
};
