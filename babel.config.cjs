module.exports = {
  presets: [
    ["@babel/preset-env", { targets: { node: "18" } }],
    ["@babel/preset-typescript", { allowDeclareFields: true }],
    ["@babel/preset-react", { runtime: "automatic", importSource: "react" }]
  ]
};
