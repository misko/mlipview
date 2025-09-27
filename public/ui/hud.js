// ui/hud.js - Energy HUD and control buttons
export function createEnergyHUD() {
  const hud = document.createElement("div");
  hud.style.position = "absolute";
  hud.style.top = "12px";
  hud.style.right = "12px";
  hud.style.padding = "8px 10px";
  hud.style.background = "rgba(15,18,24,0.7)";
  hud.style.border = "1px solid rgba(255,255,255,0.08)";
  hud.style.borderRadius = "8px";
  hud.style.color = "#f5ffd1";
  hud.style.font = "14px system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
  
  hud.innerHTML = `
    <div style="font-size:12px; color:#a4b0c0; margin-bottom:4px">Mock Energy</div>
    <div id="energyVal" style="font-weight:700; font-size:22px; letter-spacing:0.4px; margin-bottom:8px;">0.000</div>
    <div style="display:flex; gap:6px; flex-wrap:wrap;">
      <button id="btnForces">Forces: ON</button>
      <button id="btnLenMode" title="Toggle between normalized length and true magnitude">Length: Normalized</button>
      <button id="btnPlot">Plot: ON</button>
      <button id="btnMolecules" title="Select molecule">Molecules</button>
      <button id="btnFFFair" title="Switch to FAIR-Chem backend">FF: FAIR</button>
      <button id="btnFFLJ" title="Switch to Lennard-Jones backend">FF: LJ</button>
    </div>
  `;
  
  document.body.appendChild(hud);
  
  return {
    container: hud,
    energyVal: hud.querySelector("#energyVal"),
    btnForces: hud.querySelector("#btnForces"),
    btnLenMode: hud.querySelector("#btnLenMode"),
    btnPlot: hud.querySelector("#btnPlot")
    ,btnMolecules: hud.querySelector('#btnMolecules')
    ,btnFFFair: hud.querySelector('#btnFFFair')
    ,btnFFLJ: hud.querySelector('#btnFFLJ')
  };
}

export function createStateBar() {
  const stateBar = document.createElement("div");
  stateBar.style.position = "absolute";
  stateBar.style.top = "12px";
  stateBar.style.left = "50%";
  stateBar.style.transform = "translateX(-50%)";
  stateBar.style.background = "rgba(15,18,24,0.7)";
  stateBar.style.border = "1px solid rgba(255,255,255,0.08)";
  stateBar.style.borderRadius = "8px";
  stateBar.style.padding = "6px 8px";
  stateBar.style.display = "flex";
  stateBar.style.gap = "8px";
  stateBar.style.color = "#cfe3ff";
  stateBar.style.font = "12px system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
  
  stateBar.innerHTML = `
    <button id="btnRebuild">Rebuild from rotations</button>
    <button id="btnExport">Export XYZ</button>
  `;
  
  document.body.appendChild(stateBar);
  
  return {
    container: stateBar,
    btnRebuild: stateBar.querySelector("#btnRebuild"),
    btnExport: stateBar.querySelector("#btnExport")
  };
}
