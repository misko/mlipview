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
      <button id="btnRelax" title="Start/stop geometry relaxation">Relax</button>
    </div>
    <div id="relaxControls" style="margin-top:8px; display:none; padding:8px; background:rgba(0,0,0,0.3); border-radius:4px;">
      <div style="font-size:11px; color:#a4b0c0; margin-bottom:4px;">Relaxation Status</div>
      <div id="relaxStatus" style="font-size:12px; margin-bottom:6px; color:#cfe3ff;">Idle</div>
      <div style="display:flex; gap:4px; align-items:center; margin-bottom:4px;">
        <label style="font-size:11px; color:#cfe3ff;">Method:</label>
        <select id="relaxOptimizer" style="font-size:11px; background:#2a3441; color:#cfe3ff; border:1px solid #444; border-radius:2px;">
          <option value="steepest_descent">Steepest Descent</option>
          <option value="conjugate_gradient">Conjugate Gradient</option>
        </select>
      </div>
      <div style="display:flex; gap:4px; align-items:center; margin-bottom:4px;">
        <label style="font-size:11px; color:#cfe3ff;">Step Size:</label>
        <input id="relaxStepSize" type="range" min="0.001" max="0.1" step="0.001" value="0.01" style="flex:1;">
        <span id="relaxStepSizeVal" style="font-size:11px; color:#cfe3ff; min-width:40px;">0.01</span>
      </div>
      <div style="display:flex; gap:4px; align-items:center;">
        <label style="font-size:11px; color:#cfe3ff;">Max Steps:</label>
        <input id="relaxMaxSteps" type="number" min="1" max="10000" value="1000" style="flex:1; font-size:11px; background:#2a3441; color:#cfe3ff; border:1px solid #444; border-radius:2px; padding:2px;">
      </div>
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
    ,btnRelax: hud.querySelector('#btnRelax')
    ,relaxControls: hud.querySelector('#relaxControls')
    ,relaxStatus: hud.querySelector('#relaxStatus')
    ,relaxOptimizer: hud.querySelector('#relaxOptimizer')
    ,relaxStepSize: hud.querySelector('#relaxStepSize')
    ,relaxStepSizeVal: hud.querySelector('#relaxStepSizeVal')
    ,relaxMaxSteps: hud.querySelector('#relaxMaxSteps')
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
  
  // Rebuild / Export controls removed
  stateBar.innerHTML = ``;
  
  document.body.appendChild(stateBar);
  
  return { container: stateBar };
}
