// vr/vr-plot.js - 3D energy plot for VR using Babylon.js meshes
export function create3DEnergyPlot(scene) {
  const maxPoints = 100;
  const plotData = [];
  let plotMesh = null;
  let axisMeshes = [];
  
  // Create 3D plot visualization
  function createPlotMesh() {
    // Clear existing plot
    if (plotMesh) {
      plotMesh.dispose();
    }
    axisMeshes.forEach(mesh => mesh.dispose());
    axisMeshes = [];
    
    if (plotData.length < 2) return;
    
    // Create line mesh for energy plot
    const points = plotData.map((data, index) => {
      const x = (index / (maxPoints - 1)) * 10 - 5; // -5 to 5
      const y = data.energy * 0.1; // Scale energy for visibility
      const z = 0;
      return new BABYLON.Vector3(x, y, z);
    });
    
    const lines = BABYLON.MeshBuilder.CreateLines("energyPlot", {
      points: points,
      updatable: true
    }, scene);
    
    lines.color = new BABYLON.Color3(0.96, 1, 0.82); // Light green
    lines.position = new BABYLON.Vector3(0, 2, -3); // Position in front of user
    
    plotMesh = lines;
    
    // Create axes
    createAxes();
    
    // Create data point spheres
    createDataPoints(points);
  }
  
  function createAxes() {
    // X-axis (time/steps)
    const xAxis = BABYLON.MeshBuilder.CreateLines("xAxis", {
      points: [
        new BABYLON.Vector3(-5, 0, 0),
        new BABYLON.Vector3(5, 0, 0)
      ]
    }, scene);
    xAxis.color = new BABYLON.Color3(0.5, 0.5, 0.5);
    xAxis.position = new BABYLON.Vector3(0, 2, -3);
    axisMeshes.push(xAxis);
    
    // Y-axis (energy)
    const yAxis = BABYLON.MeshBuilder.CreateLines("yAxis", {
      points: [
        new BABYLON.Vector3(-5, -2, 0),
        new BABYLON.Vector3(-5, 2, 0)
      ]
    }, scene);
    yAxis.color = new BABYLON.Color3(0.5, 0.5, 0.5);
    yAxis.position = new BABYLON.Vector3(0, 2, -3);
    axisMeshes.push(yAxis);
    
    // Add axis labels using GUI
    createAxisLabels();
  }
  
  function createDataPoints(points) {
    points.forEach((point, index) => {
      const sphere = BABYLON.MeshBuilder.CreateSphere(`dataPoint_${index}`, {
        diameter: 0.1
      }, scene);
      
      sphere.position = point.add(new BABYLON.Vector3(0, 2, -3));
      
      const material = new BABYLON.StandardMaterial(`pointMat_${index}`, scene);
      material.diffuseColor = new BABYLON.Color3(1, 0.8, 0.2); // Gold
      material.emissiveColor = new BABYLON.Color3(0.3, 0.2, 0.05);
      sphere.material = material;
      
      axisMeshes.push(sphere);
    });
  }
  
  function createAxisLabels() {
    // Create label planes
    const labelTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateForMesh(
      BABYLON.MeshBuilder.CreatePlane("labelPlane", {size: 1}, scene)
    );
    
    const xLabel = new BABYLON.GUI.TextBlock("xLabel", "Bond Rotations");
    xLabel.color = "white";
    xLabel.fontSize = "60px";
    labelTexture.addControl(xLabel);
    
    labelTexture.getChildren()[0].position = new BABYLON.Vector3(0, 1, -3);
  }
  
  function addDataPoint(rotationCount, energy) {
    plotData.push({ step: rotationCount, energy: energy });
    
    // Keep only last maxPoints
    if (plotData.length > maxPoints) {
      plotData.shift();
    }
    
    createPlotMesh(); // Refresh the 3D plot
  }
  
  function toggle(visible) {
    if (plotMesh) {
      plotMesh.setEnabled(visible);
    }
    axisMeshes.forEach(mesh => mesh.setEnabled(visible));
  }
  
  return {
    addDataPoint,
    toggle,
    dispose: () => {
      if (plotMesh) plotMesh.dispose();
      axisMeshes.forEach(mesh => mesh.dispose());
    }
  };
}
