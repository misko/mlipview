// ui/energy-plot.js - Energy plotting functionality
export function createEnergyPlot() {
  const plotContainer = document.getElementById("plotContainer");
  const energyChart = document.getElementById("energyChart");
  
  console.log("[DEBUG] Plot container found:", !!plotContainer);
  console.log("[DEBUG] Energy chart canvas found:", !!energyChart);
  console.log("[DEBUG] Chart.js available:", typeof Chart !== 'undefined');

  // Return a no-op implementation if required elements are missing
  if (!plotContainer || !energyChart) {
    console.warn("[WARN] Plot elements not found, creating no-op energy plot");
    return {
      isEnabled: () => false,
      toggle: () => false,
      addInitialDataPoint: () => {},
      updateChart: () => {},
      clear: () => {},
      destroy: () => {}
    };
  }

  let showPlot = true; // Default to ON
  let chart = null;
  let lastRotationCount = 0;
  const maxDataPoints = 100;

  // Initialize plot as visible since it's on by default
  plotContainer.style.display = showPlot ? "block" : "none";

  function initializeChart() {
    console.log("[DEBUG] Initializing chart...");
    
    if (typeof Chart === 'undefined') {
      console.error("[ERROR] Chart.js is not loaded!");
      return;
    }
    
    if (!energyChart) {
      console.error("[ERROR] Energy chart canvas not found!");
      return;
    }
    
    // Destroy existing chart if it exists
    const existingChart = Chart.getChart(energyChart);
    if (existingChart) {
      console.log("[DEBUG] Destroying existing chart...");
      existingChart.destroy();
    }
    
    const ctx = energyChart.getContext('2d');
    chart = new Chart(ctx, {
      type: 'line',
      data: {
        labels: [],
        datasets: [{
          label: 'Energy',
          data: [],
          borderColor: '#f5ffd1',
          backgroundColor: 'rgba(245, 255, 209, 0.1)',
          borderWidth: 2,
          fill: false,
          tension: 0.1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          x: {
            title: {
              display: true,
              text: 'Bond Rotations',
              color: '#a4b0c0'
            },
            ticks: {
              color: '#a4b0c0'
            },
            grid: {
              color: 'rgba(164, 176, 192, 0.2)'
            }
          },
          y: {
            title: {
              display: true,
              text: 'Energy',
              color: '#a4b0c0'
            },
            ticks: {
              color: '#a4b0c0'
            },
            grid: {
              color: 'rgba(164, 176, 192, 0.2)'
            }
          }
        },
        plugins: {
          legend: {
            labels: {
              color: '#f5ffd1'
            }
          }
        },
        animation: false
      }
    });
    console.log("[DEBUG] Chart initialized:", !!chart);
  }

  function updateChart(energy, state) {
    if (!chart) return;
    
    const currentRotationCount = state.rotations.length;
    
    if (currentRotationCount > lastRotationCount) {
      lastRotationCount = currentRotationCount;
      
      chart.data.labels.push(currentRotationCount);
      chart.data.datasets[0].data.push(energy);
      
      if (chart.data.labels.length > maxDataPoints) {
        chart.data.labels.shift();
        chart.data.datasets[0].data.shift();
      }
      
      chart.update('none');
      console.log("[DEBUG] Chart updated, rotation count:", currentRotationCount, "energy:", energy);
    }
  }

  function toggle() {
    showPlot = !showPlot;
    plotContainer.style.display = showPlot ? "block" : "none";
    
    if (showPlot && !chart) {
      initializeChart();
    }
    
    return showPlot;
  }

  function addInitialDataPoint(energy, state) {
    if (chart) {
      const currentRotationCount = state.rotations.length;
      chart.data.labels.push(currentRotationCount);
      chart.data.datasets[0].data.push(energy);
      lastRotationCount = currentRotationCount;
      chart.update('none');
      console.log("[DEBUG] Added initial data point - rotation:", currentRotationCount, "energy:", energy);
    }
  }

  function clear(state) {
    if (chart) {
      chart.data.labels = [];
      chart.data.datasets[0].data = [];
      lastRotationCount = state.rotations.length;
      chart.update();
      console.log("[DEBUG] Plot data cleared, reset to rotation count:", lastRotationCount);
    }
  }

  // Initialize chart since plot is on by default
  if (showPlot) {
    console.log("[DEBUG] Initializing chart on startup...");
    initializeChart();
  }

  // Setup clear button handler
  const clearButton = document.getElementById("clearPlot");
  if (clearButton) {
    clearButton.onclick = () => {
      console.log("[DEBUG] Clear plot button clicked");
      clear(window.appState); // Will be set by main.js
    };
  } else {
    console.warn("[WARN] Clear plot button not found");
  }

  function destroy() {
    if (chart) {
      console.log("[DEBUG] Destroying chart...");
      chart.destroy();
      chart = null;
    }
  }

  return {
    isEnabled: () => showPlot,
    toggle,
    updateChart,
    addInitialDataPoint,
    clear,
    destroy,
    initializeChart
  };
}
