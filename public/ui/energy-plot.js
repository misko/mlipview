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
  let internalStep = 0; // internal monotonic counter independent of removed rotation history
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

  function recordStep(energy) {
    if (!chart) return;
    internalStep += 1;
    chart.data.labels.push(internalStep);
    chart.data.datasets[0].data.push(energy);
    if (chart.data.labels.length > maxDataPoints) {
      chart.data.labels.shift();
      chart.data.datasets[0].data.shift();
    }
    chart.update('none');
  }

  function toggle() {
    showPlot = !showPlot;
    plotContainer.style.display = showPlot ? "block" : "none";
    
    if (showPlot && !chart) {
      initializeChart();
    }
    
    return showPlot;
  }

  function addInitialDataPoint(energy) {
    if (chart) {
      internalStep = 0;
      chart.data.labels.push(internalStep);
      chart.data.datasets[0].data.push(energy);
      chart.update('none');
    }
  }

  function clear() {
    if (chart) {
      chart.data.labels = [];
      chart.data.datasets[0].data = [];
      internalStep = 0;
      chart.update();
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
  recordStep,
    addInitialDataPoint,
    clear,
    destroy,
    initializeChart
  };
}
