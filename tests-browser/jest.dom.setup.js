// Additional DOM/browser-specific setup for jsdom tests.
// Provide minimal canvas + PointerEvent polyfills if not present.

if (typeof window !== 'undefined') {
  if (!window.PointerEvent) {
    class PointerEvent extends MouseEvent {}
    window.PointerEvent = PointerEvent;
  }
  if (typeof HTMLCanvasElement !== 'undefined' && !HTMLCanvasElement.prototype.getContext) {
    HTMLCanvasElement.prototype.getContext = function(kind) {
      if (kind === '2d') {
        return {
          canvas: this,
          clearRect(){}, fillRect(){}, beginPath(){}, moveTo(){}, lineTo(){}, stroke(){},
          fillStyle: '#000', strokeStyle: '#000'
        };
      }
      if (kind === 'webgl' || kind === 'webgl2') {
        // Stub object sufficient for code paths that only feature-detect.
        return {
          getExtension(){}, clear(){}, viewport(){}, createShader(){}, shaderSource(){}, compileShader(){},
          createProgram(){}, attachShader(){}, linkProgram(){}, useProgram(){}, getShaderParameter(){ return true; },
          getProgramParameter(){ return true; }
        };
      }
      return null;
    };
  }
}
