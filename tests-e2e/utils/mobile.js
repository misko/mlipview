export async function singleTouchDrag(page, { start, end }) {
  await page.evaluate(({ start, end }) => {
    const canvas = document.getElementById('viewer');
    if (!canvas) return;

    const createTouch = (id, x, y) => {
      const init = {
        identifier: id,
        target: canvas,
        clientX: x,
        clientY: y,
        pageX: x,
        pageY: y,
        screenX: x,
        screenY: y,
        radiusX: 8,
        radiusY: 8,
        rotationAngle: 0,
        force: 0.5,
      };
      try {
        return new Touch(init);
      } catch {
        return init;
      }
    };

    const dispatch = (type, touchesArray, changedArray) => {
      const eventInit = {
        cancelable: true,
        bubbles: true,
        touches: touchesArray,
        targetTouches: touchesArray,
        changedTouches: changedArray,
      };
      let event;
      try {
        event = new TouchEvent(type, eventInit);
      } catch {
        event = document.createEvent('Event');
        event.initEvent(type, true, true);
        event.touches = eventInit.touches;
        event.targetTouches = eventInit.targetTouches;
        event.changedTouches = eventInit.changedTouches;
      }
      canvas.dispatchEvent(event);
    };

    const id = Date.now() % 10000;
    const startTouch = createTouch(id, start.x, start.y);
    dispatch('touchstart', [startTouch], [startTouch]);
    const moveTouch = createTouch(id, end.x, end.y);
    dispatch('touchmove', [moveTouch], [moveTouch]);
    dispatch('touchend', [], [moveTouch]);
  }, { start, end });
}
