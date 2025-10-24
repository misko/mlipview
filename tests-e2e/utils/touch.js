/**
 * Dispatches a synthetic touch event against the given selector within the page.
 * @param {import('@playwright/test').Page} page
 * @param {string} selector
 * @param {'touchstart'|'touchmove'|'touchend'|'touchcancel'} type
 * @param {Array<{id?:number,x:number,y:number,force?:number,radius?:number}>} points
 */
export async function dispatchTouch(page, selector, type, points) {
  await page.evaluate(
    ({ selector, type, points }) => {
      const target = document.querySelector(selector);
      if (!target) throw new Error(`dispatchTouch: selector ${selector} not found`);

      const toTouch = (point) => {
        const detail = {
          identifier: point.id ?? 0,
          target,
          clientX: point.x,
          clientY: point.y,
          screenX: point.x,
          screenY: point.y,
          pageX: point.x,
          pageY: point.y,
          radiusX: point.radius ?? 2,
          radiusY: point.radius ?? 2,
          force: point.force ?? 1,
          rotationAngle: 0,
        };
        // Some browsers require Touch constructor to exist; fall back if missing.
        if (typeof Touch !== 'undefined') return new Touch(detail);
        const touch = document.createTouch
          ? document.createTouch(window, target, detail.identifier, detail.pageX, detail.pageY, detail.screenX, detail.screenY)
          : detail;
        return touch;
      };

      const changed = points.map(toTouch);
      const active = type === 'touchend' || type === 'touchcancel' ? [] : changed;

      let eventInit = {
        cancelable: true,
        bubbles: true,
        touches: active,
        targetTouches: active,
        changedTouches: changed,
      };

      let event;
      try {
        event = new TouchEvent(type, eventInit);
      } catch {
        event = document.createEvent('Event');
        event.initEvent(type, true, true);
        event.touches = active;
        event.targetTouches = active;
        event.changedTouches = changed;
      }

      try {
        event.__mlip_synthetic_from_touch = true;
      } catch {}

      target.dispatchEvent(event);
    },
    { selector, type, points },
  );
}
