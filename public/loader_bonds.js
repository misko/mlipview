// Parse a simple INI-like ROTATABLE declaration file (e.g., ROY.BONDS)

export async function loadBondsSpec(url) {
  const res = await fetch(url, { cache: "no-cache" });
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status} ${res.statusText}`);
  const text = await res.text();
  return parseBondsSpec(text);
}

export function parseBondsSpec(text) {
  const lines = text.split(/\r?\n/);
  const items = [];
  for (const raw of lines) {
    const line = raw.trim();
    if (!line || line.startsWith("#")) continue;
    const [tag, ...rest] = line.split(/\s+/);
    if (tag.toUpperCase() !== "ROTATABLE") continue;

    const obj = { i: null, j: null, side: "j", label: "", step: 5, min: -180, max: 180 };
    for (const kv of rest) {
      const m = kv.match(/^(\w+)=(.+)$/);
      if (!m) continue;
      const k = m[1].toLowerCase();
      const v = m[2];
      if (k === "i" || k === "j") obj[k] = parseInt(v, 10);
      else if (k === "side") obj.side = v.toLowerCase() === "i" ? "i" : "j";
      else if (k === "label") obj.label = v;
      else if (k === "step") obj.step = parseFloat(v);
      else if (k === "min") obj.min = parseFloat(v);
      else if (k === "max") obj.max = parseFloat(v);
    }
    if (Number.isInteger(obj.i) && Number.isInteger(obj.j)) items.push(obj);
  }
  return items;
}
