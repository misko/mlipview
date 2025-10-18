// REST endpoints removed: enforce WS usage only
const removedMsg = '[REST removed] Use WebSocket API via fairchem_ws_client.js';
export async function getEndpoint(kind){ throw new Error(removedMsg+' (getEndpoint:'+kind+')'); }
export function getEndpointSync(kind){ throw new Error(removedMsg+' (getEndpointSync:'+kind+')'); }
export const __endpointMaps = {};