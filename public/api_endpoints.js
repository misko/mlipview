// Modern Ray Serve endpoints only (legacy removed from tests).
const modern = { simple:'/serve/simple', relax:'/serve/relax', md:'/serve/md', health:'/serve/health' };
export async function getEndpoint(kind){ if(!(kind in modern)) throw new Error('unknown endpoint '+kind); return modern[kind]; }
export function getEndpointSync(kind){ if(!(kind in modern)) throw new Error('unknown endpoint '+kind); return modern[kind]; }
export const __endpointMaps = { modern };