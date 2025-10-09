// Modern Ray Serve endpoints only (legacy removed from tests).
// Include cache-reuse variants to avoid resending atoms when geometry hasn't changed.
const modern = {
	simple: '/serve/simple',
	simple_from_cache: '/serve/simple_from_cache',
	relax: '/serve/relax',
	relax_from_cache: '/serve/relax_from_cache',
	md: '/serve/md',
	md_from_cache: '/serve/md_from_cache',
	health: '/serve/health'
};
export async function getEndpoint(kind){ if(!(kind in modern)) throw new Error('unknown endpoint '+kind); return modern[kind]; }
export function getEndpointSync(kind){ if(!(kind in modern)) throw new Error('unknown endpoint '+kind); return modern[kind]; }
export const __endpointMaps = { modern };