// Shared test utilities for server availability & base URL resolution.
const http = require('http');

function getApiBase(){
  return global.__MLIP_API_URL || process.env.MLIP_API_URL || 'http://127.0.0.1:8000';
}

let __cachedUp = false;
async function haveServer(timeoutMs=1500, attempts=4){
  if(__cachedUp) return true;
  const base = getApiBase();
  for(let i=0;i<attempts;i++){
    const ok = await new Promise(resolve=>{
      const req = http.request(base + '/serve/health', { method:'GET', agent:false }, res => {
        // Consume & immediately destroy to prevent socket reuse/keep-alive
  res.resume();
  const success = res.statusCode === 200;
  // Defer destroy to next tick without relying on setImmediate (compat)
  setTimeout(()=>{ try{req.destroy();}catch{} resolve(success); },0);
      });
      req.on('error', ()=> { try{req.destroy();}catch{} resolve(false); });
      req.setTimeout(timeoutMs, ()=>{ try{req.destroy();}catch{} resolve(false); });
      req.end();
    });
    if(ok){ __cachedUp = true; return true; }
    await new Promise(r=>setTimeout(r, 250 * (i+1))); // linear backoff
  }
  return false;
}

module.exports = { haveServer, getApiBase };
