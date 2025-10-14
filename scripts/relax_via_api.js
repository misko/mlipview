#!/usr/bin/env node
// Simple helper: perform a /relax call via the FastAPI backend and print JSON.
// Usage: node scripts/relax_via_api.js --calc uma --steps 5 --coords "[[0,0,0],[0.96,0,0],[-0.24,0.93,0]]" --atomic-numbers "[8,1,1]"
// Defaults to a water-like geometry (O,H,H ordering) and UMA calculator.

import fs from 'fs';
import { MAX_STEP } from '../public/util/constants.js';

function parseArg(name, def){
  const idx = process.argv.indexOf('--'+name);
  if(idx===-1) return def;
  return process.argv[idx+1];
}

const baseUrl = parseArg('url','http://127.0.0.1:8000');
const calc = parseArg('calc','uma');
const steps = parseInt(parseArg('steps','5'),10);
const fmax = parseFloat(parseArg('fmax','0.05'));
const returnTrace = parseArg('return-trace','false') === 'true';
const anStr = parseArg('atomic-numbers','[8,1,1]');
const coordsStr = parseArg('coords','[[0,0,0],[0.96,0,0],[-0.24,0.93,0]]');

let atomic_numbers, coordinates;
try { atomic_numbers = JSON.parse(anStr); } catch { console.error('Bad atomic-numbers JSON'); process.exit(1); }
try { coordinates = JSON.parse(coordsStr); } catch { console.error('Bad coords JSON'); process.exit(1); }

async function main(){
  const body = { atomic_numbers, coordinates, steps, calculator: calc, fmax, return_trace: returnTrace, max_step: MAX_STEP };
  // Updated to use prefixed modern endpoint
  const resp = await fetch(baseUrl + '/serve/relax', {
    method:'POST',
    headers:{ 'Content-Type':'application/json' },
    body: JSON.stringify(body)
  });
  if(!resp.ok){
    const txt = await resp.text();
    console.error('HTTP error', resp.status, txt);
    process.exit(2);
  }
  const data = await resp.json();
  console.log(JSON.stringify(data,null,2));
}

main().catch(e=>{ console.error(e); process.exit(3); });