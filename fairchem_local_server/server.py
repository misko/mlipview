from ase.io.jsonio import encode
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .model_runtime import health_snapshot
from .models import (
    MDFromCacheIn,
    MDIn,
    MDResult,
    RelaxCalculatorName,
    RelaxFromCacheIn,
    RelaxIn,
    RelaxResult,
    SimpleFromCacheIn,
    SimpleIn,
)
from .services import (
    md_step,
    md_step_from_cache,
    relax,
    relax_from_cache,
    simple_calculate,
    simple_calculate_from_cache,
)

app = FastAPI(title="UMA-small ASE HTTP server (slim)", debug=True)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    snap = health_snapshot()
    snap["status"] = "ok"
    return snap


@app.post("/simple_calculate")
def simple_calculate_ep(inp: SimpleIn):
    return simple_calculate(inp)


@app.post("/simple_calculate_from_cache")
def simple_calculate_from_cache_ep(inp: SimpleFromCacheIn):
    return simple_calculate_from_cache(inp)


@app.post("/relax", response_model=RelaxResult)
def relax_ep(inp: RelaxIn):
    return relax(inp)


@app.post("/relax_from_cache", response_model=RelaxResult)
def relax_from_cache_ep(inp: RelaxFromCacheIn):
    return relax_from_cache(inp)


@app.post("/md", response_model=MDResult)
def md_ep(inp: MDIn):
    return md_step(inp)


@app.post("/md_from_cache", response_model=MDResult)
def md_from_cache_ep(inp: MDFromCacheIn):
    return md_step_from_cache(inp)


__all__ = [
    "app",
    "encode",
    "SimpleIn",
    "RelaxIn",
    "MDIn",
    "RelaxResult",
    "MDResult",
    "RelaxCalculatorName",
]
