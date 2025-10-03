from ase.io.jsonio import encode
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .model_runtime import health_snapshot
from .models import MDIn, MDResult, RelaxCalculatorName, RelaxIn, RelaxResult, SimpleIn
from .services import md_step, relax, simple_calculate

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


@app.post("/relax", response_model=RelaxResult)
def relax_ep(inp: RelaxIn):
    return relax(inp)


@app.post("/md", response_model=MDResult)
def md_ep(inp: MDIn):
    return md_step(inp)


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
