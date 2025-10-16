"""Thin compatibility layer for tests that import
`fairchem_local_server.server`.

Historically this module exposed FastAPI routes directly. The new
Serve-owns-HTTP entrypoint is `serve_app.py`. Some tests import
`server.health()`, so we provide that shim here delegating to
`model_runtime.health_snapshot()`.
"""

from __future__ import annotations

from .model_runtime import health_snapshot
from .models import (
    MDIn,
    MDResult,
    PrecomputedValues,
    RelaxCalculatorName,
    RelaxIn,
    RelaxResult,
    SimpleIn,
)


def health():
    snap = health_snapshot()
    snap["status"] = "ok"
    return snap


def encode(_atoms):
    """Legacy no-op encode shim.

    Older tests imported encode() from this module for its side-effects in
    legacy FastAPI server. With Ray Serve entrypoint now in serve_app.py,
    preserve a
    harmless no-op so imports continue to work without altering behavior.
    """
    return None


__all__ = [
    "health",
    "RelaxCalculatorName",
    "SimpleIn",
    "RelaxIn",
    "RelaxResult",
    "MDIn",
    "MDResult",
    "PrecomputedValues",
    "encode",
]
