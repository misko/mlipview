"""Lightweight structured logging helper."""

from __future__ import annotations

import json
import time
from typing import Any


def log_event(event: str, **fields: Any):  # pragma: no cover (formatting)
    rec = {"ts": round(time.time(), 3), "event": event}
    rec.update(fields)
    try:
        print(json.dumps(rec, separators=(",", ":")))
    except Exception:
        print(f"[log-fallback] {event} {fields}")


__all__ = ["log_event"]
