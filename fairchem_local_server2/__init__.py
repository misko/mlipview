"""WebSocket-based UMA server (v2).

This package implements a single-connection-per-session design:
- Serve ingress exposes a single /ws endpoint maintained for the session.
- Per-connection state holds positions, velocities, forces, cell, and params.
- A pool of ASE CPU workers performs MD/relax steps, backed by UMA GPU batch.

Entrypoint: see `ws_app.deploy()`.
"""

from .ws_app import deploy  # re-export for convenience

__all__ = ["deploy"]
