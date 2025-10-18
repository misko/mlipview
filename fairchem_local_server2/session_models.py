from __future__ import annotations

from enum import Enum
from typing import List, Literal, Optional

from pydantic import BaseModel, validator


class SimType(str, Enum):
    md = "md"
    relax = "relax"


class SimulationParams(BaseModel):
    # Common
    calculator: Literal["uma", "lj"] = "uma"
    # MD-specific
    temperature: float = 298.0
    timestep_fs: float = 1.0
    friction: float = 0.02
    # Relax-specific
    fmax: float = 0.05
    max_step: float = 0.2
    optimizer: Optional[str] = "bfgs"


class ClientAction(BaseModel):
    # sequence must monotonically increase from client
    seq: int
    # acknowledge last server sequence processed by client
    ack: Optional[int] = None
    type: Literal[
        "init_system",
        "update_positions",
        "start_simulation",
        "stop_simulation",
        "ping",
    ]

    # Common system fields
    atomic_numbers: Optional[List[int]] = None
    positions: Optional[List[List[float]]] = None  # N x 3
    velocities: Optional[List[List[float]]] = None  # N x 3
    cell: Optional[List[List[float]]] = None  # 3 x 3

    # start_simulation
    simulation_type: Optional[SimType] = None
    simulation_params: Optional[SimulationParams] = None
    # Correlation counters from client
    user_interaction_count: Optional[int] = None
    sim_step: Optional[int] = None

    @validator("positions")
    def _validate_pos(cls, v):  # type: ignore
        if v is None:
            return v
        if not all(len(row) == 3 for row in v):
            raise ValueError("positions must be Nx3")
        return v


class ServerResult(BaseModel):
    # increasing sequence from server to client
    seq: int
    # server can optionally mirror last seen client seq
    client_seq: Optional[int] = None

    # System state outputs
    positions: List[List[float]]
    forces: List[List[float]]
    velocities: List[List[float]]
    cell: Optional[List[List[float]]] = None

    # Optional metadata
    message: Optional[str] = None
    # Echo/correlation fields
    user_interaction_count: Optional[int] = None
    sim_step: Optional[int] = None


class SessionState(BaseModel):
    atomic_numbers: List[int] = []
    positions: List[List[float]] = []
    velocities: List[List[float]] = []
    forces: List[List[float]] = []
    cell: Optional[List[List[float]]] = None

    running: bool = False
    sim_type: Optional[SimType] = None
    params: SimulationParams = SimulationParams()

    server_seq: int = 0
    client_seq: int = 0
    client_ack: int = 0
    # Frontend correlation counters (last seen)
    user_interaction_count: int = 0
    sim_step: int = 0
