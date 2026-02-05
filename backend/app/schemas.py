from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


class PathwayState(BaseModel):
    id: str
    label: str
    inhibition: float = Field(ge=0.0, le=1.0)
    activity: float = Field(ge=0.0, le=1.0)


class ElectricalState(BaseModel):
    heart_rate_bpm: float = Field(ge=20.0, le=220.0)
    trace_points: List[float]


class MechanicalState(BaseModel):
    phase: str
    phase_radians: float
    pump_metric: float = Field(ge=0.0, le=1.0)


class SimulationState(BaseModel):
    seed: int
    step_index: int = Field(ge=0)
    time_ms: int = Field(ge=0)
    pathways: Dict[str, PathwayState]
    electrical: ElectricalState
    mechanical: MechanicalState


class ServiceMeta(BaseModel):
    placeholder_backend: bool = False
    model_name: str = "CardioTwin Converted MATLAB Dynamics"
    notes: str = (
        "Frontend contract adapter backed by converted MATLAB signaling/electrical kernels."
    )


class SimulationEnvelope(BaseModel):
    meta: ServiceMeta
    state: SimulationState


class StepRequest(BaseModel):
    steps: int = Field(default=1, ge=1, le=200)
    dt_ms: int = Field(default=60, ge=10, le=1000)


class InterventionRequest(BaseModel):
    pathway_ids: List[str] = Field(min_length=1)
    inhibition_strength: float = Field(ge=0.0, le=1.0)
    step_index: Optional[int] = Field(default=None, ge=0)

    @field_validator("pathway_ids")
    @classmethod
    def non_empty_ids(cls, value: List[str]) -> List[str]:
        cleaned = [item.strip() for item in value if item.strip()]
        if not cleaned:
            raise ValueError("pathway_ids must contain at least one pathway id")
        return cleaned


class ResetRequest(BaseModel):
    seed: Optional[int] = 42
