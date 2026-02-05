from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from app.schemas import (
    InterventionRequest,
    ResetRequest,
    ServiceMeta,
    SimulationEnvelope,
    StepRequest,
)
from app.simulation.engine import ConvertedSimulationEngine

app = FastAPI(
    title="CardioTwin Simulation API",
    version="0.1.0",
    summary="Converted-model API contract for pathway->electrical->mechanical simulation",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

engine = ConvertedSimulationEngine(seed=42)
meta = ServiceMeta()


def envelope() -> SimulationEnvelope:
    return SimulationEnvelope(meta=meta, state=engine.get_state())


@app.get("/health")
def health_check() -> dict[str, str]:
    return {"status": "ok", "backend": "converted-model"}


@app.get("/state", response_model=SimulationEnvelope)
def get_state() -> SimulationEnvelope:
    return envelope()


@app.post("/simulate/step", response_model=SimulationEnvelope)
def simulate_step(request: StepRequest) -> SimulationEnvelope:
    engine.simulate_step(steps=request.steps, dt_ms=request.dt_ms)
    return envelope()


@app.post("/interventions", response_model=SimulationEnvelope)
def interventions(request: InterventionRequest) -> SimulationEnvelope:
    try:
        engine.apply_intervention(
            pathway_ids=request.pathway_ids,
            inhibition_strength=request.inhibition_strength,
        )
    except KeyError as error:
        raise HTTPException(status_code=400, detail=f"Unknown pathway id: {error.args[0]}") from error
    return envelope()


@app.post("/reset", response_model=SimulationEnvelope)
def reset(request: ResetRequest) -> SimulationEnvelope:
    engine.reset(seed=request.seed or 42)
    return envelope()
