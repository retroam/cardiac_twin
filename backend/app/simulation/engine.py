from __future__ import annotations

import math
from typing import Dict

from app.schemas import (
    ElectricalState,
    MechanicalState,
    PathwayState,
    SimulationState,
)
from app.simulation.model import CardioTwinModel
from app.simulation.placeholders import (
    PATHWAY_LABELS,
    clamp,
    derive_heart_rate,
    derive_pump_metric,
)


PATHWAY_PARAM_MAPPING: Dict[str, tuple[str, str, float]] = {
    "beta_adrenergic": ("Atot", "direct", 1.0),
    "calcium_handling": ("FSK", "inverse", 1.0),
    "ion_channel_current": ("IBMX", "inverse", 1.0),
    "contractility": ("Atot", "direct", 0.8),
}

ACTIVITY_BOUNDS = {
    "beta_adrenergic": (0.55, 1.0),  # cAMPtot
    "calcium_handling": (1.0, 4.0),  # PLBp
    "ion_channel_current": (0.018, 0.058),  # PKACI
    "contractility": (0.6, 2.2),  # TnIp
}


def _normalize(value: float, lower: float, upper: float) -> float:
    if upper <= lower:
        return 0.0
    return clamp((value - lower) / (upper - lower), 0.0, 1.0)


def _normalize_trace(trace: list[float]) -> list[float]:
    if not trace:
        return [0.0]
    low = min(trace)
    high = max(trace)
    spread = max(high - low, 1e-9)
    return [((value - low) / spread) * 2.0 - 1.0 for value in trace]


class ConvertedSimulationEngine:
    """Contract adapter on top of the converted MATLAB simulation core."""

    def __init__(self, seed: int = 42, *, time_acceleration: float = 5.0) -> None:
        self._seed = seed
        self._time_acceleration = time_acceleration
        self._model = CardioTwinModel()
        self._pathway_inhibitions = {pathway: 0.0 for pathway in PATHWAY_LABELS}
        self._sync_model_interventions()

    def _sync_model_interventions(self) -> None:
        intervention_values = {"Ltot": 0.0, "Atot": 0.0, "FSK": 0.0, "IBMX": 0.0}
        for pathway, inhibition in self._pathway_inhibitions.items():
            key, mode, scale = PATHWAY_PARAM_MAPPING[pathway]
            if mode == "direct":
                value = scale * inhibition
            else:
                value = scale * (1.0 - inhibition)
            intervention_values[key] += value

        self._model.set_interventions(intervention_values)

    def _derive_activities(self, state: dict[str, object]) -> Dict[str, float]:
        pathway_activity = state["pathway_activity"]
        if not isinstance(pathway_activity, dict):
            return {pathway: 0.0 for pathway in PATHWAY_LABELS}

        c_amp = float(pathway_activity.get("cAMPtot", 0.0))
        pkaci = float(pathway_activity.get("PKACI", 0.0))
        plbp = float(pathway_activity.get("PLBp", 0.0))
        tnip = float(pathway_activity.get("TnIp", 0.0))

        return {
            "beta_adrenergic": _normalize(c_amp, *ACTIVITY_BOUNDS["beta_adrenergic"]),
            "calcium_handling": _normalize(plbp, *ACTIVITY_BOUNDS["calcium_handling"]),
            "ion_channel_current": _normalize(pkaci, *ACTIVITY_BOUNDS["ion_channel_current"]),
            "contractility": _normalize(tnip, *ACTIVITY_BOUNDS["contractility"]),
        }

    def _state_from_model(self) -> SimulationState:
        model_state = self._model.get_state()
        activities = self._derive_activities(model_state)
        heart_rate = derive_heart_rate(activities)
        pump_metric = derive_pump_metric(activities)
        phase_radians = (float(model_state["mechanical_phase"]) % 1.0) * (2.0 * math.pi)
        cycle_phase = "systole" if math.sin(phase_radians) >= 0 else "diastole"

        pathways = {
            pathway: PathwayState(
                id=pathway,
                label=PATHWAY_LABELS[pathway],
                inhibition=round(self._pathway_inhibitions[pathway], 4),
                activity=round(activities[pathway], 4),
            )
            for pathway in PATHWAY_LABELS
        }

        trace = _normalize_trace([float(value) for value in model_state["electrical_trace"]])
        electrical = ElectricalState(
            heart_rate_bpm=round(heart_rate, 3),
            trace_points=[round(value, 5) for value in trace],
        )
        mechanical = MechanicalState(
            phase=cycle_phase,
            phase_radians=round(phase_radians, 5),
            pump_metric=round(pump_metric, 4),
        )

        return SimulationState(
            seed=self._seed,
            step_index=int(model_state["step_index"]),
            time_ms=int(round(float(model_state["time"]) * 1000.0)),
            pathways=pathways,
            electrical=electrical,
            mechanical=mechanical,
        )

    def get_state(self) -> SimulationState:
        return self._state_from_model()

    def reset(self, seed: int = 42) -> SimulationState:
        self._seed = seed
        self._model.reset()
        self._pathway_inhibitions = {pathway: 0.0 for pathway in PATHWAY_LABELS}
        self._sync_model_interventions()
        return self.get_state()

    def apply_intervention(self, pathway_ids: list[str], inhibition_strength: float) -> SimulationState:
        for pathway_id in pathway_ids:
            if pathway_id not in self._pathway_inhibitions:
                raise KeyError(pathway_id)
            self._pathway_inhibitions[pathway_id] = float(inhibition_strength)

        self._sync_model_interventions()
        return self.get_state()

    def simulate_step(self, steps: int, dt_ms: int) -> SimulationState:
        dt_seconds = (dt_ms / 1000.0) * self._time_acceleration
        self._model.step(dt=max(dt_seconds, 1e-4), steps=steps)
        return self.get_state()
