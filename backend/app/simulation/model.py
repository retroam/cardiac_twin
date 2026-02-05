from __future__ import annotations

import threading
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from converted_model.matlab.scenarios import (
    create_coupled_params,
    load_y0full,
    run_coupled_segment,
)

TRACE_POINTS = 300

INTERVENTION_KEY_MAP = {
    "iso": "Ltot",
    "ltot": "Ltot",
    "agonist": "Ltot",
    "beta_blocker": "Atot",
    "blocker": "Atot",
    "bb": "Atot",
    "atot": "Atot",
    "fsk": "FSK",
    "ibmx": "IBMX",
    # Compatibility mappings from initial placeholder keys.
    "ras_mapk": "Ltot",
    "pi3k_akt": "Atot",
    "calcium_handling": "FSK",
    "oxidative_stress": "IBMX",
}

PARAM_VECTOR_INDEX = {
    "Ltot": 0,
    "Atot": 1,
    "FSK": 2,
    "IBMX": 3,
}


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _normalize_interventions(interventions: Mapping[str, float]) -> dict[str, float]:
    normalized: dict[str, float] = {}
    for raw_key, raw_value in interventions.items():
        key = raw_key.strip().lower()
        mapped = INTERVENTION_KEY_MAP.get(key, raw_key)
        normalized[mapped] = _clamp(float(raw_value), 0.0, 1_000.0)
    return normalized


class CardioTwinModel:
    """Backend model backed by converted MATLAB coupled ODE scripts."""

    def __init__(self, *, blocker_profile: str = "car") -> None:
        self._lock = threading.Lock()
        self._blocker_profile = blocker_profile
        self._models_dir = Path(__file__).resolve().parents[3] / "models"
        self._base_params = create_coupled_params(self._blocker_profile, self._models_dir)
        self._y0full = load_y0full(self._models_dir)
        self.reset()

    def reset(self) -> None:
        with self._lock:
            self.time = 0.0
            self.step_index = 0
            self.state_vector = self._y0full.copy()
            self.active_interventions: dict[str, float] = {}
            self.electrical_trace: list[float] = []
            self.mechanical_phase = 0.0
            self.pump_metric = 0.0
            self.pathway_activity = self._extract_pathway_activity()
            self._append_trace_point()
            self._update_mechanics()

    def set_interventions(self, interventions: Mapping[str, float]) -> None:
        with self._lock:
            self.active_interventions = _normalize_interventions(interventions)

    def step(
        self,
        dt: float,
        steps: int = 1,
        transient_interventions: Mapping[str, float] | None = None,
    ) -> None:
        with self._lock:
            transient = _normalize_interventions(transient_interventions or {})

            for _ in range(steps):
                params = self._params_with_interventions(transient)
                t_next = self.time + dt
                _, y_segment = run_coupled_segment(self.state_vector, params, self.time, t_next)

                self.state_vector = y_segment[-1].copy()
                self.time = t_next
                self.step_index += 1
                self.pathway_activity = self._extract_pathway_activity()
                self._append_trace_point()
                self._update_mechanics()

    def get_state(self) -> dict[str, object]:
        with self._lock:
            return {
                "step_index": self.step_index,
                "time": self.time,
                "pathway_activity": dict(self.pathway_activity),
                "electrical_trace": list(self.electrical_trace),
                "mechanical_phase": self.mechanical_phase,
                "pump_metric": self.pump_metric,
                "active_interventions": dict(self.active_interventions),
            }

    def _params_with_interventions(self, transient: Mapping[str, float]) -> np.ndarray:
        params = self._base_params.copy()
        merged = dict(self.active_interventions)
        merged.update(transient)
        for key, value in merged.items():
            if key in PARAM_VECTOR_INDEX:
                params[PARAM_VECTOR_INDEX[key]] = float(value)
        return params

    def _extract_pathway_activity(self) -> dict[str, float]:
        # Uses converted MATLAB signaling states directly.
        return {
            "cAMPtot": float(self.state_vector[8]),
            "PKACI": float(self.state_vector[14]),
            "PLBp": float(self.state_vector[26]),
            "TnIp": float(self.state_vector[28]),
        }

    def _append_trace_point(self) -> None:
        vm = float(self.state_vector[47])
        self.electrical_trace.append(vm)
        if len(self.electrical_trace) > TRACE_POINTS:
            self.electrical_trace = self.electrical_trace[-TRACE_POINTS:]

    def _update_mechanics(self) -> None:
        cai = float(self.state_vector[46])
        ca_nsr = float(self.state_vector[42])

        # The pacing protocol in daeODEc is 1 Hz.
        self.mechanical_phase = self.time % 1.0

        # Stroke-volume proxy derived from calcium states.
        pump = 25.0 + 50_000.0 * cai + 0.30 * ca_nsr
        self.pump_metric = _clamp(pump, 0.0, 100.0)
