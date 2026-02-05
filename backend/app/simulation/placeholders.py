from __future__ import annotations

import math
from typing import Dict, List


PATHWAY_LABELS: Dict[str, str] = {
    "beta_adrenergic": "Beta Adrenergic Signaling",
    "calcium_handling": "Calcium Handling",
    "ion_channel_current": "Ion Channel Current",
    "contractility": "Myocyte Contractility",
}


UPSTREAM_COUPLING: Dict[str, Dict[str, float]] = {
    "calcium_handling": {"beta_adrenergic": 0.32},
    "ion_channel_current": {"beta_adrenergic": 0.18, "calcium_handling": 0.24},
    "contractility": {"calcium_handling": 0.56, "beta_adrenergic": 0.16},
}


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(value, upper))


def next_pathway_activity(
    activities: Dict[str, float],
    inhibitions: Dict[str, float],
    alpha: float = 0.22,
) -> Dict[str, float]:
    """Reduced-order deterministic transition placeholder.

    Current activity relaxes toward an inhibition-adjusted target with upstream coupling.
    This is intentionally simple for API/UI integration only.
    """

    new_activities: Dict[str, float] = {}
    for pathway, old_value in activities.items():
        direct = 1.0 - inhibitions[pathway]
        coupling_penalty = 0.0
        for upstream, weight in UPSTREAM_COUPLING.get(pathway, {}).items():
            coupling_penalty += (1.0 - activities[upstream]) * weight
        target = clamp(direct - coupling_penalty, 0.0, 1.0)
        new_activities[pathway] = clamp(old_value + alpha * (target - old_value), 0.0, 1.0)
    return new_activities


def derive_heart_rate(activities: Dict[str, float]) -> float:
    beta = activities["beta_adrenergic"]
    ion = activities["ion_channel_current"]
    calcium = activities["calcium_handling"]
    bpm = 72.0 + (beta - 0.5) * 24.0 - (1.0 - ion) * 21.0 - (1.0 - calcium) * 8.0
    return clamp(bpm, 40.0, 145.0)


def derive_pump_metric(activities: Dict[str, float]) -> float:
    contractility = activities["contractility"]
    calcium = activities["calcium_handling"]
    beta = activities["beta_adrenergic"]
    return clamp(0.58 * contractility + 0.24 * calcium + 0.18 * beta, 0.0, 1.0)


def _ecg_waveform_at(phase: float, qrs_gain: float, st_shift: float) -> float:
    p = 0.12 * math.exp(-((phase - 0.18) ** 2) / (2 * 0.028**2))
    q = -0.14 * math.exp(-((phase - 0.36) ** 2) / (2 * 0.010**2))
    r = 1.05 * qrs_gain * math.exp(-((phase - 0.40) ** 2) / (2 * 0.013**2))
    s = -0.27 * qrs_gain * math.exp(-((phase - 0.43) ** 2) / (2 * 0.012**2))
    t = 0.33 * math.exp(-((phase - 0.70) ** 2) / (2 * 0.060**2))
    return p + q + r + s + t + st_shift


def generate_ecg_trace(
    heart_rate_bpm: float,
    activities: Dict[str, float],
    phase_offset: float,
    num_points: int = 260,
) -> List[float]:
    ion = activities["ion_channel_current"]
    calcium = activities["calcium_handling"]

    qrs_gain = clamp(0.8 + 0.55 * ion, 0.6, 1.45)
    st_shift = -0.16 * (1.0 - calcium)

    cycles_per_second = heart_rate_bpm / 60.0
    values: List[float] = []
    for i in range(num_points):
        t = i / max(1, num_points - 1)
        phase = (phase_offset + cycles_per_second * t) % 1.0
        values.append(_ecg_waveform_at(phase, qrs_gain, st_shift))

    min_val = min(values)
    max_val = max(values)
    spread = max(max_val - min_val, 1e-6)
    return [((value - min_val) / spread) * 2.0 - 1.0 for value in values]
