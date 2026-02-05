"""Python conversions of ISO blocker validation workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.io import loadmat

from converted_model.matlab.dae_ode import dae_ode
from converted_model.matlab.dae_odec import dae_odec
from converted_model.matlab.dae_params import dae_params, dae_paramsc

DEFAULT_KL = 0.5
DEFAULT_KR = 10.0
DEFAULT_KG = 2.4131
DEFAULT_ALPHA_L = 1.0 / 32.0
DEFAULT_GAMMA_L = 0.3762
DEFAULT_ALPHA_A = 1.0
DEFAULT_GAMMA_A = 1.0

BLOCKER_TO_INDEX_1BASED = {
    "car": 17,
    "met": 18,
    "pro": 20,
}


@dataclass(frozen=True)
class ValidationSweepResult:
    blocker: str
    doses: np.ndarray
    camp: np.ndarray
    calcium_traces: list[np.ndarray]


def _models_dir(models_dir: str | Path | None) -> Path:
    if models_dir is None:
        return Path(__file__).resolve().parents[3] / "models"
    return Path(models_dir).resolve()


def _load_vec(path: Path, key: str) -> np.ndarray:
    data = loadmat(path)
    if key not in data:
        raise KeyError(f"Expected key `{key}` in {path}")
    return np.asarray(data[key], dtype=float).reshape(-1)


def load_calibration(models_dir: str | Path | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    base = _models_dir(models_dir)
    klcalc = _load_vec(base / "KLcalc.mat", "KLcalc")
    fitalpha = _load_vec(base / "fitalpha.mat", "fitalpha")
    ki = _load_vec(base / "Ki.mat", "Ki")
    return klcalc, fitalpha, ki


def load_y0full(models_dir: str | Path | None = None) -> np.ndarray:
    base = _models_dir(models_dir)
    data = loadmat(base / "y0full.mat")
    if "y0full" in data:
        y0 = np.asarray(data["y0full"], dtype=float).reshape(-1)
    elif "y0" in data:
        y0 = np.asarray(data["y0"], dtype=float).reshape(-1)
    else:
        raise KeyError("Could not find `y0full` or `y0` in y0full.mat")

    if y0.shape[0] != 49:
        raise ValueError(f"Expected 49-state y0full, got shape {y0.shape}")
    return y0


def _blocker_params(blocker: str, klcalc: np.ndarray, fitalpha: np.ndarray) -> tuple[float, float]:
    idx_1 = BLOCKER_TO_INDEX_1BASED[blocker.lower()]
    idx_0 = idx_1 - 1
    return float(klcalc[idx_0]), float(fitalpha[idx_0])


def create_signaling_params(
    blocker: str,
    models_dir: str | Path | None = None,
    *,
    KL: float = DEFAULT_KL,
    KR: float = DEFAULT_KR,
    KG: float = DEFAULT_KG,
    alpha_L: float = DEFAULT_ALPHA_L,
    gamma_L: float = DEFAULT_GAMMA_L,
    gamma_A: float = DEFAULT_GAMMA_A,
) -> np.ndarray:
    klcalc, fitalpha, _ = load_calibration(models_dir)
    KA, alpha_A = _blocker_params(blocker, klcalc, fitalpha)
    return dae_params(KR, KL, KA, KG, alpha_L, alpha_A, gamma_L, gamma_A)


def create_coupled_params(
    blocker: str,
    models_dir: str | Path | None = None,
    *,
    KL: float = DEFAULT_KL,
    KR: float = DEFAULT_KR,
    KG: float = DEFAULT_KG,
    alpha_L: float = DEFAULT_ALPHA_L,
    gamma_L: float = DEFAULT_GAMMA_L,
    gamma_A: float = DEFAULT_GAMMA_A,
) -> np.ndarray:
    klcalc, fitalpha, _ = load_calibration(models_dir)
    KA, alpha_A = _blocker_params(blocker, klcalc, fitalpha)
    return dae_paramsc(KR, KL, KA, KG, alpha_L, alpha_A, gamma_L, gamma_A)


def run_signaling_warmup(y0: np.ndarray, params: np.ndarray, duration_ms: float = 60.0 * 60.0 * 1000.0) -> np.ndarray:
    def rhs(t: float, state: np.ndarray) -> np.ndarray:
        dydt, _ = dae_ode(t, state, params)
        return dydt

    sol = solve_ivp(
        rhs,
        t_span=(0.0, duration_ms),
        y0=np.asarray(y0, dtype=float),
        method="BDF",
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"Signaling warmup failed: {sol.message}")
    return np.asarray(sol.y[:, -1], dtype=float)


def run_coupled_segment(
    y0: np.ndarray,
    params: np.ndarray,
    t_start: float,
    t_end: float,
    *,
    rel_tol: float = 1e-5,
    max_step: float = 2e-2,
) -> tuple[np.ndarray, np.ndarray]:
    def rhs(t: float, state: np.ndarray) -> np.ndarray:
        dydt, _ = dae_odec(t, state, params)
        return dydt

    sol = solve_ivp(
        rhs,
        t_span=(float(t_start), float(t_end)),
        y0=np.asarray(y0, dtype=float),
        method="LSODA",
        rtol=rel_tol,
        atol=1e-8,
        max_step=max_step,
    )
    if not sol.success:
        raise RuntimeError(f"Coupled segment failed: {sol.message}")
    return np.asarray(sol.t, dtype=float), np.asarray(sol.y.T, dtype=float)


def run_iso_blocker_validation(
    blocker: str,
    models_dir: str | Path | None = None,
    *,
    iso: float = 0.1,
    doses: np.ndarray | None = None,
    baseline_end_s: float = 2.0 * 60.0,
    challenge_end_s: float = 6.0 * 60.0,
) -> ValidationSweepResult:
    """Converted workflow from ISO_PRO/MET/CAR validation scripts."""
    if blocker.lower() not in BLOCKER_TO_INDEX_1BASED:
        raise ValueError(f"Unsupported blocker `{blocker}`; choose one of {sorted(BLOCKER_TO_INDEX_1BASED)}")

    if doses is None:
        doses = 10.0 ** np.arange(-5.0, 2.5 + 1e-12, 0.25)
    doses = np.asarray(doses, dtype=float)
    if baseline_end_s <= 0:
        raise ValueError("baseline_end_s must be positive")
    if challenge_end_s <= baseline_end_s:
        raise ValueError("challenge_end_s must be greater than baseline_end_s")

    y0full = load_y0full(models_dir)
    params_base = create_coupled_params(blocker, models_dir)

    camp = np.zeros_like(doses)
    calcium_traces: list[np.ndarray] = []

    for idx, dose in enumerate(doses):
        p = params_base.copy()
        y_start = y0full.copy()

        # 0-2 min baseline segment.
        p[0] = 0.0
        p[1] = 0.0
        _, y1 = run_coupled_segment(y_start, p, 0.0, baseline_end_s)

        # 2-6 min ISO + blocker challenge.
        p[0] = float(iso)
        p[1] = float(dose)
        t2, y2 = run_coupled_segment(y1[-1], p, baseline_end_s, challenge_end_s)

        y_combined = np.vstack([y1, y2])
        camp[idx] = y_combined[-1, 8]  # cAMPtot
        calcium_traces.append(y_combined[:, 46].copy())  # Cai

    return ValidationSweepResult(blocker=blocker.lower(), doses=doses, camp=camp, calcium_traces=calcium_traces)
