"""Index layouts for converted MATLAB parameter/state vectors."""

from __future__ import annotations

# Shared signaling state order from daeODE/daeODEc (29 states).
SIGNALING_STATE_NAMES = (
    "Ri",
    "G",
    "b1AR_S464",
    "b1AR_S301",
    "GsaGTPtot",
    "GsaGDP",
    "Gsby",
    "AC_GsaGTP",
    "cAMPtot",
    "PDEp",
    "RC_I",
    "RCcAMP_I",
    "RCcAMPcAMP_I",
    "RcAMPcAMP_I",
    "PKACI",
    "PKACI_PKI",
    "RC_II",
    "RCcAMP_II",
    "RCcAMPcAMP_II",
    "RcAMPcAMP_II",
    "PKACII",
    "PKACII_PKI",
    "I1p_PP1",
    "I1ptot",
    "LCCap",
    "LCCbp",
    "PLBp",
    "PLMp",
    "TnIp",
)

COUPLED_STATE_NAMES = SIGNALING_STATE_NAMES + (
    "m",
    "h",
    "jo",
    "v",
    "w",
    "x",
    "yo",
    "z",
    "rto",
    "sto",
    "ssto",
    "rss",
    "sss",
    "Ca_nsr",
    "Ca_jsr",
    "Nai",
    "Ki",
    "Cai",
    "Vm",
    "trelo",
)

SIGNALING_STATE_INDEX = {name: idx for idx, name in enumerate(SIGNALING_STATE_NAMES)}
COUPLED_STATE_INDEX = {name: idx for idx, name in enumerate(COUPLED_STATE_NAMES)}

# Parameter vector indices for intervention controls.
PARAM_INDEX = {
    "Ltot": 0,
    "Atot": 1,
    "FSK": 2,
    "IBMX": 3,
}
