"""Direct Python conversions of `daePARAMS.m` and `daePARAMSc.m`."""

from __future__ import annotations

import numpy as np


def _signaling_params(
    KR: float,
    KL: float,
    KA: float,
    KG: float,
    alpha_L: float,
    alpha_A: float,
    gamma_L: float,
    gamma_A: float,
) -> list[float]:
    # Drug concentrations
    Ltot = 0.0
    Atot = 0.0
    FSK = 0.0
    IBMX = 0.0

    # ETC receptor module / Gs module
    b1ARtot = 0.0132
    Gstot = 3.83
    kf_bARK = 1.1e-6
    kr_bARK = 2.2e-6
    k2_bARK = 5.5e-06
    Km_bARK = 6.26e-04
    kf_PKA = 3.6e-6
    kr_PKA = 2.2e-6
    k_G_act = 16e-3
    k_G_hyd = 0.8e-3
    k_G_reassoc = 1.21

    kf_Ra = 1.0
    kr_Ra = KR
    kf_LRa1 = 1.0
    kr_LRa1 = alpha_L * KL
    kf_LRi = 1.0
    kr_LRi = KL
    kf_LRa2 = 1.0
    kr_LRa2 = alpha_L * KR
    kf_RaG = 1.0
    kr_RaG = KG
    kf_LRaG1 = 1.0
    kr_LRaG1 = alpha_L * gamma_L * KL
    kf_LRaG2 = 1.0
    kr_LRaG2 = gamma_L * KG
    kf_ARa1 = 1.0
    kr_ARa1 = alpha_A * KA
    kf_ARi = 1.0
    kr_ARi = KA
    kf_ARa2 = 1.0
    kr_ARa2 = alpha_A * KR
    kf_ARaG1 = 1.0
    kr_ARaG1 = alpha_A * gamma_A * KA
    kf_ARaG2 = 1.0
    kr_ARaG2 = gamma_A * KG

    # AC module
    ACtot = 49.7e-3
    ATP = 5e3
    k_AC_basal = 0.2e-3
    Km_AC_basal = 1.03e3

    Kd_AC_Gsa = 0.2250
    kf_AC_Gsa = 1.0
    kr_AC_Gsa = Kd_AC_Gsa

    k_AC_Gsa = 8.5e-3
    Km_AC_Gsa = 315.0

    Kd_AC_FSK = 44.0
    k_AC_FSK = 7.3e-3
    Km_AC_FSK = 860.0

    PDEtot = 38.9e-3
    k_cAMP_PDE = 5e-3
    k_cAMP_PDEp = 2.0 * k_cAMP_PDE
    Km_PDE_cAMP = 1.3

    Kd_PDE_IBMX = 30.0
    k_PKA_PDE = 7.5e-3
    k_PP_PDE = 1.5e-3

    # PKA module
    PKAItot = 0.59
    PKAIItot = 0.025
    PKItot = 0.18
    kf_RC_cAMP = 1.0
    kf_RCcAMP_cAMP = 1.0
    kf_RcAMPcAMP_C = 4.375
    kf_PKA_PKI = 1.0
    kr_RC_cAMP = 1.64
    kr_RCcAMP_cAMP = 9.14
    kr_RcAMPcAMP_C = 1.0
    kr_PKA_PKI = 2e-4
    epsilon = 10.0

    # PP1 module
    PP1tot = 0.89
    I1tot = 0.3
    k_PKA_I1 = 60e-3
    Km_PKA_I1 = 1.0
    Vmax_PP2A_I1 = 14.0e-3
    Km_PP2A_I1 = 1.0

    Ki_PP1_I1 = 1.0e-3
    kf_PP1_I1 = 1.0
    kr_PP1_I1 = Ki_PP1_I1

    # LCC module
    LCCtot = 0.025
    PKACII_LCCtot = 0.025
    PP1_LCC = 0.025
    PP2A_LCC = 0.025
    k_PKA_LCC = 54e-3
    Km_PKA_LCC = 21.0
    k_PP1_LCC = 8.52e-3
    Km_PP1_LCC = 3.0
    k_PP2A_LCC = 10.1e-3
    Km_PP2A_LCC = 3.0

    # PLB module
    PLBtot = 106.0
    k_PKA_PLB = 54e-3
    Km_PKA_PLB = 21.0
    k_PP1_PLB = 8.5e-3
    Km_PP1_PLB = 7.0

    # PLM module
    PLMtot = 48.0
    k_PKA_PLM = 54e-3
    Km_PKA_PLM = 21.0
    k_PP1_PLM = 8.5e-3
    Km_PP1_PLM = 7.0

    # TnI module
    TnItot = 70.0
    PP2A_TnI = 0.67
    k_PKA_TnI = 54e-3
    Km_PKA_TnI = 21.0
    k_PP2A_TnI = 10.1e-3
    Km_PP2A_TnI = 4.1

    return [
        Ltot,
        Atot,
        FSK,
        IBMX,
        kf_Ra,
        kr_Ra,
        kf_LRa1,
        kr_LRa1,
        kf_LRi,
        kr_LRi,
        kf_LRa2,
        kr_LRa2,
        kf_RaG,
        kr_RaG,
        kf_LRaG1,
        kr_LRaG1,
        kf_LRaG2,
        kr_LRaG2,
        kf_ARa1,
        kr_ARa1,
        kf_ARi,
        kr_ARi,
        kf_ARa2,
        kr_ARa2,
        kf_ARaG1,
        kr_ARaG1,
        kf_ARaG2,
        kr_ARaG2,
        b1ARtot,
        Gstot,
        kf_bARK,
        kr_bARK,
        k2_bARK,
        Km_bARK,
        kf_PKA,
        kr_PKA,
        k_G_act,
        k_G_hyd,
        k_G_reassoc,
        ACtot,
        ATP,
        k_AC_basal,
        Km_AC_basal,
        kf_AC_Gsa,
        kr_AC_Gsa,
        k_AC_Gsa,
        Km_AC_Gsa,
        Kd_AC_FSK,
        k_AC_FSK,
        Km_AC_FSK,
        PDEtot,
        k_cAMP_PDE,
        k_cAMP_PDEp,
        Km_PDE_cAMP,
        Kd_PDE_IBMX,
        k_PKA_PDE,
        k_PP_PDE,
        PKAIItot,
        PKItot,
        kf_RC_cAMP,
        kf_RCcAMP_cAMP,
        kf_RcAMPcAMP_C,
        kf_PKA_PKI,
        kr_RC_cAMP,
        kr_RCcAMP_cAMP,
        kr_RcAMPcAMP_C,
        kr_PKA_PKI,
        epsilon,
        PP1tot,
        I1tot,
        k_PKA_I1,
        Km_PKA_I1,
        Vmax_PP2A_I1,
        Km_PP2A_I1,
        kf_PP1_I1,
        kr_PP1_I1,
        LCCtot,
        PKACII_LCCtot,
        PP1_LCC,
        PP2A_LCC,
        k_PKA_LCC,
        Km_PKA_LCC,
        k_PP1_LCC,
        Km_PP1_LCC,
        k_PP2A_LCC,
        Km_PP2A_LCC,
        PLBtot,
        k_PKA_PLB,
        Km_PKA_PLB,
        k_PP1_PLB,
        Km_PP1_PLB,
        PLMtot,
        k_PKA_PLM,
        Km_PKA_PLM,
        k_PP1_PLM,
        Km_PP1_PLM,
        TnItot,
        PP2A_TnI,
        k_PKA_TnI,
        Km_PKA_TnI,
        k_PP2A_TnI,
        Km_PP2A_TnI,
    ]


def dae_params(
    KR: float,
    KL: float,
    KA: float,
    KG: float,
    alpha_L: float,
    alpha_A: float,
    gamma_L: float,
    gamma_A: float,
) -> np.ndarray:
    """Converted `daePARAMS.m` parameters."""
    return np.array(
        _signaling_params(KR, KL, KA, KG, alpha_L, alpha_A, gamma_L, gamma_A),
        dtype=float,
    )


def dae_paramsc(
    KR: float,
    KL: float,
    KA: float,
    KG: float,
    alpha_L: float,
    alpha_A: float,
    gamma_L: float,
    gamma_A: float,
) -> np.ndarray:
    """Converted `daePARAMSc.m` parameters (signaling + EC coupling)."""
    base = _signaling_params(KR, KL, KA, KG, alpha_L, alpha_A, gamma_L, gamma_A)

    # EC coupling parameters.
    Vmyo = 20.8e-6
    Vnsr = 9.88e-7
    Vjsr = 9.3e-8
    ACap = 1.534e-4
    Temp = 310.0

    Nao = 140.0
    Ko = 5.4
    Cao = 1.8

    g_Na = 8.0
    g_to = 0.35
    g_ss = 0.07
    g_kibar = 0.24
    g_kp = 0.008

    f = 300.0
    g = 2e3
    gammao = 5187.5
    omega = 10.0
    pCa = 5.823e-9 * 3.0
    pK = 1.078e-11 * 3.0
    Nlcc = 3e5
    I_Ca05 = -0.458

    k_NaCa = 1483.0
    Km_Na = 87.5
    Km_Ca = 1.38
    k_sat = 0.1
    eta = 0.35
    ibarnak = 1.1
    Km_Nai = 10.0
    Km_Ko = 1.5
    ibarpca = 1.15
    Km_pca = 0.5e-3
    g_Cab = 2.8e-3
    g_Nab = 1.18e-3
    Pns = 0.0
    Km_ns = 1.2e-3

    I_upbar = 4.7
    Km_up = 3e-4
    nsrbar = 15.0
    tauon = 2e-3
    tauoff = 2e-3
    gmaxrel = 60e3
    dcaith = 0.18e-3
    Km_rel = 0.8e-3
    CSQNth = 8.75
    CSQNbar = 15.0
    Km_csqn = 0.8
    tau_tr = 5.7e-4
    TRPNbar = 0.07
    CMDNbar = 0.05
    INDObar = 0.07
    Km_trpn = 0.5128e-3
    Km_cmdn = 2.38e-3
    Km_indo = 8.44e-4

    ec = [
        Vmyo,
        Vnsr,
        Vjsr,
        ACap,
        Temp,
        Nao,
        Ko,
        Cao,
        g_Na,
        g_to,
        g_ss,
        g_kibar,
        g_kp,
        f,
        g,
        gammao,
        omega,
        pCa,
        pK,
        Nlcc,
        I_Ca05,
        k_NaCa,
        Km_Na,
        Km_Ca,
        k_sat,
        eta,
        ibarnak,
        Km_Nai,
        Km_Ko,
        ibarpca,
        Km_pca,
        g_Cab,
        g_Nab,
        Pns,
        Km_ns,
        I_upbar,
        Km_up,
        nsrbar,
        tauon,
        tauoff,
        gmaxrel,
        dcaith,
        Km_rel,
        CSQNth,
        CSQNbar,
        Km_csqn,
        tau_tr,
        TRPNbar,
        CMDNbar,
        INDObar,
        Km_trpn,
        Km_cmdn,
        Km_indo,
    ]

    return np.array(base + ec, dtype=float)
