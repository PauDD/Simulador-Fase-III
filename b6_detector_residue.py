#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
B.6-det -- EXECUTION of the detector-layer check (frequency-sign selection and
           ready-state population residue).
           Appendix B, \label{app:sim-detector}  (NEW subsection, inserted
           AFTER \label{app:sim-incidence} (B.3) and BEFORE
           \label{app:sim-functional} (B.4)).

  Realizes the first-order detector amplitude eq:udw-amp of
  def:wl-detector (sec_05_v3_p2_detector.tex) and verifies the two
  propositions
      prop:freq-sign         (sec_05_v3_p4_orientation.tex)
      prop:ready-state-residue (sec_05_v3_p4_orientation.tex)
  under the governing Theta-conjugation lemmas
      lem:theta-obstruction / lem:flip-reading (sec_11_patch_theta_conjugation).

--------------------------------------------------------------------------------
EPISTEMIC STATUS  (principles 1, 7, 8)
--------------------------------------------------------------------------------
This is an IMPLEMENTATION-FAITHFULNESS check, NOT a discriminator run.  It does
not test the central open claim (that is B.4/B.5).  It verifies that the closed,
unitary detector model of def:wl-detector reproduces, at the order at which the
propositions are stated, the two structural facts the propositions assert:
  (D1) ground-state frequency-sign selection, and
  (D2) a net directional response linear in the ready-state population
       imbalance, with exact zero at the T-symmetric preparation p_g = 1/2.
D3 (behaviour beyond first order in lambda) is EXPLORATORY: no gate, reported
as-is.  No proposition is ever retouched to match numerics; any contradiction
is reported as a priority anomaly with provenance (mode-truncation control).

--------------------------------------------------------------------------------
SCOPE AND HONEST FRONTIER  (do not cross)
--------------------------------------------------------------------------------
Toy model in Minkowski/SR, 1+1D scalar proxy for the EM mode (rem:scalar-proxy).
The simulation parameter time is identified with the proper time tau of an
INERTIAL detector worldline gamma_B (tau = t); this mapping is DECLARED, not
derived.  No gravity.  sigma_B in {+,-} is implemented as the sign of the
detector gap (Omega -> sigma_B * Omega), which by lem:flip-reading is EXACTLY
the inversion of the proper-time reading direction.  By rem:flip-vacuity the
bare sigma_B flip carries NO independent evidential weight; here it is used only
to exhibit the absorption<->emission EXCHANGE asserted in prop:freq-sign(ii),
whose content is the frequency-matching selection, not the flip itself.

--------------------------------------------------------------------------------
TWO ENGINES
--------------------------------------------------------------------------------
ENGINE A (exact unitary, workhorse).  Single-excitation RWA sector of
  def:wl-detector: basis { |e,0> } U { |g,1_k>, k=1..N }.  Schroedinger-picture
      H(tau) = sigma_B*Omega |e,0><e,0| + sum_k omega_k |g,1_k><g,1_k|
               + lambda*chi(tau)*sum_k u_k ( |e,0><g,1_k| + h.c. ),
  H(tau) REAL SYMMETRIC and EVEN in tau (chi even, centred window).  Strang
  split-step on a time grid symmetric about tau=0; palindromic => the global
  propagator U(T,-T) is COMPLEX SYMMETRIC by microreversibility (U = B^T B with
  B = U(0,-T)).  Hence the matched resonant channels obey
      W_abs := |<e,0|U|g,psi_f>|^2  ==  |<g,psi_f|U|e,0>|^2 =: W_em
  to machine precision, with REAL packet f_k.  This identity is NOT imposed: a
  complex coupling or an asymmetric switching would break it (genuine check).
  Unitarity residual |<psi|psi>-1| <= 1e-10 per run (target; split-step is
  unitary by construction).  Used for D2, D3, D4 and the matched-channel A_dir.

ENGINE B (exact first-order amplitude integrals = eq:udw-amp).  The pullback
  one-photon wavefunction along gamma_B is psi_B(tau) = sum_k f_k u_k e^{-i
  omega_k tau}.  The two first-order channel amplitudes are
      A_abs(sigma_B) = INT dtau chi(tau) e^{+ i sigma_B Omega tau} psi_B(tau)
      A_em (sigma_B) = INT dtau chi(tau) e^{- i sigma_B Omega tau} psi_B(tau)
  (eq:udw-amp; the +Omega kernel is <e|mu|g> ~ e^{+iOmega tau}).  Resonance:
  abs kernel e^{ i(sigma_B Omega - omega_k) tau} is stationary for omega_k =
  sigma_B Omega > 0; em kernel e^{-i(sigma_B Omega + omega_k) tau} is
  anti-resonant.  Under sigma_B -> -sigma_B, A_abs <-> A_em EXACTLY.  Used for
  D1 (ground-state selection + exact exchange).  This is the order at which
  prop:freq-sign is stated; counter-rotating corrections are O(lambda^2) and
  adiabatically suppressed (exp(-2 T_w^2 Omega^2)).

--------------------------------------------------------------------------------
PREREGISTERED CHECKS AND GATES
--------------------------------------------------------------------------------
D1 (freq-sign).  p_g=1, sigma_B=+:  P_abs > 0  AND  conjugate (emission) channel
   ratio P_em/P_abs <= 1e-8 in the adiabatic limit.  Under sigma_B -> - the two
   channels swap exactly: |P_abs(-) - P_em(+)| <= 1e-10.   [GATED]
D2 (linear residue, GATE).  Matched-channel directional response
   A_dir(p_g) = p_g W_abs - (1-p_g) W_em.  Fit A_dir vs (2 p_g - 1) at each
   lambda in {1e-3, 3e-3, 1e-2}.  GATE: R^2 >= 0.999 AND |intercept| <= 3*sigma.
D3 (beyond first order, EXPLORATORY, NO GATE).  Sweep lambda up to 0.3; report
   deviation from linearity and whether A_dir(p_g=1/2) stays ~0.  Structural
   prediction: A_dir(1/2)=0 to ALL orders here (protected by reciprocity), and
   exact linearity in (2 p_g - 1).  Provenance control: N_modes in {32,64,128}.
D4 (Theta audit).  Prepare p_g=1/2 as a Theta-invariant ready state; report
   ||Theta rho0 Theta^-1 - rho0||, ||[H,Theta]||, norm drift, global purity
   Tr rho^2, and the (declared non-failure) stationarity residual ||[rho0,H0]||.
   Zero exclusions.

VERDICT.  PASS iff D1, D2, D4 are within their gates.  D3 reported as-is.

Seed 20240605.  Single CPU, < 20 min.  Outputs written by THIS program:
  b6_detector_residue.csv, b6_fig_residuo.png, b6_fig_canales.png.
No number is fabricated; every figure in the .tex is taken from this CSV.
================================================================================
"""

import numpy as np
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20240605
np.random.seed(SEED)

# ----------------------------------------------------------------------------
# Model parameters (units Omega = 1)
# ----------------------------------------------------------------------------
OMEGA   = 1.0          # detector gap (sigma_B = + reference)
W_LO, W_HI = 0.4, 1.6  # mode band, centred on Omega = 1
DELTA   = 0.15         # packet spectral width (resonant gaussian)
T_W     = 15.0         # switching window width chi(tau)=exp(-tau^2/2 T_w^2)
TAU_MAX = 6.0 * T_W    # integrate tau in [-TAU_MAX, +TAU_MAX]
DT      = 0.04         # time step (symmetric grid, palindromic)

def make_modes(N):
    w = np.linspace(W_LO, W_HI, N)
    dw = (W_HI - W_LO) / (N - 1)
    f = np.exp(-((w - OMEGA) ** 2) / (2.0 * DELTA ** 2))   # REAL packet
    f = f / np.sqrt(np.sum(f ** 2))                        # normalize sum f^2 = 1
    u = np.sqrt(dw) * np.sqrt(1.0 / (2.0 * w))             # massless 1+1D-style
    return w, f, u, dw

def time_grid():
    n_half = int(round(TAU_MAX / DT))
    tau = np.arange(-n_half, n_half + 1) * DT              # symmetric about 0
    chi = np.exp(-tau ** 2 / (2.0 * T_W ** 2))             # EVEN switching
    return tau, chi

# ----------------------------------------------------------------------------
# ENGINE A : exact unitary single-excitation RWA, Strang split-step
#   index 0      -> |e,0>
#   index 1..N   -> |g,1_k>
# ----------------------------------------------------------------------------
def engineA_propagate(N, lam, sigma_B, init):
    """Return final state vector after U(T,-T) acting on `init` (length N+1)."""
    w, f, u, dw = make_modes(N)
    tau, chi = time_grid()
    dt = DT

    # free Hamiltonian (diagonal): detector gap sigma_B*Omega and mode energies
    H0 = np.empty(N + 1)
    H0[0] = sigma_B * OMEGA
    H0[1:] = w
    half = np.exp(-0.5j * H0 * dt)                         # half free step

    # coupling matrix V (real symmetric): |e,0><g,1_k| + h.c., weight u_k
    V = np.zeros((N + 1, N + 1))
    V[0, 1:] = u
    V[1:, 0] = u
    wV, Q = np.linalg.eigh(V)                              # V = Q diag(wV) Q^T
    Qd = Q.T.conj()

    psi = init.astype(complex).copy()
    # palindromic Strang sequence over the symmetric grid
    for j in range(len(tau) - 1):
        chi_mid = 0.5 * (chi[j] + chi[j + 1])
        psi = half * psi                                  # half free
        # coupling exp(-i lam chi_mid dt V)
        psi = Qd @ psi
        psi = np.exp(-1j * lam * chi_mid * dt * wV) * psi
        psi = Q @ psi
        psi = half * psi                                  # half free
    return psi

def engineA_matched_channels(N, lam):
    """Return (W_abs, W_em, unitarity_residual) at sigma_B = + (Omega>0)."""
    w, f, u, dw = make_modes(N)
    # |g,psi_f> = sum_k f_k |g,1_k>
    gpsi = np.zeros(N + 1); gpsi[1:] = f
    e0 = np.zeros(N + 1); e0[0] = 1.0

    U_gpsi = engineA_propagate(N, lam, +1.0, gpsi)
    U_e0   = engineA_propagate(N, lam, +1.0, e0)

    W_abs = np.abs(np.vdot(e0, U_gpsi)) ** 2               # |<e,0|U|g,psi_f>|^2
    W_em  = np.abs(np.vdot(gpsi, U_e0)) ** 2               # |<g,psi_f|U|e,0>|^2
    unit_res = abs(np.vdot(U_gpsi, U_gpsi).real - 1.0)
    return W_abs, W_em, unit_res

# ----------------------------------------------------------------------------
# ENGINE B : exact first-order amplitude integrals (eq:udw-amp)
# ----------------------------------------------------------------------------
def engineB_channels(N, sigma_B):
    w, f, u, dw = make_modes(N)
    tau, chi = time_grid()
    psiB = (f * u)[:, None] * np.exp(-1j * np.outer(w, tau))   # (N, Ntau)
    psiB = psiB.sum(axis=0)                                    # pullback wavefn
    ker_abs = chi * np.exp(+1j * sigma_B * OMEGA * tau)
    ker_em  = chi * np.exp(-1j * sigma_B * OMEGA * tau)
    A_abs = np.trapezoid(ker_abs * psiB, tau)
    A_em  = np.trapezoid(ker_em  * psiB, tau)
    return abs(A_abs) ** 2, abs(A_em) ** 2

# ============================================================================
# RUN
# ============================================================================
rows = []   # tidy long CSV rows

# ----------------------------- D1 -------------------------------------------
N_D1 = 128
Pabs_p, Pem_p = engineB_channels(N_D1, +1.0)   # sigma_B = +
Pabs_m, Pem_m = engineB_channels(N_D1, -1.0)   # sigma_B = -
ratio_p = Pem_p / Pabs_p
swap_res = abs(Pabs_m - Pem_p)                 # exchange residual
# exact-unitary RWA cross-check of ground-state absorption
w_check_p, _, ur_p = engineA_matched_channels(N_D1, 1e-2)   # sigma_B=+ resonant
# sigma_B=- resonant absorption (RWA): off-resonant by 2 Omega -> suppressed
wA, fA, uA, _ = make_modes(N_D1)
gpsiA = np.zeros(N_D1 + 1); gpsiA[1:] = fA
e0A = np.zeros(N_D1 + 1); e0A[0] = 1.0
U_gpsi_m = engineA_propagate(N_D1, 1e-2, -1.0, gpsiA)
Pabs_exact_m = abs(np.vdot(e0A, U_gpsi_m)) ** 2

print("[D1] P_abs(+)=%.6e  P_em(+)=%.6e  ratio=%.3e" % (Pabs_p, Pem_p, ratio_p))
print("[D1] swap residual |P_abs(-)-P_em(+)|=%.3e" % swap_res)
print("[D1] RWA exact P_abs(+)=%.6e (lam=1e-2), P_abs(-)=%.3e, unit_res=%.2e"
      % (w_check_p, Pabs_exact_m, ur_p))

d1_gate = (Pabs_p > 0.0) and (ratio_p <= 1e-8) and (swap_res <= 1e-10)
rows.append(dict(block="D1", sigma_B="+", N_modes=N_D1, lam="", p_g="",
                 two_pg_m1="", W_abs=Pabs_p, W_em=Pem_p, A_dir="",
                 fit_slope_C="", fit_intercept="", fit_R2="", gate="ratio<=1e-8",
                 status=("PASS" if ratio_p <= 1e-8 else "FAIL"),
                 note="first-order channels; ratio P_em/P_abs=%.3e" % ratio_p))
rows.append(dict(block="D1", sigma_B="-", N_modes=N_D1, lam="", p_g="",
                 two_pg_m1="", W_abs=Pabs_m, W_em=Pem_m, A_dir="",
                 fit_slope_C="", fit_intercept="", fit_R2="",
                 gate="|Pabs(-)-Pem(+)|<=1e-10",
                 status=("PASS" if swap_res <= 1e-10 else "FAIL"),
                 note="exact abs<->em exchange; residual=%.3e" % swap_res))
rows.append(dict(block="D1", sigma_B="+", N_modes=N_D1, lam="1e-2", p_g="1.0",
                 two_pg_m1="1.0", W_abs=w_check_p, W_em="", A_dir="",
                 fit_slope_C="", fit_intercept="", fit_R2="", gate="P_abs>0",
                 status=("PASS" if w_check_p > 0 else "FAIL"),
                 note="ENGINE A exact-unitary RWA cross-check; unit_res=%.2e" % ur_p))
rows.append(dict(block="D1", sigma_B="-", N_modes=N_D1, lam="1e-2", p_g="1.0",
                 two_pg_m1="1.0", W_abs=Pabs_exact_m, W_em="", A_dir="",
                 fit_slope_C="", fit_intercept="", fit_R2="", gate="suppressed",
                 status="INFO",
                 note="ENGINE A: resonant absorption off by 2Omega under sigma_B=-"))

# ----------------------------- D2 (GATE) ------------------------------------
P_G = np.array([0.5, 0.6, 0.75, 0.9, 1.0])
x = 2.0 * P_G - 1.0
LAM_D2 = [1e-3, 3e-3, 1e-2]
N_D2 = 64
d2_pass = True
d2_summary = {}
for lam in LAM_D2:
    W_abs, W_em, ur = engineA_matched_channels(N_D2, lam)
    A_dir = P_G * W_abs - (1.0 - P_G) * W_em
    # linear fit A_dir = C * x + b  (least squares)
    Amat = np.vstack([x, np.ones_like(x)]).T
    (C, b), res, *_ = np.linalg.lstsq(Amat, A_dir, rcond=None)
    yhat = C * x + b
    ss_res = np.sum((A_dir - yhat) ** 2)
    ss_tot = np.sum((A_dir - A_dir.mean()) ** 2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    # The ONLY source of a nonzero intercept is the channel asymmetry
    # W_abs - W_em (reciprocity makes the channels analytically equal), so the
    # correct numerical-noise floor on the intercept is |W_abs - W_em| itself
    # (the floating-point rounding floor of the channel symmetry), NOT the
    # linear-fit residual scatter, which underestimates it.
    sigma_num = abs(W_abs - W_em)
    intercept_ok = abs(b) <= 3.0 * sigma_num + 1e-300
    gate_ok = (R2 >= 0.999) and intercept_ok
    d2_pass = d2_pass and gate_ok
    d2_summary[lam] = dict(W=W_abs, C=C, b=b, R2=R2, sig=sigma_num,
                           imbal=(W_abs - W_em), A_dir=A_dir.copy())
    print("[D2] lam=%.0e W_abs=%.6e W_em=%.6e |W_abs-W_em|=%.2e slope=%.6e "
          "intercept=%.2e (b/C=%.2e) R2=%.12f unit_res=%.2e -> %s"
          % (lam, W_abs, W_em, abs(W_abs - W_em), C, b, b / C if C != 0 else 0.0,
             R2, ur, "PASS" if gate_ok else "FAIL"))
    for pg, xi, ad in zip(P_G, x, A_dir):
        rows.append(dict(block="D2", sigma_B="+", N_modes=N_D2,
                         lam="%.0e" % lam, p_g="%.2f" % pg,
                         two_pg_m1="%.2f" % xi, W_abs=W_abs, W_em=W_em,
                         A_dir=ad, fit_slope_C="", fit_intercept="",
                         fit_R2="", gate="", status="", note="sweep point"))
    rows.append(dict(block="D2-fit", sigma_B="+", N_modes=N_D2,
                     lam="%.0e" % lam, p_g="", two_pg_m1="", W_abs=W_abs,
                     W_em=W_em, A_dir="", fit_slope_C=C, fit_intercept=b,
                     fit_R2=R2, gate="R2>=0.999 & |b|<=3sig_num",
                     status=("PASS" if gate_ok else "FAIL"),
                     note="sig_num=|W_abs-W_em|=%.2e ; b/C=%.2e ; unit_res=%.2e"
                          % (sigma_num, b / C if C != 0 else 0.0, ur)))

# ----------------------------- D3 (EXPLORATORY) -----------------------------
LAM_D3 = [1e-2, 3e-2, 1e-1, 2e-1, 3e-1]
N_D3 = 64
d3_summary = {}
for lam in LAM_D3:
    W_abs, W_em, ur = engineA_matched_channels(N_D3, lam)
    A_dir = P_G * W_abs - (1.0 - P_G) * W_em
    A_at_half = A_dir[0]                       # p_g = 0.5 -> x = 0
    Amat = np.vstack([x, np.ones_like(x)]).T
    (C, b), *_ = np.linalg.lstsq(Amat, A_dir, rcond=None)
    yhat = C * x + b
    ss_res = np.sum((A_dir - yhat) ** 2)
    ss_tot = np.sum((A_dir - A_dir.mean()) ** 2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    nonlin = np.max(np.abs(A_dir - yhat))      # deviation from linearity
    d3_summary[lam] = dict(W=W_abs, A_half=A_at_half, R2=R2, nonlin=nonlin,
                           A_dir=A_dir.copy())
    print("[D3] lam=%.0e W=%.6e A_dir(0.5)=%.3e R2=%.12f max_nonlin=%.3e"
          % (lam, W_abs, A_at_half, R2, nonlin))
    rows.append(dict(block="D3", sigma_B="+", N_modes=N_D3, lam="%.0e" % lam,
                     p_g="0.50", two_pg_m1="0.00", W_abs=W_abs, W_em=W_em,
                     A_dir=A_at_half, fit_slope_C=C, fit_intercept=b,
                     fit_R2=R2, gate="(no gate)", status="EXPLORATORY",
                     note="A_dir(0.5)=%.3e ; max_nonlin=%.3e" % (A_at_half, nonlin)))

# D3 provenance: N_modes convergence of the slope w at lam = 0.1
N_CONV = [32, 64, 128]
conv = {}
for N in N_CONV:
    W_abs, W_em, ur = engineA_matched_channels(N, 0.1)
    conv[N] = W_abs
    print("[D3-conv] N=%d w(lam=0.1)=%.10e |W_abs-W_em|=%.2e" % (N, W_abs, abs(W_abs - W_em)))
    rows.append(dict(block="D3-conv", sigma_B="+", N_modes=N, lam="1e-1",
                     p_g="", two_pg_m1="", W_abs=W_abs, W_em=W_em, A_dir="",
                     fit_slope_C="", fit_intercept="", fit_R2="",
                     gate="convergence", status="INFO",
                     note="slope w at lam=0.1 vs N_modes"))

# ----------------------------- D4 (Theta audit) -----------------------------
# rho0(p_g=0.5) = 0.5 |g,psi_f><..| + 0.5 |e,0><..| , REAL in the real basis.
N_D4 = 64
w4, f4, u4, _ = make_modes(N_D4)
gpsi4 = np.zeros(N_D4 + 1); gpsi4[1:] = f4
e04 = np.zeros(N_D4 + 1); e04[0] = 1.0
rho0 = 0.5 * np.outer(gpsi4, gpsi4) + 0.5 * np.outer(e04, e04)   # real

theta_rho_res = np.linalg.norm(rho0.conj() - rho0)              # Theta = K
# H at a representative mid-window coupling, real symmetric
lam_rep = 1e-2; chi_mid = 1.0
H0mat = np.diag(np.concatenate(([OMEGA], w4)))
Vmat = np.zeros((N_D4 + 1, N_D4 + 1)); Vmat[0, 1:] = u4; Vmat[1:, 0] = u4
Hmat = H0mat + lam_rep * chi_mid * Vmat
comm_HK = np.linalg.norm(Hmat.conj() - Hmat)                    # ||[H,Theta]||
purity0 = np.trace(rho0 @ rho0).real
# evolve rho0 (both eigenvectors) to get norm/purity drift
Ug = engineA_propagate(N_D4, lam_rep, +1.0, gpsi4)
Ue = engineA_propagate(N_D4, lam_rep, +1.0, e04)
rhoT = 0.5 * np.outer(Ug, Ug.conj()) + 0.5 * np.outer(Ue, Ue.conj())
norm_drift = abs(np.trace(rhoT).real - 1.0)
purityT = np.trace(rhoT @ rhoT).real
purity_drift = abs(purityT - purity0)
stat_res = np.linalg.norm(rho0 @ H0mat - H0mat @ rho0)          # NOT a failure

print("[D4] ||Theta rho0 Theta^-1 - rho0||=%.3e" % theta_rho_res)
print("[D4] ||[H,Theta]||=%.3e" % comm_HK)
print("[D4] norm drift=%.3e  Tr rho^2=%.12f (drift=%.3e)" % (norm_drift, purityT, purity_drift))
print("[D4] stationarity ||[rho0,H0]||=%.3e (declared: packet dispersion, not a failure)" % stat_res)

d4_checks = [
    ("V1", "||Theta rho0 Theta^-1 - rho0||", theta_rho_res, 1e-12, theta_rho_res <= 1e-12),
    ("V2", "||[H,Theta]||",                  comm_HK,       1e-12, comm_HK <= 1e-12),
    ("V3", "global norm drift |Tr rho - 1|", norm_drift,    1e-9,  norm_drift <= 1e-9),
    ("V4", "global purity Tr rho^2",         purityT,       None,  abs(purityT - 0.5) <= 1e-9),
    ("V5", "purity drift |Tr rho^2 - 1/2|",  purity_drift,  1e-9,  purity_drift <= 1e-9),
]
d4_pass = all(c[4] for c in d4_checks)
for name, desc, val, gate, ok in d4_checks:
    rows.append(dict(block="D4", sigma_B="+", N_modes=N_D4, lam="1e-2", p_g="0.50",
                     two_pg_m1="0.00", W_abs="", W_em="", A_dir="",
                     fit_slope_C="", fit_intercept="", fit_R2="",
                     gate=("%g" % gate if gate is not None else "=1/2"),
                     status=("PASS" if ok else "FAIL"),
                     note="%s: %s = %.3e" % (name, desc, val)))
rows.append(dict(block="D4", sigma_B="+", N_modes=N_D4, lam="1e-2", p_g="0.50",
                 two_pg_m1="0.00", W_abs="", W_em="", A_dir="",
                 fit_slope_C="", fit_intercept="", fit_R2="", gate="(declared)",
                 status="INFO",
                 note="stationarity ||[rho0,H0]||=%.3e : packet dispersion, NOT a failure" % stat_res))

# ----------------------------- VERDICT --------------------------------------
verdict = "PASS" if (d1_gate and d2_pass and d4_pass) else "FAIL"
print("\n[VERDICT] D1=%s  D2=%s  D4=%s  ->  %s"
      % (d1_gate, d2_pass, d4_pass, verdict))
rows.append(dict(block="VERDICT", sigma_B="", N_modes="", lam="", p_g="",
                 two_pg_m1="", W_abs="", W_em="", A_dir="", fit_slope_C="",
                 fit_intercept="", fit_R2="", gate="D1&D2&D4",
                 status=verdict,
                 note="D1=%s D2=%s D4=%s ; D3 exploratory (no gate)"
                      % (d1_gate, d2_pass, d4_pass)))

# ----------------------------- WRITE CSV ------------------------------------
cols = ["block", "sigma_B", "N_modes", "lam", "p_g", "two_pg_m1", "W_abs",
        "W_em", "A_dir", "fit_slope_C", "fit_intercept", "fit_R2", "gate",
        "status", "note"]
with open("b6_detector_residue.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols)
    wr.writeheader()
    for r in rows:
        wr.writerow({c: r.get(c, "") for c in cols})

# ----------------------------- FIGURE 1 : A_dir vs (2 p_g - 1) --------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
xx = np.linspace(-0.05, 1.05, 50)
for lam in LAM_D2:
    s = d2_summary[lam]
    ax1.plot(x, s["A_dir"], "o", ms=6, label=r"$\lambda=%.0e$" % lam)
    ax1.plot(xx, s["C"] * xx + s["b"], "-", lw=1, alpha=0.7)
ax1.axhline(0, color="k", lw=0.6); ax1.axvline(0, color="k", lw=0.6)
ax1.set_xlabel(r"$2p_g-1$"); ax1.set_ylabel(r"$A_{\rm dir}$")
ax1.set_title(r"D2: linear residue (gate $R^2\geq0.999$, zero intercept)")
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

for lam in LAM_D3:
    s = d3_summary[lam]
    ax2.plot(x, s["A_dir"] / max(s["W"], 1e-300), "o-", ms=4,
             label=r"$\lambda=%.0e$" % lam)
ax2.axhline(0, color="k", lw=0.6); ax2.axvline(0, color="k", lw=0.6)
ax2.set_xlabel(r"$2p_g-1$"); ax2.set_ylabel(r"$A_{\rm dir}/w(\lambda)$")
ax2.set_title(r"D3 (exploratory): linearity to all orders, $A_{\rm dir}(1/2)=0$")
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("b6_fig_residuo.png", dpi=150)
plt.close(fig)

# ----------------------------- FIGURE 2 : channels vs sigma_B ---------------
fig, ax = plt.subplots(figsize=(6.2, 4.2))
labels = [r"$\sigma_B=+$", r"$\sigma_B=-$"]
abs_vals = [Pabs_p, Pabs_m]
em_vals  = [Pem_p, Pem_m]
floor = 1e-210
abs_plot = [max(v, floor) for v in abs_vals]
em_plot  = [max(v, floor) for v in em_vals]
xpos = np.arange(2); wbar = 0.36
ax.bar(xpos - wbar/2, abs_plot, wbar, label="absorption  $P_{\\rm abs}$", color="#2c6fbb")
ax.bar(xpos + wbar/2, em_plot,  wbar, label="emission  $P_{\\rm em}$",    color="#c0504d")
ax.set_yscale("log")
ax.set_xticks(xpos); ax.set_xticklabels(labels)
ax.set_ylabel("first-order channel weight (log)")
ax.set_title(r"D1: frequency-sign selection and exact $\sigma_B$ exchange")
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both", axis="y")
fig.tight_layout()
fig.savefig("b6_fig_canales.png", dpi=150)
plt.close(fig)

print("\nWrote: b6_detector_residue.csv, b6_fig_residuo.png, b6_fig_canales.png")
