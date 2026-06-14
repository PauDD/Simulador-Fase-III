#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
B.2  --  Numerical implementation and VERIFICATION of the T-symmetric matter
         sector (gradient OFF).   Appendix B, \label{app:sim-tsym}.

  Realizes  Definition  def:tsym-sector  of  sec_11_discriminator.tex.

--------------------------------------------------------------------------------
SCOPE AND HONEST FRONTIER  (do not cross; principles 1 and 8)
--------------------------------------------------------------------------------
This program implements and validates the CLOSED, UNITARY, T-SYMMETRIC matter
sector with the thermodynamic gradient switched OFF.  It establishes the
BASELINE only.

It is NOT the discriminator (protocol G1).  It does NOT:
    * flip sigma_B,
    * run the registration-substrate audit,
    * measure the directional asymmetry D(off).
Those belong to B.3 / B.4 / B.5 and remain "specified, not executed".

Therefore the numeric-vs-theory comparison below is a VERIFICATION that the
implementation faithfully realizes the sector DEFINED in def:tsym-sector.
It is NOT, and is not presented as, a confirmation of new physics.  A PASS here
means only: "the closed T-symmetric sector is implemented correctly" (unitary,
no entropy production, recurrence instead of irreversible decay).

--------------------------------------------------------------------------------
THE MODEL  (minimal standard T-symmetric model; def:tsym-sector leaves the
            concrete Hamiltonian free, so this choice is declared explicitly)
--------------------------------------------------------------------------------
A "which-path" qubit S (the matter superposition) couples to N bath spins via a
REAL-SYMMETRIC Hamiltonian:

        H  =  sigma_z^S  (x)  sum_{k=1}^{N} (g_k / 2) sigma_x^{(k)},     g_k = g0 * k.

  * H is real and symmetric  =>  with the antiunitary Theta = K (complex
    conjugation in the computational basis), K H K^{-1} = H* = H, hence
    [H, Theta] = 0   (Eq. eq:tsym-commutator).      ==> T-SYMMETRIC.
  * Evolution is CLOSED and UNITARY: psi(t) = exp(-iHt) psi(0).  No dissipator,
    no GKSL generator, no bath that carries away which-path information.
  * Initial state |+>_S (x) |0...0>_bath is REAL  =>  Theta-invariant; it is a
    pure global state (required to even speak of global unitarity / S_vN) and is
    NOT a designated low-entropy entropy-sink "ready" state of the einselection
    kind.  This is the honest reading of def:tsym-sector(iii); the subtlety
    (a pure state has zero entropy) is addressed in the summary, not hidden.
  * No geometric input is present in B.2: there is NO sigma_B and NO incidence
    operation here.  That is by design (honest frontier above).

Couplings g_k = g0 * k are COMMENSURATE (integer multiples).  This is a
deliberate verification choice: it makes the recurrence time analytically exact
(gcd of the integers is 1), so the spectral prediction can be checked to solver
precision.  The qualitative physics -- undamped revivals, no N-suppression of
the revival -- is robust to incommensurate couplings (which would give
quasi-periodic Poincare recurrence at longer, but still finite-revival, times).

GRADIENT SWITCH:
    OFF = the closed T-symmetric model above (the B.2 case).
    ON  = a minimal T-BREAKING control (qubit GKSL pure-dephasing, L=sqrt(kappa)
          sigma_z): it produces monotone decay and POSITIVE entropy production.
          Its ONLY role is to validate that the switch works (OFF: recurrence,
          zero entropy production; ON: decay, positive entropy production).
          It is NOT tuned to reproduce the P5/P8 floor V_rec ~ N^{-1.8}; that
          exponent is the IMPORTED effective benchmark of Appendix A
          (\cite{Phase2}) and is plotted only as a reference line for contrast.

--------------------------------------------------------------------------------
WHAT IS COMPUTED, AND THE THEORETICAL PREDICTION IT IS COMPARED AGAINST
--------------------------------------------------------------------------------
Reduced qubit coherence (closed form, derived from the spectrum):
        rho_01(t) = (1/2) * prod_{k=1}^{N} cos(g_k t).
  Spectral derivation: cos(g_k t) = 1/2 (e^{i g_k t}+e^{-i g_k t}); the product
  expands into frequencies  omega_eps = sum_k eps_k g_k  (eps_k = +-1), which are
  exactly the differences of the two qubit-branch spectra +-W of H.
  * Initial (Gaussian) envelope from the spectral variance Var(omega)=sum g_k^2:
        |rho_01(t)| ~ (1/2) exp(-(1/2) Gamma^2 t^2),  Gamma^2 = sum_{k=1}^N g_k^2.
    Dephasing time (drop to (1/2)/e):  t* = sqrt(2)/Gamma.
  * Full revival |rho_01| = 1/2 when every |cos(g_k t)| = 1, i.e. g_k t in pi Z;
    with g_k = g0 k (gcd=1) the first full revival is at  t_rec = pi / g0,
    INDEPENDENT of N, and the revival HEIGHT is 1/2 for ALL N (no N-suppression).

  (1) Unitarity:        ||psi(t)||, global purity Tr(rho^2), global S_vN.
                        Prediction: 1, 1, 0 (constant to solver tolerance).
  (2) T-symmetry /      net entropy production over a recurrence period and
      zero entropy      microreversibility |rho_01(t)| = |rho_01(-t)|.
      production:       Prediction: net production = 0; microreversibility exact.
  (3) Recurrence:       |rho_01(t)| numeric vs closed form; revival height at
                        t_rec; dephasing time t*; recurrence time t_rec.
                        Prediction: matches; NO monotone decay; NO N^{-1.8} floor.
  (4) N sweep:          revival height vs N (flat 1/2), recurrence time vs N
                        (flat pi/g0), dephasing time vs N (~ N^{-3/2}).

OUTPUTS: console + CSV table (numeric, theory, deviation, tolerance, PASS/FAIL
per observable per N) and four PNG figures (a)-(d).

HONESTY: the program is EXECUTED and reports REAL numbers.  If any observable
FAILS -- e.g. nonzero entropy production with the bath off, or coherence that
decays without recurrence -- it is reported as a FINDING (a gradient has crept
into the sector, the Trap-1 concern at the sector level), NOT whitewashed.
================================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------- 
# Reproducibility and global parameters
# ----------------------------------------------------------------------------- 
np.random.seed(20240601)          # fixed seed (model is deterministic; set anyway)

OUTDIR   = "/mnt/user-data/outputs"
os.makedirs(OUTDIR, exist_ok=True)

G0       = 1.0                    # base coupling; g_k = G0 * k
N_LIST   = [2, 4, 6, 8, 10]       # bath sizes (full Hilbert dim 2^(N+1) <= 2048)
N_T      = 2001                   # time-grid points; odd => t = pi lands exactly
T_MAX    = 2.0 * np.pi            # covers two recurrences (pi and 2pi)
N_REPR   = 8                      # representative N for the |rho_01|(t) figure
KAPPA_ON = 0.35                   # ON-control dephasing rate (switch validator)

# Tolerances (explicit; deviations are measured, not eyeballed)
TOL_NORM    = 1e-10               # ||psi|| - 1
TOL_PURITY  = 1e-10               # |Tr(rho^2) - 1|
TOL_GLOBS   = 1e-10               # global S_vN
TOL_MICRO   = 1e-10               # | |rho01(t)| - |rho01(-t)| |
TOL_NETS    = 1e-8                # net reduced-entropy production over a period
TOL_COH     = 1e-8               # | |rho01_num| - |rho01_theory| |
TOL_REVIVAL = 1e-6               # | revival height - 1/2 |
TOL_TREC    = 5e-3               # | t_rec_num - pi/g0 | (grid-resolution limited)
TOL_TSTAR   = 0.15               # relative error on dephasing time t*

# Pauli matrices
I2 = np.eye(2)
SX = np.array([[0.0, 1.0], [1.0, 0.0]])
SZ = np.array([[1.0, 0.0], [0.0, -1.0]])


# ----------------------------------------------------------------------------- 
# Operator construction
# ----------------------------------------------------------------------------- 
def kron_list(ops):
    out = np.array([[1.0]])
    for op in ops:
        out = np.kron(out, op)
    return out


def build_hamiltonian(N, g0=G0):
    r"""H = sigma_z^S (x) sum_k (g_k/2) sigma_x^{(k)},  g_k = g0*k.  Real symmetric."""
    dimB = 2 ** N
    W = np.zeros((dimB, dimB))            # bath operator sum_k (g_k/2) sigma_x^{(k)}
    for k in range(1, N + 1):
        g_k = g0 * k
        ops = [I2] * N
        ops[k - 1] = (g_k / 2.0) * SX
        W = W + kron_list(ops)
    H = np.kron(SZ, W)                     # qubit (x) bath
    assert np.allclose(H, H.T), "H must be symmetric"
    assert np.allclose(H, H.real), "H must be real"
    return H


def initial_state(N):
    r"""|+>_S (x) |0...0>_bath  (real => Theta-invariant, pure global state)."""
    dimB = 2 ** N
    psi = np.zeros(2 * dimB, dtype=complex)
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    psi[0 * dimB + 0] = inv_sqrt2          # qubit |0>, bath |0...0>
    psi[1 * dimB + 0] = inv_sqrt2          # qubit |1>, bath |0...0>
    return psi


# ----------------------------------------------------------------------------- 
# Closed (T-symmetric) evolution via exact diagonalization
# ----------------------------------------------------------------------------- 
def evolve_closed(N, tgrid, g0=G0):
    r"""
    Returns numeric observables on tgrid for the closed sector:
      rho01[t]  (complex),  norm[t],  purity_glob[t],  Sred[t]  (reduced qubit S_vN).
    Genuine numerical computation: build H, diagonalize, evolve state vector,
    trace out the bath numerically.
    """
    dimB = 2 ** N
    H = build_hamiltonian(N, g0)
    E, V = np.linalg.eigh(H)               # real symmetric: real E, real ortho V
    psi0 = initial_state(N)
    c = V.T @ psi0                         # expansion coefficients (real here)

    # psi(t) = V @ (exp(-iEt) * c)  vectorized over the time grid
    phases = np.exp(-1j * np.outer(tgrid, E))        # (n_t, dim)
    coeff = phases * c[None, :]                       # (n_t, dim)
    Psi = coeff @ V.T                                 # (n_t, dim)

    norm = np.sum(np.abs(Psi) ** 2, axis=1).real      # ||psi(t)||^2
    purity_glob = norm ** 2                           # Tr(rho_glob^2) for pure state

    Psi_mat = Psi.reshape(len(tgrid), 2, dimB)        # (n_t, 2, 2^N)
    # rho_S = Psi_mat @ Psi_mat^dagger  (2x2), traced over bath
    rhoS = np.einsum("tab,tcb->tac", Psi_mat, np.conjugate(Psi_mat))
    rho01 = rhoS[:, 0, 1]

    # reduced von Neumann entropy from 2x2 eigenvalues
    Sred = np.empty(len(tgrid))
    for i in range(len(tgrid)):
        ev = np.linalg.eigvalsh(rhoS[i])
        ev = np.clip(ev.real, 0.0, 1.0)
        Sred[i] = -np.sum([p * np.log(p) for p in ev if p > 1e-15])
    return rho01, norm, purity_glob, Sred


def global_entropy_exact(N, t, g0=G0):
    r"""Exact global S_vN at a single time t (small N): rho_glob = |psi><psi|."""
    H = build_hamiltonian(N, g0)
    E, V = np.linalg.eigh(H)
    c = V.T @ initial_state(N)
    psi = V @ (np.exp(-1j * E * t) * c)
    rho = np.outer(psi, np.conjugate(psi))
    ev = np.linalg.eigvalsh(rho)
    ev = np.clip(ev.real, 0.0, 1.0)
    return float(-np.sum([p * np.log(p) for p in ev if p > 1e-15]))


# ----------------------------------------------------------------------------- 
# Theory (closed form derived from the spectrum)
# ----------------------------------------------------------------------------- 
def rho01_theory(N, tgrid, g0=G0):
    out = 0.5 * np.ones_like(tgrid)
    for k in range(1, N + 1):
        out = out * np.cos(g0 * k * tgrid)
    return out


def envelope_theory(N, tgrid, g0=G0):
    Gamma2 = sum((g0 * k) ** 2 for k in range(1, N + 1))   # = g0^2 * N(N+1)(2N+1)/6
    return 0.5 * np.exp(-0.5 * Gamma2 * tgrid ** 2)


def gamma_rate(N, g0=G0):
    return np.sqrt(sum((g0 * k) ** 2 for k in range(1, N + 1)))


def t_recurrence_theory(g0=G0):
    return np.pi / g0


def t_star_theory(N, g0=G0):
    return np.sqrt(2.0) / gamma_rate(N, g0)


# ----------------------------------------------------------------------------- 
# ON control: minimal T-breaking qubit GKSL pure-dephasing (switch validator)
# ----------------------------------------------------------------------------- 
def evolve_on_control(tgrid, kappa=KAPPA_ON):
    r"""
    Lindblad pure dephasing, L = sqrt(kappa) sigma_z, H_S = 0:
        d rho/dt = kappa (sigma_z rho sigma_z - rho).
    Closed form rho01(t) = rho01(0) exp(-2 kappa t).  Integrated with solve_ivp
    for a genuine numerical solve; compared to the closed form.
    Breaks microreversibility and produces positive entropy.
    """
    rho0 = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)   # |+><+|

    def rhs(t, y):
        rho = y.reshape(2, 2)
        D = kappa * (SZ @ rho @ SZ - rho)
        return D.reshape(-1)

    sol = solve_ivp(rhs, (tgrid[0], tgrid[-1]), rho0.reshape(-1),
                    t_eval=tgrid, rtol=1e-9, atol=1e-11, method="RK45")
    rhos = sol.y.T.reshape(-1, 2, 2)
    rho01 = rhos[:, 0, 1]
    Sred = np.empty(len(tgrid))
    for i in range(len(tgrid)):
        ev = np.clip(np.linalg.eigvalsh(rhos[i]).real, 0.0, 1.0)
        Sred[i] = -np.sum([p * np.log(p) for p in ev if p > 1e-15])
    rho01_closed = 0.5 * np.exp(-2.0 * kappa * tgrid)
    return rho01, Sred, rho01_closed


# ----------------------------------------------------------------------------- 
# Extraction helpers
# ----------------------------------------------------------------------------- 
def measure_recurrence(tgrid, coh_abs):
    r"""FIRST full revival after the initial decay. The window (0.5 pi, 1.5 pi)
    isolates the first recurrence; full revivals recur at every multiple of pi
    (theory t_rec = pi/g0), so the global max over [0,2pi] is ambiguous between
    pi and 2pi -- we report the first."""
    mask = (tgrid > 0.5 * np.pi) & (tgrid < 1.5 * np.pi)
    idx = np.where(mask)[0]
    j = idx[np.argmax(coh_abs[idx])]
    return tgrid[j], coh_abs[j]


def measure_dephasing_time(tgrid, coh_abs):
    r"""First time |rho01| drops to (1/2)/e."""
    target = 0.5 / np.e
    below = np.where(coh_abs <= target)[0]
    if len(below) == 0:
        return np.nan
    return tgrid[below[0]]


def verdict(dev, tol):
    return "PASS" if dev <= tol else "FAIL"


# ----------------------------------------------------------------------------- 
# Main verification loop
# ----------------------------------------------------------------------------- 
def main():
    tgrid = np.linspace(0.0, T_MAX, N_T)
    t_pi_idx = np.argmin(np.abs(tgrid - np.pi))     # index closest to pi

    rows = []
    store = {}      # per-N stored arrays for plotting

    print("=" * 78)
    print("B.2  VERIFICATION OF THE T-SYMMETRIC MATTER SECTOR (gradient OFF)")
    print("Faithfulness check of def:tsym-sector -- NOT a test of new physics.")
    print("=" * 78)

    for N in N_LIST:
        rho01, norm, purity_glob, Sred = evolve_closed(N, tgrid)
        coh = np.abs(rho01)
        th = rho01_theory(N, tgrid)
        env = envelope_theory(N, tgrid)

        # microreversibility: |rho01(-t)| via time reversal
        rho01_neg, _, _, _ = evolve_closed(N, -tgrid)
        coh_neg = np.abs(rho01_neg)

        store[N] = dict(t=tgrid, coh=coh, th=th, env=env, Sred=Sred,
                        norm=norm, purity=purity_glob)

        # ---- (1) Unitarity ----
        dev_norm = float(np.max(np.abs(np.sqrt(norm) - 1.0)))
        dev_pur = float(np.max(np.abs(purity_glob - 1.0)))
        Sg = global_entropy_exact(N, np.pi)            # exact global S_vN at t=pi
        dev_sg = abs(Sg)

        # ---- (2) T-symmetry / zero entropy production ----
        dev_micro = float(np.max(np.abs(coh - coh_neg)))
        net_entropy_prod = float(Sred[t_pi_idx] - Sred[0])   # over one recurrence

        # ---- (3) Recurrence ----
        dev_coh = float(np.max(np.abs(coh - np.abs(th))))
        t_rec_num, revival_h = measure_recurrence(tgrid, coh)
        dev_rev = abs(revival_h - 0.5)
        dev_trec = abs(t_rec_num - t_recurrence_theory())
        t_star_num = measure_dephasing_time(tgrid, coh)
        t_star_th = t_star_theory(N)
        rel_tstar = abs(t_star_num - t_star_th) / t_star_th if t_star_num == t_star_num else np.nan

        # ---- assemble table rows ----
        rows += [
            dict(N=N, observable="unitarity_norm", numeric=1.0 - dev_norm,
                 theory=1.0, deviation=dev_norm, tol=TOL_NORM,
                 verdict=verdict(dev_norm, TOL_NORM)),
            dict(N=N, observable="unitarity_global_purity", numeric=float(purity_glob.mean()),
                 theory=1.0, deviation=dev_pur, tol=TOL_PURITY,
                 verdict=verdict(dev_pur, TOL_PURITY)),
            dict(N=N, observable="global_S_vN", numeric=Sg,
                 theory=0.0, deviation=dev_sg, tol=TOL_GLOBS,
                 verdict=verdict(dev_sg, TOL_GLOBS)),
            dict(N=N, observable="microreversibility", numeric=dev_micro,
                 theory=0.0, deviation=dev_micro, tol=TOL_MICRO,
                 verdict=verdict(dev_micro, TOL_MICRO)),
            dict(N=N, observable="net_entropy_production_per_period", numeric=net_entropy_prod,
                 theory=0.0, deviation=abs(net_entropy_prod), tol=TOL_NETS,
                 verdict=verdict(abs(net_entropy_prod), TOL_NETS)),
            dict(N=N, observable="coherence_vs_theory", numeric=dev_coh,
                 theory=0.0, deviation=dev_coh, tol=TOL_COH,
                 verdict=verdict(dev_coh, TOL_COH)),
            dict(N=N, observable="revival_height", numeric=revival_h,
                 theory=0.5, deviation=dev_rev, tol=TOL_REVIVAL,
                 verdict=verdict(dev_rev, TOL_REVIVAL)),
            dict(N=N, observable="recurrence_time", numeric=t_rec_num,
                 theory=t_recurrence_theory(), deviation=dev_trec, tol=TOL_TREC,
                 verdict=verdict(dev_trec, TOL_TREC)),
            dict(N=N, observable="dephasing_time_tstar", numeric=t_star_num,
                 theory=t_star_th, deviation=rel_tstar, tol=TOL_TSTAR,
                 verdict=verdict(rel_tstar, TOL_TSTAR)),
        ]

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTDIR, "b2_tsym_results.csv")
    df.to_csv(csv_path, index=False)

    # pretty console print
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.float_format", lambda x: f"{x:.3e}")
    print(df.to_string(index=False))
    print(f"\n[saved] {csv_path}")

    # ---- ON control (switch validator) ----
    rho01_on, Sred_on, rho01_on_closed = evolve_on_control(tgrid)
    dev_on = float(np.max(np.abs(np.abs(rho01_on) - rho01_on_closed)))
    net_entropy_on = float(Sred_on[t_pi_idx] - Sred_on[0])
    print("\n--- ON control (T-breaking GKSL dephasing; switch validator) ---")
    print(f"  closed-form match  max| |rho01_num| - 1/2 e^(-2kt) | = {dev_on:.2e}")
    print(f"  net reduced-entropy production over [0,pi]          = {net_entropy_on:.4f}  (>0 => arrow ON)")
    print(f"  contrast: OFF net entropy production over [0,pi]    = {df[df.observable=='net_entropy_production_per_period'].deviation.max():.2e}  (~0 => no arrow)")

    # =========================== FIGURES ===========================
    # (a) |rho01|(t) numeric vs theory envelope, representative N
    Nr = N_REPR
    d = store[Nr]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(d["t"], d["coh"], lw=1.4, label=r"$|\rho_{01}(t)|$ numeric")
    ax.plot(d["t"], np.abs(d["th"]), "k:", lw=1.0, label="closed-form theory")
    ax.plot(d["t"], d["env"], "r--", lw=1.2,
            label=r"Gaussian envelope $\frac{1}{2}\,e^{-\Gamma^2 t^2/2}$")
    ax.axvline(np.pi, color="green", ls="-.", lw=1.0)
    ax.annotate(r"recurrence $t_{\rm rec}=\pi/g_0$", xy=(np.pi, 0.5),
                xytext=(np.pi - 1.6, 0.42), color="green",
                arrowprops=dict(arrowstyle="->", color="green"))
    ax.set_xlabel("t"); ax.set_ylabel(r"$|\rho_{01}(t)|$")
    ax.set_title(rf"(a) Reduced coherence: dephasing $\to$ full revival (N={Nr}, gradient OFF)")
    ax.set_ylim(-0.02, 0.55); ax.legend(loc="upper center", fontsize=9)
    fig.tight_layout(); fa = os.path.join(OUTDIR, "b2_fig_a_coherence.png")
    fig.savefig(fa, dpi=140); plt.close(fig)

    # (b) reduced S_vN(t): OFF (returns to 0) vs ON (monotone to ln2); global S=0
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(tgrid, store[Nr]["Sred"], lw=1.4, label=rf"$S_{{\rm red}}(t)$ OFF (N={Nr}): returns to 0")
    ax.plot(tgrid, Sred_on, "r-", lw=1.4, label=r"$S_{\rm red}(t)$ ON control: monotone $\to \ln 2$")
    ax.axhline(np.log(2), color="grey", ls=":", lw=0.9); ax.text(0.1, np.log(2)+0.01, r"$\ln 2$", color="grey")
    ax.axhline(0.0, color="green", ls="-.", lw=0.9)
    ax.text(3.6, 0.02, r"global $S_{\rm vN}=0$ (constant, pure state)", color="green", fontsize=9)
    ax.set_xlabel("t"); ax.set_ylabel(r"$S_{\rm vN}$")
    ax.set_title("(b) Entropy: OFF non-monotone, zero net production; ON monotone production")
    ax.legend(loc="center right", fontsize=9)
    fig.tight_layout(); fb = os.path.join(OUTDIR, "b2_fig_b_entropy.png")
    fig.savefig(fb, dpi=140); plt.close(fig)

    # (c) recurrence time, revival height, and dephasing time vs N
    Ns = np.array(N_LIST)
    trec = np.array([df[(df.N == n) & (df.observable == "recurrence_time")].numeric.values[0] for n in Ns])
    revh = np.array([df[(df.N == n) & (df.observable == "revival_height")].numeric.values[0] for n in Ns])
    tstar_num = np.array([df[(df.N == n) & (df.observable == "dephasing_time_tstar")].numeric.values[0] for n in Ns])
    tstar_th = np.array([t_star_theory(n) for n in Ns])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.plot(Ns, trec, "o-", label=r"$t_{\rm rec}$ numeric")
    ax1.axhline(np.pi, color="k", ls=":", label=r"theory $\pi/g_0$ (N-independent)")
    ax1b = ax1.twinx()
    ax1b.plot(Ns, revh, "s--", color="darkorange", label="revival height")
    ax1b.axhline(0.5, color="darkorange", ls=":", lw=0.8)
    ax1b.set_ylabel("revival height", color="darkorange"); ax1b.set_ylim(0, 0.6)
    ax1.set_xlabel("N"); ax1.set_ylabel(r"$t_{\rm rec}$")
    ax1.set_title("(c1) recurrence time & height vs N (both flat)")
    ax1.legend(loc="center left", fontsize=8)

    ax2.loglog(Ns, tstar_num, "o", label=r"$t_*$ numeric")
    ax2.loglog(Ns, tstar_th, "k-", label=r"$t_*=\sqrt{2}/\Gamma,\ \Gamma^2=\sum g_k^2$")
    ax2.set_xlabel("N"); ax2.set_ylabel(r"dephasing time $t_*$")
    ax2.set_title(r"(c2) dephasing time $\sim N^{-3/2}$")
    ax2.legend(fontsize=9)
    fig.tight_layout(); fc = os.path.join(OUTDIR, "b2_fig_c_scaling.png")
    fig.savefig(fc, dpi=140); plt.close(fig)

    # (d) contrast: OFF revival height (flat) vs P5 effective benchmark N^{-1.8}
    fig, ax = plt.subplots(figsize=(8, 4.7))
    ax.plot(Ns, revh, "o-", color="C0", lw=1.6,
            label="OFF: revival height (this work, real data)")
    A = 0.5 * (2 ** 1.8)
    p5 = A * Ns.astype(float) ** (-1.8)
    ax.plot(Ns, p5, "r^--", lw=1.4,
            label=r"P5/P8 effective floor $\sim N^{-1.8}$ (imported, Phase II)")
    on_floor = 0.5 * np.exp(-2.0 * KAPPA_ON * np.pi)
    ax.axhline(on_floor, color="purple", ls=":", lw=1.2,
               label=rf"ON control floor at $t=\pi$ (={on_floor:.3f}, single qubit)")
    ax.set_yscale("log"); ax.set_xlabel("N")
    ax.set_ylabel(r"residual coherence at $t=\pi$")
    ax.set_title("(d) Contrast: OFF undamped revival vs ON / effective N-suppression")
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout(); fd = os.path.join(OUTDIR, "b2_fig_d_contrast.png")
    fig.savefig(fd, dpi=140); plt.close(fig)

    # =========================== SUMMARY ===========================
    overall = "PASS" if (df.verdict == "PASS").all() else "FAIL"
    n_fail = int((df.verdict == "FAIL").sum())
    print("\n" + "=" * 78)
    print("SUMMARY (per-observable verdict; deviations measured, not eyeballed)")
    print("=" * 78)
    for obs in df.observable.unique():
        sub = df[df.observable == obs]
        worst = sub.loc[sub.deviation.idxmax()]
        print(f"  {obs:38s} worst dev = {worst.deviation:.2e} (N={int(worst.N)})  "
              f"tol={worst.tol:.1e}  [{'PASS' if (sub.verdict=='PASS').all() else 'FAIL'}]")
    print("-" * 78)
    print(f"  OVERALL: {overall}   ({n_fail} failing rows of {len(df)})")
    print("  Interpretation (honest frontier): a PASS verifies that the CLOSED,")
    print("  UNITARY, T-SYMMETRIC sector is implemented faithfully (def:tsym-sector):")
    print("  global unitarity exact, zero net entropy production with the bath OFF,")
    print("  microreversibility, and recurrence (NOT irreversible decay; NO N^{-1.8}")
    print("  floor). This is a faithfulness check, NOT evidence of new physics.")
    print("  The discriminator D(off), the sigma_B flip and the registration audit")
    print("  are NOT performed here; they remain specified-not-executed (B.3-B.5).")
    print("=" * 78)

    print("\n[figures saved]")
    for f in (fa, fb, fc, fd):
        print("  ", f)

    return df


if __name__ == "__main__":
    main()
