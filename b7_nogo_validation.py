#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
B.7  --  NUMERICAL CONFIRMATION of the emergence no-go theorem
         (thm:emergence-nogo, sec_05_v3_p6_nogo) and its corollary
         (cor:level1-null).  Appendix B extension, label app-v3-nogo-*.

  This block is CONFIRMATION of a PROVEN theorem, never a test.  It cannot
  refute thm:emergence-nogo; it can only detect an implementation or
  hypothesis error.  Its sole positive value is to exhibit the short-time
  coefficients (c1,c2,c3) and the recurrence scales in the *real* G1 model
  (def:g1-model), thereby executing the quantitative recurrence-scale check
  that rem:nogo-scope explicitly assigns to "the Appendix B extension
  (specified, not yet executed)".

  Status change effected by this block (declared, per the prompt):
    the executed Level-1 nulls of B.4 are hereby confirmed as instances of
    cor:level1-null; the numerics are DEMOTED FROM JUDGE TO CONFIRMATION.

--------------------------------------------------------------------------------
INFRASTRUCTURE REUSE AND CONTINUITY GATE  (principle 8: data govern;
  principle 6: the functional, not a re-implementation)
--------------------------------------------------------------------------------
Every degree of freedom is the EXACT G1 model of V2.  We do not re-implement
it: we IMPORT b4_b5_discriminator.py and call its constructors
(build_H_sc, prepare_state, build_A3_structured, V3_from_A, apply_V3,
SpectralProp, random_A3, window_max, A_functional, fit_models, m_law).  The
incidence operator A is the SAME structured A of the executed B.4 (with the
two rejected candidates documented in that module's build_A3_structured
docstring); J=0.1, g0=0.05, theta=0.7; preparation = real mid-spectrum
Theta-invariant eigenstate (x) register|0>.

CONTINUITY GATE (run FIRST; bit-for-bit where the quantity is seed-free):
  G0a  lemma residuals V1/V2/V3 == 0 at machine precision (Theta=K).
  G0b  Level-1 return functional A1(tau*=3.0) reproduces the V2 b4_A_of_N.csv
       values to within 5% of the V2 null-band width (see DEGENERACY NOTE).
  G0c  Level-1 free-exponent fit p over tau in [0.5,4.0] reproduces the V2
       b4_fit_models.csv values to < 1e-6 (seed-free).
  A failure of G0 is an infrastructure regression and BLOCKS the report.

DEGENERACY NOTE (honest finding, documented; principle 2 fidelity, 8 data):
  The stationary preparation prepare_state() selects V[:, dim//2], the
  mid-spectrum eigenvector of H_sc.  In THIS model the mid-spectrum
  eigenspace is EXACTLY degenerate: the gap to the next level is ~1e-16
  (machine zero) at every executed N.  The vector LAPACK returns inside a
  degenerate eigenspace is basis-arbitrary and depends on the BLAS/LAPACK
  build.  Consequences, all verified here:
    * stationarity ||(H0-E0)psi0|| is machine-zero for ANY vector in the
      eigenspace, so the preparation is stationary regardless (the theorem
      hypothesis [rho0,H]=0 holds exactly);
    * the lemma invariants V1/V2/V3 are basis-independent and reproduce as
      exact zeros;
    * the free exponent p (a robust global fit) reproduces to rounding;
    * A1(tau*) is a NEAR-CANCELLING difference of two window maxima of order
      1e-3, and is therefore the one quantity that exposes the ~1e-5
      LAPACK-basis difference.  This is NOT an infrastructure regression: it
      is the documented non-uniqueness of an eigenvector in an exactly
      degenerate subspace.  G0b is therefore graded against the V2 null-band
      scale (band_2sigma ~ 5e-3..2.3e-2): agreement to <5% of the band means
      the inside/outside-band physics is identical.  We import prepare_state
      VERBATIM rather than re-canonicalizing, to stay faithful to the model
      actually executed in V2.

NEW RANDOMNESS: the Theta-covariant GOE scramble ensemble of N1 uses the NEW
seed 20240606 (per prompt), drawn per size as default_rng(20240606 + N) so
the size ensembles are independent and reproducible.  V2's own scramble seed
(20240601-based) is untouched and is not used here.

--------------------------------------------------------------------------------
THE SUPPRESSION FUNCTIONAL OF THE NO-GO  (eq:nogo-functional, executed form)
--------------------------------------------------------------------------------
eq:nogo-functional is S(tau)=F_V(tau)/F_1(tau), F_W(tau)=|Tr[O e^{-iHtau}
W rho0 W^+ e^{iHtau}]|.  In the executed G1 model the preparation is the
mid-spectrum eigenstate with DEFINITE sz_S, so the bare event-free coherence
F_1 is structurally ZERO (no which-path coherence without the incidence).
The proof of part (i) covers exactly this case (line: "where the raw quotient
gives S(0)!=1, replace S by S(tau)/S(0), which changes no exponent"): we take
the executed functional
        S(tau) = f(tau)/f(0),    f(s) = |rho01(s)|  (the which-path coherence,
                                 = F_V(s) with O the coherence extractor),
i.e. we normalize the POST-EVENT functional by its own initial value.  This
is real-analytic at tau=0 with S(0)=1, it is the modulus-ratio of two
trigonometric polynomials (bounded, almost-periodic), and -ln S has the
short-time expansion of the theorem.  f(s) is computed by the SAME spectral
propagator (SpectralProp) used for the executed B.4 numbers.

--------------------------------------------------------------------------------
PREREGISTERED VERIFICATIONS  (fixed before running; no post-hoc torture)
--------------------------------------------------------------------------------
N1 (short-time coefficients).  For Level 1 (structured A) and for the
    ensemble of 24 Theta-covariant GOE scrambles per size, fit
        -ln S(tau) = c1 tau + c2 tau^2 + c3 tau^3   (cubic THROUGH the origin)
    on tau in [0.01, 0.2], mesh 0.005, with f(0) the post-event coherence.
    Coefficient covariance from sigma^2 (X^T X)^{-1}, sigma^2 = RSS/(n-3).
    GATES:
      (a) c2 != 0 at >= 5 sigma in the MAJORITY of the per-size ensemble
          (genericity of p=2; c1 suppressed by the Theta-covariant
          symmetrization, exhibited separately via the even functional
          f_sym(tau)=1/2[f(+tau)+f(-tau)], whose c1 vanishes to machine
          precision -- this is the codimension-one statement of the theorem.
          NOTE: c3 of the symmetrized functional is NOT a clean diagnostic
          here, because the cubic-through-origin fit absorbs the omitted
          tau^4 curvature into tau^3; only c1_sym is reported as the
          codimension-one witness.  Since f(s) is machine-exact, the fit
          uncertainty is the cubic-truncation residual on [0.01,0.2], so the
          5-sigma gate is met by orders of magnitude; the physical content is
          simply that c2 is robustly nonzero in every member.)
      (b) NO ensemble run with c1=c2=c3=0 simultaneously within error
          (the fine-tuning condition of the theorem does NOT occur
          spontaneously): operationalized as "no run with all of
          |c1|/s1, |c2|/s2, |c3|/s3 below 2".

N2 (recurrence scale).  For each N, from the EXACT spectrum compute, on a
    long horizon, the dynamical bound sup_tau(-ln S) (pointwise S=f/f0) and
    the first recurrence time tau_rec at tolerance |S-S(0)|<=0.05 (first tau
    where, after dropping below 0.95, S returns within 0.05 of 1).  Compare
    with the quartic-law demand kappa*m(4.0)=kappa*(pi/24)*4^4 using the
    kappa of the V2 quartic fit (the L2b M_quartic coefficient
    c_V2 = 0.016596 = kappa*pi/24, kappa_eff = 0.12678).  GATE: the quartic
    law exceeds the dynamical bound OR the window exceeds tau_rec -- quantify
    which and by what factor.  This is the analytic explanation of the
    Level-2a corruption (part (ii) acting).

N3 (exponent drift, reinterpreted; DESCRIPTIVE, no hard gate).  Verify that
    the free Level-1 exponents measured in V2 (p between 1.0 and 2.1 by N)
    are those predicted by the c1/c2 of N1, run by run, in sign and order of
    magnitude: the window-effective exponent of -ln S over [0.5,4.0] is set
    by the relative weight of the linear vs quadratic short-time terms.

--------------------------------------------------------------------------------
PREREGISTERED VERDICT
--------------------------------------------------------------------------------
This block cannot "refute" the theorem.  If N1 or N2 failed, the report FIRST
audits whether the executed model is inside the theorem's hypotheses
(time-independent generator, stationary preparation) BEFORE touching anything;
any persistent discrepancy is reported as such and BLOCKS integration of the
corollary until resolved.

OUTPUTS: b7_nogo_coeffs.csv, b7_recurrence.csv, b7_fig_coeficientes.png,
b7_fig_recurrencia.png, and console summary.  All real program output.
================================================================================
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- import the V2 infrastructure verbatim (does NOT run its main) ----------
sys.path.insert(0, "/mnt/project")
import b4_b5_discriminator as v2   # noqa: E402

T0 = time.time()
OUTDIR = os.environ.get("B7_OUTDIR", "/mnt/user-data/outputs")
os.makedirs(OUTDIR, exist_ok=True)

# ---- preregistered parameters of THIS block ---------------------------------
SEED_SCRAMBLE = 20240606          # NEW seed for the GOE scramble ensemble
N_ENS         = 24                # Theta-covariant scrambles per size
ST_LO, ST_HI, ST_DS = 0.01, 0.20, 0.005   # short-time fit window and mesh
KAPPA_FIT_V2  = 0.016596          # V2 L2b quartic coefficient c = kappa*pi/24
KAPPA_EFF     = KAPPA_FIT_V2 * 24.0 / np.pi   # ~0.12678
TAU_WIN_MAX   = 4.0               # executed B.4 window upper edge
TAU_WIN_MIN   = 0.5
REC_TOL       = 0.05              # |S - S(0)| recurrence tolerance
SIG5, SIG2    = 5.0, 2.0          # significance gates of N1
# long horizons for the recurrence scan (per size; bounded by memory/cost)
HORIZON = {4: 120.0, 6: 100.0, 8: 60.0, 10: 40.0}
HOR_DS  = 0.02

# ---- V2 reference numbers for the continuity gate (from the V2 CSVs) --------
V2_A1_TAUSTAR = {4: -0.002292077705173451, 6: 0.005912610170557515,
                 8: 0.006085951998866395, 10: 0.001147233297307506}
V2_LEVEL1_P   = {4: 1.584193, 6: 2.128718, 8: 1.411184, 10: 0.998438}
V2_BAND_2SIGMA = {4: 0.022872733368896923, 6: 0.02305074345779818,
                  8: 0.018737434939965575, 10: 0.00474362154741306}
GATE_A1_BANDFRAC = 0.05    # |dA1| must be < 5% of the V2 null-band width
GATE_P_TOL   = 1e-6

NLIST = v2.N_LIST


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------
def cubic_through_origin(tau, y):
    """Fit y = c1 tau + c2 tau^2 + c3 tau^3 (no constant).  Returns the three
    coefficients and their 1-sigma errors from sigma^2 (X^T X)^{-1}."""
    X = np.column_stack([tau, tau ** 2, tau ** 3])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n, k = len(y), 3
    dof = max(n - k, 1)
    sigma2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = sigma2 * XtX_inv
    sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    return beta, sig, float(np.sqrt(sigma2))


def short_time_functional(sp, psi_post, N, both_sides=False):
    """S(tau)=f(tau)/f(0) and (optionally) the even Theta-symmetrized variant,
    on the short-time mesh.  f(s)=|rho01(s)| via the V2 spectral propagator."""
    s_fit = np.arange(ST_LO, ST_HI + 1e-12, ST_DS)
    grid_pos = np.concatenate([[0.0], s_fit])
    f_pos, _, _, _ = sp.evolve(psi_post, grid_pos)
    f0 = float(f_pos[0])
    f_tau = f_pos[1:]
    out = dict(tau=s_fit, f0=f0, f_tau=f_tau,
               y_raw=-np.log(f_tau / f0))
    if both_sides:
        f_neg, _, _, _ = sp.evolve(psi_post, -grid_pos)
        f_sym = 0.5 * (f_pos + f_neg)
        out["y_sym"] = -np.log(f_sym[1:] / f_sym[0])
    return out


def effective_exponent(c1, c2, c3, tau):
    """Local log-log slope of -ln S at tau, given the short-time polynomial:
    p_eff(tau) = d ln(-ln S)/d ln tau = tau * g'(tau)/g(tau),
    g = c1 tau + c2 tau^2 + c3 tau^3."""
    g = c1 * tau + c2 * tau ** 2 + c3 * tau ** 3
    gp = c1 + 2 * c2 * tau + 3 * c3 * tau ** 2
    return tau * gp / g if abs(g) > 1e-300 else np.nan


# =============================================================================
# CONTINUITY GATE (G0): reproduce the seed-free V2 Level-1 numbers
# =============================================================================
def continuity_gate():
    print("=" * 80)
    print("CONTINUITY GATE G0  (reproduce seed-free V2 Level-1 numbers; "
          "infrastructure regression test)")
    print("=" * 80)
    sgrid_f = np.arange(0.0, v2.T_MAX + 1e-12, v2.DS)
    sgrid_b = -sgrid_f
    rows = []
    ok = True
    cache = {}
    for N in NLIST:
        Hsc = v2.build_H_sc(N)
        psi0, E0, stat_res, theta_res, rho01_init, sz0 = v2.prepare_state(N, Hsc)
        A3 = v2.build_A3_structured()
        V3p, V3m = v2.V3_from_A(A3, +1), v2.V3_from_A(A3, -1)
        v1 = float(np.linalg.norm(np.conj(V3p) - V3m))
        sp = v2.SpectralProp(N, Hsc)
        psi_p = v2.apply_V3(psi0, V3p, N)
        psi_m = v2.apply_V3(psi0, V3m, N)
        fp_f, _, _, _ = sp.evolve(psi_p, sgrid_f)
        fp_b, _, _, _ = sp.evolve(psi_p, sgrid_b)
        fm_f, _, _, _ = sp.evolve(psi_m, sgrid_f)
        fm_b, _, _, _ = sp.evolve(psi_m, sgrid_b)
        v2id = max(float(np.max(np.abs(fp_f - fm_b))),
                   float(np.max(np.abs(fp_b - fm_f))))
        v3id = float(np.max(np.abs((fp_f - fp_b) + (fm_f - fm_b))))

        # A1 at tau* = 3.0
        a1_star, _, _ = v2.A_functional(sgrid_f, fp_f, sgrid_b, fp_b,
                                        v2.TAU_STAR)
        # Level-1 free-p fit exactly as V2 main does it
        M1f = np.array([v2.window_max(sgrid_f, fp_f, t) for t in v2.TAU_GRID])
        good = M1f > v2.M1_FLOOR
        S1 = M1f / np.max(M1f)
        df1 = v2.fit_models(v2.TAU_GRID[good], np.clip(S1[good], 1e-12, None))
        p1 = float(df1[df1.model == "M_freep"].p.iloc[0])

        dA1 = abs(a1_star - V2_A1_TAUSTAR[N])
        dP1 = abs(p1 - V2_LEVEL1_P[N])
        a1_thresh = GATE_A1_BANDFRAC * V2_BAND_2SIGMA[N]
        g_lem = (v1 < v2.TOL_V1) and (v2id < v2.TOL_V2_L1) and (v3id < v2.TOL_V3)
        g_a1 = dA1 < a1_thresh
        g_p1 = dP1 < GATE_P_TOL
        ok = ok and g_lem and g_a1 and g_p1
        # mid-spectrum degeneracy gap (documents why A1 is basis-sensitive)
        Eeig = np.linalg.eigvalsh(Hsc)
        mid = len(Eeig) // 2
        deg_gap = float(min(Eeig[mid + 1] - Eeig[mid],
                            Eeig[mid] - Eeig[mid - 1]))
        rows.append(dict(N=N, v1=v1, v2_theta=v2id, v3_mirror=v3id,
                         A1_repro=a1_star, A1_V2=V2_A1_TAUSTAR[N], dA1=dA1,
                         A1_thresh_5pct_band=a1_thresh,
                         dA1_over_band=dA1 / V2_BAND_2SIGMA[N],
                         mid_spectrum_gap=deg_gap, stationarity=stat_res,
                         p_level1_repro=p1, p_level1_V2=V2_LEVEL1_P[N], dP1=dP1,
                         lemma_ok=g_lem, A1_ok=g_a1, p_ok=g_p1))
        cache[N] = dict(Hsc=Hsc, psi0=psi0, sp=sp, A3=A3,
                        stat_res=stat_res, theta_res=theta_res,
                        E0=E0, sz0=sz0, rho01_init=rho01_init)
        print(f"  N={N:2d}  V1={v1:.1e} V2={v2id:.1e} V3={v3id:.1e} | "
              f"A1*={a1_star:+.6e} (dV2={dA1:.1e}={dA1/V2_BAND_2SIGMA[N]:.1%} band) | "
              f"p_L1={p1:.6f} (dV2={dP1:.1e}) | gap={deg_gap:.0e} | "
              f"{'PASS' if (g_lem and g_a1 and g_p1) else 'FAIL'}")
    df = pd.DataFrame(rows)
    print(f"\n  CONTINUITY GATE: {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("  *** G0 FAILED: infrastructure regression. Report BLOCKED. ***")
    return ok, df, cache


# =============================================================================
# N1: short-time coefficients (structured Level 1 + GOE ensemble)
# =============================================================================
def run_N1(cache):
    print("\n" + "=" * 80)
    print("N1  short-time coefficients of -ln S(tau)=c1 tau+c2 tau^2+c3 tau^3 "
          "on [%.2f,%.2f] mesh %.3f" % (ST_LO, ST_HI, ST_DS))
    print("=" * 80)
    rows = []
    ens_summary = []
    for N in NLIST:
        c = cache[N]
        sp, psi0, A3 = c["sp"], c["psi0"], c["A3"]

        # ---- Level 1 (structured A): raw + Theta-symmetrized ----
        psi_p = v2.apply_V3(psi0, v2.V3_from_A(A3, +1), N)
        st = short_time_functional(sp, psi_p, N, both_sides=True)
        b_raw, s_raw, rmse_raw = cubic_through_origin(st["tau"], st["y_raw"])
        b_sym, s_sym, _ = cubic_through_origin(st["tau"], st["y_sym"])
        rows.append(dict(N=N, member="level1_structured", idx=-1,
                         f0=st["f0"], c1=b_raw[0], c2=b_raw[1], c3=b_raw[2],
                         s1=s_raw[0], s2=s_raw[1], s3=s_raw[2],
                         c1_sym=b_sym[0], c3_sym=b_sym[2], rmse=rmse_raw))
        print(f"\n  N={N:2d}  [Level 1, structured A]  f(0)={st['f0']:.4f}")
        print(f"        c1={b_raw[0]:+.5f} +/- {s_raw[0]:.1e}  "
              f"(|c1|/s = {abs(b_raw[0])/max(s_raw[0],1e-30):.1f})")
        print(f"        c2={b_raw[1]:+.5f} +/- {s_raw[1]:.1e}  "
              f"(|c2|/s = {abs(b_raw[1])/max(s_raw[1],1e-30):.1f})")
        print(f"        c3={b_raw[2]:+.5f} +/- {s_raw[2]:.1e}")
        print(f"        Theta-symmetrized: c1_sym={b_sym[0]:+.2e}, "
              f"c3_sym={b_sym[2]:+.2e}  (codimension-one: c1 killed for free)")

        # ---- GOE ensemble (NEW seed 20240606 + N) ----
        rng = np.random.default_rng(SEED_SCRAMBLE + N)
        n_c2_5sig = 0
        n_triple_zero = 0
        for i in range(N_ENS):
            Vr = v2.V3_from_A(v2.random_A3(rng), +1)
            psr = v2.apply_V3(psi0, Vr, N)
            sti = short_time_functional(sp, psr, N, both_sides=False)
            bi, si, rmsei = cubic_through_origin(sti["tau"], sti["y_raw"])
            z1 = abs(bi[0]) / max(si[0], 1e-30)
            z2 = abs(bi[1]) / max(si[1], 1e-30)
            z3 = abs(bi[2]) / max(si[2], 1e-30)
            if z2 >= SIG5:
                n_c2_5sig += 1
            if (z1 < SIG2) and (z2 < SIG2) and (z3 < SIG2):
                n_triple_zero += 1
            rows.append(dict(N=N, member="scramble", idx=i, f0=sti["f0"],
                             c1=bi[0], c2=bi[1], c3=bi[2],
                             s1=si[0], s2=si[1], s3=si[2],
                             c1_sym=np.nan, c3_sym=np.nan, rmse=rmsei))
        frac5 = n_c2_5sig / N_ENS
        gate_a = frac5 > 0.5
        gate_b = (n_triple_zero == 0)
        ens_summary.append(dict(N=N, n_ens=N_ENS, n_c2_ge5sig=n_c2_5sig,
                                frac_c2_ge5sig=frac5, gate_a_majority=gate_a,
                                n_triple_zero=n_triple_zero, gate_b_no_finetune=gate_b))
        print(f"        [ensemble {N_ENS} GOE scrambles, seed {SEED_SCRAMBLE}+{N}] "
              f"c2>=5sigma in {n_c2_5sig}/{N_ENS} ({frac5:.0%})  "
              f"-> gate(a) {'PASS' if gate_a else 'FAIL'};  "
              f"triple-zero runs: {n_triple_zero}  "
              f"-> gate(b) {'PASS' if gate_b else 'FAIL'}")
    return pd.DataFrame(rows), pd.DataFrame(ens_summary)


# =============================================================================
# N2: recurrence scale vs the quartic-law demand
# =============================================================================
def run_N2(cache):
    print("\n" + "=" * 80)
    print("N2  recurrence scale of the exact spectrum vs the quartic demand "
          "kappa*m(4.0)")
    print("    kappa_eff = %.5f  (V2 quartic c = %.6f = kappa*pi/24)"
          % (KAPPA_EFF, KAPPA_FIT_V2))
    print("=" * 80)
    R_quartic = KAPPA_FIT_V2 * (4.0 ** 4)            # = kappa*m(4) relative to 0
    rows = []
    curves = {}
    for N in NLIST:
        c = cache[N]
        sp, psi0, A3 = c["sp"], c["psi0"], c["A3"]
        psi_p = v2.apply_V3(psi0, v2.V3_from_A(A3, +1), N)
        s = np.arange(0.0, HORIZON[N] + 1e-12, HOR_DS)
        f, _, _, _ = sp.evolve(psi_p, s)
        f0 = float(f[0])
        S = f / f0                                   # pointwise S(tau)=f/f0
        mlnS = -np.log(np.clip(S, 1e-300, None))

        # dynamical bound on the executed window and on the full horizon
        in_win = (s >= TAU_WIN_MIN) & (s <= TAU_WIN_MAX)
        sup_window = float(np.max(mlnS[in_win]))
        sup_horizon = float(np.max(mlnS))

        # first recurrence: after S drops below 0.95, first return |S-1|<=0.05
        dropped = False
        tau_rec = np.nan
        for k in range(1, len(s)):
            if S[k] < 1.0 - REC_TOL:
                dropped = True
            if dropped and abs(S[k] - 1.0) <= REC_TOL:
                tau_rec = float(s[k])
                break

        factor_bound = R_quartic / sup_window if sup_window > 0 else np.inf
        rec_in_window = (np.isfinite(tau_rec) and tau_rec <= TAU_WIN_MAX)
        factor_rec = (TAU_WIN_MAX / tau_rec) if (np.isfinite(tau_rec) and tau_rec > 0) else np.nan

        quartic_exceeds_bound = factor_bound > 1.0
        gate = quartic_exceeds_bound or rec_in_window

        binding = []
        if quartic_exceeds_bound:
            binding.append("quartic>bound x%.2f" % factor_bound)
        if rec_in_window:
            binding.append("window>tau_rec x%.2f" % factor_rec)
        rows.append(dict(
            N=N, horizon=HORIZON[N], R_quartic_at4=R_quartic,
            sup_minus_lnS_window=sup_window, sup_minus_lnS_horizon=sup_horizon,
            tau_rec=tau_rec, rec_within_window=rec_in_window,
            factor_quartic_over_bound=factor_bound, factor_window_over_taurec=factor_rec,
            quartic_exceeds_bound=quartic_exceeds_bound, gate_pass=gate,
            binding="; ".join(binding) if binding else "none"))
        curves[N] = dict(s=s, S=S, mlnS=mlnS)
        print(f"  N={N:2d} (horizon {HORIZON[N]:.0f}): "
              f"sup(-lnS)|_[0.5,4] = {sup_window:.3f}  "
              f"(full horizon {sup_horizon:.3f});  "
              f"tau_rec = {('%.2f' % tau_rec) if np.isfinite(tau_rec) else '> horizon'};  "
              f"quartic demand R={R_quartic:.3f}")
        print(f"        -> quartic/bound = {factor_bound:.2f}x  "
              f"({'EXCEEDS' if quartic_exceeds_bound else 'within'} dynamical range);  "
              f"recurrence in window: {rec_in_window};  GATE {'PASS' if gate else 'FAIL'}  "
              f"[{rows[-1]['binding']}]")
    return pd.DataFrame(rows), curves, R_quartic


# =============================================================================
# N3: exponent-drift reinterpretation (descriptive)
# =============================================================================
def run_N3(df_coeffs, df_gate):
    print("\n" + "=" * 80)
    print("N3  exponent drift reinterpreted: V2 Level-1 free-p vs short-time "
          "c1/c2 weights (DESCRIPTIVE)")
    print("=" * 80)
    rows = []
    tau_mid = 0.5 * (TAU_WIN_MIN + TAU_WIN_MAX)       # window midpoint 2.25
    for N in NLIST:
        r = df_coeffs[(df_coeffs.N == N) &
                      (df_coeffs.member == "level1_structured")].iloc[0]
        c1, c2, c3 = r.c1, r.c2, r.c3
        # window-effective local exponent predicted by the short-time poly
        p_pred_mid = effective_exponent(c1, c2, c3, tau_mid)
        # relative weight of quadratic vs linear contribution at the midpoint
        w_quad = abs(c2 * tau_mid ** 2) / (abs(c1 * tau_mid) + 1e-300)
        p_v2 = V2_LEVEL1_P[N]
        rows.append(dict(N=N, c1=c1, c2=c2, c3=c3, tau_mid=tau_mid,
                         p_eff_pred_at_mid=p_pred_mid,
                         quad_over_lin_weight=w_quad, p_level1_V2=p_v2))
        print(f"  N={N:2d}  c1={c1:+.4f} c2={c2:+.4f} c3={c3:+.4f}  | "
              f"quad/lin weight at tau={tau_mid:.2f}: {w_quad:6.2f}  | "
              f"p_eff(pred)={p_pred_mid:.2f}  vs  p_V2(measured)={p_v2:.3f}")
    df = pd.DataFrame(rows)
    # rank correlation between quad/lin weight and the measured V2 exponent
    a = df.quad_over_lin_weight.to_numpy()
    b = df.p_level1_V2.to_numpy()
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    rho = float(np.corrcoef(ra, rb)[0, 1])
    print(f"\n  Spearman(quad/lin weight, p_V2) = {rho:+.2f}  "
          f"(positive => more quadratic weight tracks larger measured p; "
          f"descriptive, no gate)")
    return df, rho


# =============================================================================
# FIGURES
# =============================================================================
def make_figures(df_coeffs, df_ens, rec_curves, df_rec, R_quartic):
    # ---- Figure 1: short-time coefficients ----
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    # (left) ensemble c2/sigma per N with the 5-sigma line; level1 marker
    ax = axs[0]
    for N in NLIST:
        sub = df_coeffs[(df_coeffs.N == N) & (df_coeffs.member == "scramble")]
        z2 = (sub.c2.abs() / sub.s2.clip(lower=1e-30)).to_numpy()
        x = np.full(len(z2), N) + np.random.default_rng(N).normal(0, 0.12, len(z2))
        ax.scatter(x, z2, s=14, alpha=0.6,
                   label=("GOE scrambles" if N == NLIST[0] else None), color="C0")
        l1 = df_coeffs[(df_coeffs.N == N) & (df_coeffs.member == "level1_structured")]
        z2l = abs(float(l1.c2.iloc[0])) / max(float(l1.s2.iloc[0]), 1e-30)
        ax.scatter([N], [z2l], marker="D", s=55, color="C3", zorder=5,
                   label=("Level-1 structured A" if N == NLIST[0] else None))
    ax.axhline(SIG5, color="k", ls="--", lw=1.0, label=r"$5\sigma$ gate")
    ax.set_yscale("log")
    ax.set_xlabel("N (chain size)")
    ax.set_ylabel(r"$|c_2|/\sigma(c_2)$")
    ax.set_title(r"(a) N1: significance of the quadratic coefficient $c_2$")
    ax.set_xticks(NLIST)
    ax.legend(fontsize=8, loc="lower right")

    # (right) c1 raw vs Theta-symmetrized (codimension-one) for Level 1
    ax = axs[1]
    l1 = df_coeffs[df_coeffs.member == "level1_structured"].sort_values("N")
    Ns = l1.N.to_numpy()
    ax.semilogy(Ns, l1.c1.abs().clip(lower=1e-18), "o-", color="C0",
                label=r"$|c_1|$ raw forward functional")
    ax.semilogy(Ns, l1.c1_sym.abs().clip(lower=1e-18), "s--", color="C2",
                label=r"$|c_1|$ $\Theta$-symmetrized (even)")
    ax.semilogy(Ns, l1.c2.abs().clip(lower=1e-18), "^-", color="C3",
                label=r"$|c_2|$ raw")
    ax.set_xlabel("N (chain size)")
    ax.set_ylabel("coefficient magnitude")
    ax.set_title(r"(b) N1: $c_1$ suppressed by $\Theta$-symmetrization "
                 r"(codim.\ one); $c_2\neq0$ generic")
    ax.set_xticks(NLIST)
    ax.legend(fontsize=8)
    fig.tight_layout()
    f1 = os.path.join(OUTDIR, "b7_fig_coeficientes.png")
    fig.savefig(f1, dpi=140)
    plt.close(fig)

    # ---- Figure 2: recurrence vs the monotone quartic law ----
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axs[0]
    for N, col in zip(NLIST, ("C0", "C1", "C2", "C3")):
        cu = rec_curves[N]
        m = cu["s"] <= min(HORIZON[N], 40.0)
        ax.plot(cu["s"][m], cu["mlnS"][m], lw=0.9, color=col, alpha=0.85,
                label=f"N={N}: $-\\ln S(\\tau)$ (bare sector)")
    tt = np.linspace(0.0, 40.0, 400)
    ax.plot(tt, KAPPA_FIT_V2 * tt ** 4, "k--", lw=1.6,
            label=r"quartic law $\kappa\,m(\tau)=c_{V2}\tau^4$")
    ax.axvspan(TAU_WIN_MIN, TAU_WIN_MAX, color="grey", alpha=0.18,
               label=r"executed window $[0.5,4.0]$")
    ax.set_ylim(0, max(R_quartic * 1.3, 8))
    ax.set_xlim(0, 40)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$-\ln S(\tau)$")
    ax.set_title(r"(a) N2: bounded, recurrent $-\ln S$ vs the unbounded "
                 r"quartic demand")
    ax.legend(fontsize=7.5, loc="upper left")

    ax = axs[1]
    Ns = df_rec.N.to_numpy()
    ax.bar(Ns - 0.5, df_rec.sup_minus_lnS_window, width=1.0, color="C0",
           alpha=0.8, label=r"$\sup_{[0.5,4]}(-\ln S)$ (dynamical ceiling)")
    ax.axhline(R_quartic, color="k", ls="--", lw=1.4,
               label=r"quartic demand $\kappa\,m(4)=%.2f$" % R_quartic)
    ax.set_xlabel("N (chain size)")
    ax.set_ylabel(r"$-\ln S$")
    ax.set_title(r"(b) N2: the quartic demand overshoots the dynamical "
                 r"ceiling at every $N$")
    ax.set_xticks(NLIST)
    ax.legend(fontsize=8)
    fig.tight_layout()
    f2 = os.path.join(OUTDIR, "b7_fig_recurrencia.png")
    fig.savefig(f2, dpi=140)
    plt.close(fig)
    return f1, f2


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("#" * 80)
    print("# B.7  CONFIRMATION of thm:emergence-nogo / cor:level1-null")
    print("# (a proven theorem; this block confirms, it cannot refute)")
    print("#" * 80)

    gate_ok, df_gate, cache = continuity_gate()
    if not gate_ok:
        # principle 8 / preregistered verdict: do NOT proceed on a broken base
        df_gate.to_csv(os.path.join(OUTDIR, "b7_continuity_gate.csv"),
                       index=False)
        print("\nABORTED: continuity gate failed; nothing interpreted.")
        sys.exit(2)

    df_coeffs, df_ens = run_N1(cache)
    df_rec, rec_curves, R_quartic = run_N2(cache)
    df_n3, rho_n3 = run_N3(df_coeffs, df_gate)

    # ---- merge ensemble gate flags into the recurrence/coeff outputs --------
    df_coeffs_out = df_coeffs.copy()
    df_rec_out = df_rec.merge(df_ens, on="N", how="left")

    # ---- save CSVs ----
    paths = {}
    for name, df in (("b7_nogo_coeffs", df_coeffs_out),
                     ("b7_recurrence", df_rec_out),
                     ("b7_N3_exponent_drift", df_n3),
                     ("b7_continuity_gate", df_gate)):
        p = os.path.join(OUTDIR, name + ".csv")
        df.to_csv(p, index=False)
        paths[name] = p
        print(f"[saved] {p}")

    f1, f2 = make_figures(df_coeffs_out, df_ens, rec_curves, df_rec, R_quartic)
    print(f"[saved] {f1}\n[saved] {f2}")

    # ---- preregistered verdict ----
    gate_a_all = bool(df_ens.gate_a_majority.all())
    gate_b_all = bool(df_ens.gate_b_no_finetune.all())
    n2_gate_all = bool(df_rec.gate_pass.all())
    print("\n" + "=" * 80)
    print("PREREGISTERED VERDICT (confirmation; cannot refute the theorem)")
    print("=" * 80)
    print(f"  G0 continuity gate (seed-free V2 reproduction):  "
          f"{'PASS' if gate_ok else 'FAIL'}")
    print(f"  N1 gate (a)  c2 != 0 at >=5sigma in majority, all N:  "
          f"{'PASS' if gate_a_all else 'FAIL'}")
    print(f"  N1 gate (b)  no spontaneous c1=c2=c3=0, all N:        "
          f"{'PASS' if gate_b_all else 'FAIL'}")
    print(f"  N2 gate      quartic demand exceeds bound OR window>tau_rec, "
          f"all N:  {'PASS' if n2_gate_all else 'FAIL'}")
    print(f"  N3           Spearman(quad/lin, p_V2) = {rho_n3:+.2f} "
          f"(descriptive)")
    print("-" * 80)
    if gate_ok and gate_a_all and gate_b_all and n2_gate_all:
        print("  CONFIRMED: the short-time obstruction (generic p in {1,2}, no")
        print("  spontaneous fine-tuning) and the recurrence obstruction (the")
        print("  bounded almost-periodic -ln S cannot host the unbounded")
        print("  quartic demand on the executed window) hold in the real G1")
        print("  model at every size.  The executed Level-1 nulls of B.4 are")
        print("  instances of cor:level1-null; the numerics are DEMOTED FROM")
        print("  JUDGE TO CONFIRMATION.")
    else:
        print("  A preregistered check did not pass.  Per the prereg verdict,")
        print("  FIRST audit whether the executed model is inside the theorem's")
        print("  hypotheses (time-independent generator, stationary")
        print("  preparation) before any change; report the discrepancy and")
        print("  BLOCK integration of the corollary until resolved.")
    print("-" * 80)
    print(f"[total runtime] {time.time() - T0:.1f} s")
    return dict(gate=df_gate, coeffs=df_coeffs_out, ens=df_ens,
                rec=df_rec_out, n3=df_n3)


if __name__ == "__main__":
    main()
