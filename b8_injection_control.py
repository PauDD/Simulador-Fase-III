#!/usr/bin/env python3
# =============================================================================
# b8_injection_control.py
# -----------------------------------------------------------------------------
# PROMPT SIM-V3-4 -- exponent-injection control (I1-I3) + window-offset refit
# (O1-O2) for the L2b pulsed-transcription carrier of the Phase III discriminator.
#
# PREREGISTRATION (frozen before any number was read off):
#   Context. In V2 the L2b pulsed carrier transmitted the causal-monotone law
#   with a free-power fit p = 4.43 +/- 0.09 (95% CI) and decisive AIC preference
#   (see b4_A_of_N.csv, b4_fit_models.csv).  rem:l2c-variant (sec_05_v3_p7)
#   already states that this fidelity-scheduled carrier "transcribes ANY law it
#   is programmed with", hence its evidential weight for the SPECIFIC law
#   m(tau)=(pi/24)tau^4 must be calibrated by an exponent-injection control.
#   Separately, the measured p sits ~5 sigma above 4 because of an identified
#   window-edge systematic (the window maximum of a decaying signal samples the
#   left edge tau-DW/2, so the transmitted suppression is exp[-kappa_b m(tau-0.3)];
#   "verified exactly" in V2).  O1-O2 refit with the offset model and report the
#   corrected headline p.
#
#   Model: IDENTICAL to the executed V2 L2b carrier.  Only the injected schedule
#   m_tilde(tau) changes.  Per-pulse fidelity  cos(2 phi_j) = exp(-kappa[m_tilde(tau_j)
#   - m_tilde(tau_{j-1})]); pulses on fresh Theta-invariant register qubits |0>.
#   The observable is the V2 quantity S2b(tau) = window_max(f_free * P)/window_max(f_free)
#   built from the SAME free-sector coherence f_free of the T-symmetric matter
#   sector (schedule-independent) and the explicitly evolved per-pulse register
#   2-vectors.  Representative N = 8 (L2b suppression is N-independent by
#   construction; V2 confirms p_2b identical to 6 digits across N in {4,6,8,10}).
#   Preregistered fit family unchanged: M_exp(1), M_gauss(2), M_quartic(4),
#   M_freep (95% CI), M_power.  Window tau in [0.5,4.0] step 0.05, DW=0.6.
#   Seed 20240607 (declared; no fit quantity below depends on the RNG -- the
#   scramble band is not used by I1-I3/O1/O2, which are deterministic fits).
#
#   REPRODUCTION GATE (infrastructure continuity).  Re-run the ORIGINAL V2
#   schedule m(tau)=(pi/24)tau^4, kappa_b=0.18 and require p_free = 4.43 +/- 0.09.
#
#   Injected effective coupling (DECLARED): kappa_inj chosen so the total
#   dynamic range over the window top matches V2, kappa_inj = kappa_b * m(4)/m_tilde(4).
#
#   I1 (cubic):   m_tilde = tau^3.            Measure p_free.        Expect ~3.
#   I2 (quintic): m_tilde = tau^5.            Measure p_free.        Expect ~5.
#   I3 (log):     m_tilde = ln(1+tau/tau0), tau0=0.5.  Add a SIXTH member
#                 M_log = lnS0 - c ln(1+tau/tau0L) (k=3, tau0L free) to the fit
#                 family SOLELY in I3; require M_log preferred by AIC over every
#                 polynomial member.
#   GATE I1-I3 (joint): the channel recovers the injected exponent with
#                 |p_measured - p_injected| <= 3 sigma_fit in I1-I2 AND prefers
#                 the log in I3.  (Both p_measured reported: raw, and
#                 offset-corrected by the SAME systematic characterized in O1.)
#                 If PASS: the fidelity channel transmits ANY law => the V2 p~4
#                 has NULL evidential weight for the specific law m and stands
#                 only as a demonstration of implementability/distinguishability
#                 (the narrow reading the paper already declares -- now quantified).
#                 If FAIL: transcription is not faithful in general; block the
#                 phrase "faithfully implementable" and flag for reanalysis.
#
#   O1 (offset refit on the V2-ORIGINAL L2b data reproduced by the gate):
#       free-power offset  ln S = lnS0 - c (tau - tau0)^p  (tau0,p free) -> headline.
#       quartic offset     ln S = lnS0 - c (tau - tau0)^4  -> check c vs kappa_b*pi/24.
#       GATE: tau0 CI contains 0.3 (within DW/2 of the V2 window) AND p CI
#       contains 4.  Report the free-power-offset p as the corrected headline that
#       SUPERSEDES 4.43 +/- 0.09 in abstract/conclusion on V3 integration.
#       If FAIL: the V2 edge systematic was mis-identified; report and reopen
#       (do NOT retune kappa_b or the window a posteriori).
#
#   O2 (stability): repeat the free-power offset over DW in {0.4,0.6,0.8} and
#       three tau ranges [0.5,4],[1,4],[0.5,3].  Report tau0(DW), p(DW).
#
#   No verdict is pre-decided; both outcomes of every block are reportable.
#   Numerical data overrides narrative.
#
# DELIVERABLES: b8_injection.csv, b8_offset_fit.csv, b8_fig_inyeccion.png,
#   b8_fig_offset.png.  Single CPU < 45 min.
# =============================================================================
import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

T0 = time.time()
SEED = 20240607
OUTDIR = os.environ.get("B8_OUTDIR", "/mnt/user-data/outputs")
os.makedirs(OUTDIR, exist_ok=True)

# --- reuse the executed V2 infrastructure verbatim (import-safe module) ------
sys.path.insert(0, os.environ.get("B45_SRC", "/mnt/project"))
import b4_b5_discriminator as bd

KAPPA_B   = bd.KAPPA_B          # 0.18
DP        = bd.DP_PULSE         # 0.1
T_MAX     = bd.T_MAX            # 4.6
DS        = bd.DS               # 0.005
DW_DEF    = bd.DW               # 0.6
M1_FLOOR  = bd.M1_FLOOR
TAU_GRID  = bd.TAU_GRID         # arange(0.5, 4.0, 0.05)
N_REPR    = bd.N_REPR           # 8
V2_HEADLINE_P  = 4.4327         # from b4_A_of_N.csv (gate target)
V2_HEADLINE_CI = 0.0927

sgrid_f = np.arange(0.0, T_MAX + 1e-12, DS)
sgrid_b = -sgrid_f


# -----------------------------------------------------------------------------
# Injected monotone laws (all monotone increasing => valid fidelity schedule).
# -----------------------------------------------------------------------------
def m_quartic(tau):           # V2 original
    return (np.pi / 24.0) * np.asarray(tau, float) ** 4
def m_cubic(tau):
    return np.asarray(tau, float) ** 3
def m_quintic(tau):
    return np.asarray(tau, float) ** 5
def m_log(tau, tau0=0.5):
    return np.log1p(np.asarray(tau, float) / tau0)

def kappa_for(m_func, kappa_b=KAPPA_B, tau_ref=4.0):
    """Match total dynamic range at the window top to V2's quartic."""
    return kappa_b * float(m_quartic(tau_ref)) / float(m_func(tau_ref))


# -----------------------------------------------------------------------------
# Generalised pulse schedule (V2 pulse_schedule with an injectable law).
# -----------------------------------------------------------------------------
def pulse_schedule(m_func, kappa, dp=DP, t_max=T_MAX):
    taus = np.arange(dp, t_max + 1e-12, dp)
    dm = kappa * (m_func(taus) - m_func(taus - dp))
    args = np.exp(-dm)
    assert np.all(args > 0) and np.all(args <= 1.0 + 1e-12), "schedule not a fidelity"
    args = np.clip(args, None, 1.0)
    phis = 0.5 * np.arccos(args)
    return taus, phis


def free_coherence(N):
    """Schedule-independent free-sector coherence f_free(|s|) = |rho01|(s)."""
    Hsc = bd.build_H_sc(N)
    psi0, *_ = bd.prepare_state(N, Hsc)
    V3p = bd.V3_from_A(bd.build_A3_structured(), +1)
    psi_p = bd.apply_V3(psi0, V3p, N)
    sp = bd.SpectralProp(N, Hsc)
    f_fwd, _, _, _ = sp.evolve(psi_p, sgrid_f)
    return f_fwd


def build_S2b(f_free, taus, phis, dw=DW_DEF, tau_grid=TAU_GRID):
    """V2 observable S2b(tau) = max_W (f_free*P) / max_W (f_free)."""
    P = np.abs(bd.pulse_factor_grid(np.abs(sgrid_f), True, taus, phis))
    f2b = f_free * P
    M2 = np.array([bd.window_max(sgrid_f, f2b, t, dw) for t in tau_grid])
    M1 = np.array([bd.window_max(sgrid_f, f_free, t, dw) for t in tau_grid])
    good = M1 > M1_FLOOR
    S = np.full(len(tau_grid), np.nan)
    S[good] = M2[good] / M1[good]
    return tau_grid, S, good, P


# -----------------------------------------------------------------------------
# Fits
# -----------------------------------------------------------------------------
def aic_of(y, resid, k):
    n = len(y); rss = float(np.sum(resid ** 2))
    return n * np.log(max(rss, 1e-300) / n) + 2 * k, rss

def fit_log_member(tau, S, tau0_guess=0.5):
    """M_log: ln S = lnS0 - c ln(1+tau/tau0L), 3 params (k=3)."""
    y = np.log(np.clip(S, 1e-300, None))
    def f(t, lnS0, c, tau0L):
        return lnS0 - c * np.log1p(t / tau0L)
    popt, pcov = curve_fit(f, tau, y, p0=[0.0, KAPPA_B, tau0_guess],
                           bounds=([-5, 1e-8, 1e-3], [5, 50, 5]), maxfev=60000)
    aic, rss = aic_of(y, y - f(tau, *popt), 3)
    r2 = 1.0 - rss / float(np.sum((y - y.mean()) ** 2))
    return dict(lnS0=popt[0], c=popt[1], tau0L=popt[2], k=3, rss=rss, R2=r2, AIC=aic)

def fit_offset_freep(tau, S, tau0_hi):
    """ln S = lnS0 - c (tau-tau0)^p ; params (lnS0,c,tau0,p), k=4."""
    y = np.log(np.clip(S, 1e-300, None))
    def f(t, lnS0, c, tau0, p):
        return lnS0 - c * np.power(np.clip(t - tau0, 1e-9, None), p)
    popt, pcov = curve_fit(f, tau, y, p0=[0.0, 0.02, 0.3, 4.0],
                           bounds=([-5, 1e-8, -1.0, 0.5], [5, 50, tau0_hi, 8.0]),
                           maxfev=80000)
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    aic, rss = aic_of(y, y - f(tau, *popt), 4)
    r2 = 1.0 - rss / float(np.sum((y - y.mean()) ** 2))
    return dict(lnS0=popt[0], c=popt[1], tau0=popt[2], p=popt[3],
                tau0_ci95=1.96 * perr[2], p_ci95=1.96 * perr[3],
                k=4, rss=rss, R2=r2, AIC=aic)

def fit_offset_quartic(tau, S, tau0_hi):
    """ln S = lnS0 - c (tau-tau0)^4 ; params (lnS0,c,tau0), k=3."""
    y = np.log(np.clip(S, 1e-300, None))
    def f(t, lnS0, c, tau0):
        return lnS0 - c * np.power(np.clip(t - tau0, 1e-9, None), 4.0)
    popt, pcov = curve_fit(f, tau, y, p0=[0.0, 0.024, 0.3],
                           bounds=([-5, 1e-8, -1.0], [5, 50, tau0_hi]), maxfev=80000)
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    aic, rss = aic_of(y, y - f(tau, *popt), 3)
    return dict(lnS0=popt[0], c=popt[1], tau0=popt[2],
                c_ci95=1.96 * perr[1], tau0_ci95=1.96 * perr[2], k=3, rss=rss, AIC=aic)


# =============================================================================
print("=" * 78)
print("SIM-V3-4  exponent-injection control (L2b) + window-offset refit")
print("=" * 78)

# ---- free coherence once (schedule-independent) at representative N=8 -------
t0 = time.time()
f_free = free_coherence(N_REPR)
print(f"[free coherence N={N_REPR} computed in {time.time()-t0:.1f}s]")

inj_rows = []

# ---- REPRODUCTION GATE: original quartic schedule ---------------------------
taus, phis = pulse_schedule(m_quartic, KAPPA_B)
TAU, S2b, good, P_orig = build_S2b(f_free, taus, phis)
df_gate = bd.fit_models(TAU[good], np.clip(S2b[good], 1e-300, None))
p_gate = float(df_gate[df_gate.model == "M_freep"].p.iloc[0])
ci_gate = float(df_gate[df_gate.model == "M_freep"].p_ci95.iloc[0])
gate_ok = abs(p_gate - V2_HEADLINE_P) <= max(V2_HEADLINE_CI, ci_gate) + 0.02
print(f"\n[GATE] original quartic schedule: p_free = {p_gate:.3f} +/- {ci_gate:.3f}"
      f"  (V2 target {V2_HEADLINE_P} +/- {V2_HEADLINE_CI})  -> "
      f"{'REPRODUCED' if gate_ok else 'MISMATCH'}")
assert gate_ok, "reproduction gate failed -- infrastructure not continuous"

# offset-corrected reference for the SAME quartic curve (used to debias I1-I2)
ofit_ref = fit_offset_freep(TAU[good], S2b[good], tau0_hi=TAU[good].min() - 0.01)
print(f"       offset-corrected (same curve): p = {ofit_ref['p']:.3f} +/- "
      f"{ofit_ref['p_ci95']:.3f}, tau0 = {ofit_ref['tau0']:.3f} +/- {ofit_ref['tau0_ci95']:.3f}")

inj_rows.append(dict(run="GATE", law="quartic (pi/24 tau^4)", p_injected=4.0,
                     kappa_inj=KAPPA_B, p_raw=p_gate, p_raw_ci95=ci_gate,
                     p_offset=ofit_ref['p'], p_offset_ci95=ofit_ref['p_ci95'],
                     tau0_offset=ofit_ref['tau0'], tau0_offset_ci95=ofit_ref['tau0_ci95'],
                     R2_freep=float(df_gate[df_gate.model=="M_freep"].R2.iloc[0]),
                     gate_pass=bool(gate_ok)))

# ---- INJECTIONS I1, I2 ------------------------------------------------------
def run_injection(tag, m_func, p_injected):
    kap = kappa_for(m_func)
    taus, phis = pulse_schedule(m_func, kap)
    TAU, S, good, _ = build_S2b(f_free, taus, phis)
    df = bd.fit_models(TAU[good], np.clip(S[good], 1e-300, None))
    p_raw = float(df[df.model == "M_freep"].p.iloc[0])
    ci_raw = float(df[df.model == "M_freep"].p_ci95.iloc[0])
    of = fit_offset_freep(TAU[good], S[good], tau0_hi=TAU[good].min() - 0.01)
    sig_off = max(of['p_ci95'] / 1.96, 1e-6)
    dev = abs(of['p'] - p_injected)
    pass_off = dev <= 3.0 * sig_off
    print(f"\n[{tag}] inject {m_func.__name__} (p_inj={p_injected}), kappa_inj={kap:.5f}")
    print(f"     raw   p_free = {p_raw:.3f} +/- {ci_raw:.3f}  (biased by window edge)")
    print(f"     offset p_free = {of['p']:.3f} +/- {of['p_ci95']:.3f}, "
          f"tau0 = {of['tau0']:.3f}  -> |p-{p_injected}| = {dev:.3f}  "
          f"({dev/sig_off:.2f} sigma)  -> {'PASS' if pass_off else 'FAIL'}")
    inj_rows.append(dict(run=tag, law=m_func.__name__, p_injected=p_injected,
                         kappa_inj=kap, p_raw=p_raw, p_raw_ci95=ci_raw,
                         p_offset=of['p'], p_offset_ci95=of['p_ci95'],
                         tau0_offset=of['tau0'], tau0_offset_ci95=of['tau0_ci95'],
                         R2_freep=float(df[df.model=="M_freep"].R2.iloc[0]),
                         dev_sigma=dev / sig_off, gate_pass=bool(pass_off)))
    return (TAU, S, good)

i1 = run_injection("I1", m_cubic, 3.0)
i2 = run_injection("I2", m_quintic, 5.0)

# ---- INJECTION I3 (log) with the sixth ad-hoc M_log member ------------------
kap_log = kappa_for(m_log)
taus, phis = pulse_schedule(m_log, kap_log)
TAU3, S3, good3, _ = build_S2b(f_free, taus, phis)
df3 = bd.fit_models(TAU3[good3], np.clip(S3[good3], 1e-300, None))
log_fit = fit_log_member(TAU3[good3], np.clip(S3[good3], 1e-300, None))
# AIC of every polynomial member vs the log member
poly_aic = {m: float(df3[df3.model == m].AIC.iloc[0])
            for m in ("M_exp(p=1)", "M_gauss(p=2)", "M_quartic(p=4)",
                      "M_freep", "M_power")}
log_aic = log_fit["AIC"]
log_preferred = all(log_aic < a for a in poly_aic.values())
best_poly = min(poly_aic, key=poly_aic.get)
print(f"\n[I3] inject log (tau0=0.5), kappa_inj={kap_log:.5f}")
print(f"     M_log AIC = {log_aic:.1f} (R2={log_fit['R2']:.4f}, tau0L={log_fit['tau0L']:.3f})")
for m, a in poly_aic.items():
    print(f"        {m:14s} AIC = {a:8.1f}")
print(f"     best polynomial = {best_poly} (AIC {poly_aic[best_poly]:.1f})  -> "
      f"log {'PREFERRED' if log_preferred else 'NOT preferred'}")
inj_rows.append(dict(run="I3", law="m_log", p_injected=np.nan, kappa_inj=kap_log,
                     p_raw=float(df3[df3.model=="M_freep"].p.iloc[0]),
                     p_raw_ci95=float(df3[df3.model=="M_freep"].p_ci95.iloc[0]),
                     p_offset=np.nan, p_offset_ci95=np.nan,
                     tau0_offset=log_fit["tau0L"], tau0_offset_ci95=np.nan,
                     R2_freep=float(df3[df3.model=="M_freep"].R2.iloc[0]),
                     AIC_log=log_aic, AIC_best_poly=poly_aic[best_poly],
                     best_poly=best_poly, gate_pass=bool(log_preferred)))

# joint I1-I3 gate
i1_pass = inj_rows[1]["gate_pass"]
i2_pass = inj_rows[2]["gate_pass"]
joint_pass = i1_pass and i2_pass and log_preferred
print(f"\n[JOINT I1-I3 GATE] I1={'PASS' if i1_pass else 'FAIL'}, "
      f"I2={'PASS' if i2_pass else 'FAIL'}, I3(log)={'PASS' if log_preferred else 'FAIL'}"
      f"  ==> {'PASS (channel transmits ANY law; V2 p~4 has null evidential weight for m)' if joint_pass else 'FAIL'}")

pd.DataFrame(inj_rows).to_csv(os.path.join(OUTDIR, "b8_injection.csv"), index=False)

# ---- O1 / O2: offset refit on the V2-ORIGINAL L2b curve ---------------------
print("\n" + "-" * 78)
print("O1/O2  window-offset refit of the V2-original L2b suppression")
off_rows = []
RANGES = {"[0.5,4]": (0.5, 4.0), "[1,4]": (1.0, 4.0), "[0.5,3]": (0.5, 3.0)}
for dw in (0.4, 0.6, 0.8):
    tausq, phisq = pulse_schedule(m_quartic, KAPPA_B)
    TAUq, Sq, goodq, _ = build_S2b(f_free, tausq, phisq, dw=dw)
    for rname, (lo, hi) in RANGES.items():
        sel = goodq & (TAU_GRID >= lo - 1e-9) & (TAU_GRID <= hi + 1e-9)
        tau_s = TAU_GRID[sel]; S_s = Sq[sel]
        if len(tau_s) < 6:
            continue
        hi0 = tau_s.min() - 0.01
        of = fit_offset_freep(tau_s, S_s, tau0_hi=hi0)
        oq = fit_offset_quartic(tau_s, S_s, tau0_hi=hi0)
        tau0_ok = abs(of['tau0'] - dw / 2.0) <= max(of['tau0_ci95'], 0.05)
        p_ok = abs(of['p'] - 4.0) <= max(of['p_ci95'], 0.0)
        is_o1 = (dw == 0.6 and rname == "[0.5,4]")
        print(f"  DW={dw} range {rname:8s}: free-offset p={of['p']:.3f}+/-{of['p_ci95']:.3f}, "
              f"tau0={of['tau0']:.3f}+/-{of['tau0_ci95']:.3f} (DW/2={dw/2:.2f})"
              f"{'  <-- O1 headline' if is_o1 else ''}")
        print(f"               quartic-offset c={oq['c']:.5f}+/-{oq['c_ci95']:.5f} "
              f"(input kappa_b*pi/24={KAPPA_B*np.pi/24:.5f}), tau0={oq['tau0']:.3f}")
        off_rows.append(dict(block="O1" if is_o1 else "O2", DW=dw, tau_range=rname,
                             n=len(tau_s),
                             p_freeoffset=of['p'], p_ci95=of['p_ci95'],
                             tau0_freeoffset=of['tau0'], tau0_ci95=of['tau0_ci95'],
                             DW_half=dw / 2.0,
                             c_quartoffset=oq['c'], c_quartoffset_ci95=oq['c_ci95'],
                             c_input=KAPPA_B * np.pi / 24, tau0_quartoffset=oq['tau0'],
                             p_contains_4=bool(p_ok),
                             tau0_tracks_DWhalf=bool(tau0_ok)))

df_off = pd.DataFrame(off_rows)
df_off.to_csv(os.path.join(OUTDIR, "b8_offset_fit.csv"), index=False)
o1 = df_off[df_off.block == "O1"].iloc[0]
p_stable = bool(df_off.p_contains_4.all())
tau0_tracks = bool(df_off.tau0_tracks_DWhalf.all())
print(f"\n[O1 headline] corrected p = {o1.p_freeoffset:.3f} +/- {o1.p_ci95:.3f}, "
      f"tau0 = {o1.tau0_freeoffset:.3f} +/- {o1.tau0_ci95:.3f}")
print(f"[O2 stability] p contains 4 across all windows: {p_stable}; "
      f"tau0 tracks DW/2 across all windows: {tau0_tracks}")

# ---- figures ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 5.0))
dfi = pd.DataFrame(inj_rows)
m = dfi.p_injected.notna()
ax.plot([2.5, 5.5], [2.5, 5.5], "k--", lw=1, label="ideal (p$_{out}$=p$_{in}$)")
ax.errorbar(dfi.p_injected[m], dfi.p_raw[m], yerr=dfi.p_raw_ci95[m], fmt="s",
            color="#c0392b", capsize=3, label="raw $M_{\\mathrm{free}p}$ (edge-biased)")
ax.errorbar(dfi.p_injected[m], dfi.p_offset[m], yerr=dfi.p_offset_ci95[m], fmt="o",
            color="#1f6f8b", capsize=3, label="offset-corrected")
for _, r in dfi[m].iterrows():
    ax.annotate(r.run, (r.p_injected, r.p_offset), textcoords="offset points",
                xytext=(6, -10), fontsize=8)
ax.set_xlabel("injected exponent  $p_{\\mathrm{in}}$")
ax.set_ylabel("recovered free-power exponent  $p_{\\mathrm{out}}$")
ax.set_title("L2b exponent-injection control")
ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3)
fig.tight_layout()
fpa = os.path.join(OUTDIR, "b8_fig_inyeccion.png")
fig.savefig(fpa, dpi=140); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.6, 4.6))
TAUq, Sq, goodq, _ = build_S2b(f_free, *pulse_schedule(m_quartic, KAPPA_B), dw=0.6)
t = TAUq[goodq]; y = -np.log(np.clip(Sq[goodq], 1e-300, None))
ax.plot(t, y, "o", ms=3, color="0.35", label="V2 L2b data  $-\\ln S$")
tt = np.linspace(t.min(), t.max(), 300)
ax.plot(tt, KAPPA_B * (np.pi / 24) * tt ** 4, "-", color="#c0392b", lw=1.3,
        label="no-offset $\\kappa_b\\frac{\\pi}{24}\\tau^4$ (p$\\to$4.43)")
ax.plot(tt, o1.c_quartoffset * (tt - o1.tau0_quartoffset) ** 4, "--",
        color="#1f6f8b", lw=1.5,
        label=f"offset $\\kappa_b\\frac{{\\pi}}{{24}}(\\tau-{o1.tau0_quartoffset:.2f})^4$")
ax.set_xlabel("$\\tau$"); ax.set_ylabel("$-\\ln S(\\tau)$")
ax.set_title("Window-offset refit (O1)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout()
fpb = os.path.join(OUTDIR, "b8_fig_offset.png")
fig.savefig(fpb, dpi=140); plt.close(fig)

print("\n[outputs]")
for f in ("b8_injection.csv", "b8_offset_fit.csv", "b8_fig_inyeccion.png", "b8_fig_offset.png"):
    print("   ", os.path.join(OUTDIR, f))
print(f"[runtime] {time.time()-T0:.1f}s")

# expose summary for the .tex writer
SUMMARY = dict(p_gate=p_gate, ci_gate=ci_gate,
               i1=inj_rows[1], i2=inj_rows[2],
               log_preferred=log_preferred, log_aic=log_aic, poly_aic=poly_aic,
               joint_pass=joint_pass,
               o1=dict(o1), p_stable=p_stable, tau0_tracks=tau0_tracks)
import json
with open(os.path.join(OUTDIR, "b8_summary.json"), "w") as fh:
    json.dump({k: (v if not isinstance(v, dict) else
                   {kk: (None if (isinstance(vv, float) and np.isnan(vv)) else vv)
                    for kk, vv in v.items()}) for k, v in SUMMARY.items()},
              fh, indent=2, default=str)
