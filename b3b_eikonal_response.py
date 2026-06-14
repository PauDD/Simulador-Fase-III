#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
SIM-V3-1  --  Stationary-phase weight validation (extension of the B.3 layer).
              Numerical realization of the analytic sketches of
                  thm:eikonal-localization   (sec_05_v3_p3_eikonal.tex)
                  prop:stationary-weight     (sec_06_v3_p1_response.tex)
              for the first-order worldline-detector amplitude
                  A(Omega) = -i lambda \int dtau chi(tau) e^{i Omega tau} psi_B(tau),
                  psi_B(tau) = \int dk~ a(k) e^{i S(k)/eps} e^{i eta(k, gamma_B(tau))}.
              Signature (-,+,+,+); positive-frequency convention (k^0>0).

  Implemented in a 1+1D right-mover reduction, where the eikonal localization and
  the stationary-phase / Parseval bookkeeping are transparent and exact. This is a
  DECLARED modeling choice (principle 2), analogous to the tau=t identification
  declared in the B.4/B.5 layer: it is the regime in which the theorem's content
  is sharpest, not a hidden assumption.

--------------------------------------------------------------------------------
BINDING PRINCIPLES (1-4, 8): data govern; zero fabrication; preregistration;
zero silent exclusions; if a simulation contradicts a conceptual claim, the
simulation wins and the claim is flagged -- the theory text is NOT edited from
the sim.
--------------------------------------------------------------------------------
PREREGISTERED GATES (verbatim intent):
  V1  off-incidence |A|^2 decays faster than any power of eps (p_eff grows).
  V2  dominant tau-contribution width ~ eps^q, q=0.5 (95% CI contains 0.5).
  V3  (THE GATE) broadband weight P([Z]) ~ |a(k_Z)|^2 (J_AB)^beta, beta=1;
      |beta-1| <= 0.05.
  V4  two separated events: cross/diagonal ratio decays as in V1.
  control: numeric quadrature vs mpmath, rtol 1e-8 at control points.
VERDICT clause: if beta != 1 out of tolerance, the Hessian/measure bookkeeping
of prop:stationary-weight is reported mis-normalized; proposition flagged
pending analytic correction; theory unchanged. Both outcomes reportable.

Seed 20240604. eps mesh {0.1,0.05,0.02,0.01,0.005,0.002}. Budget < 30 min CPU.
================================================================================
"""
import os, csv, time
import numpy as np

SEED        = 20240604
EPS_MESH    = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002]
V_LIST      = [0.0, 0.2, 0.4, 0.6, 0.8]
W0          = 1.0
S_SIGMA     = 0.18
CW, CE      = 3.5, 1.0          # chi plateau half-width, edge roll-off
ACCEL       = 0.15             # proper acceleration for V2 (a*tau<=0.3 in window)
BETA_GATE   = 0.05
Q_TARGET    = 0.5
MPMATH_RTOL = 1e-8
NEARZERO    = 1e-9
N_CONTROL   = 6

OUTDIR = "/mnt/user-data/outputs"
os.makedirs(OUTDIR, exist_ok=True)
rng = np.random.default_rng(SEED)
EXCLUSIONS = []


# ---------------------------------------------------------------- primitives
def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def chi(tau):
    tau = np.asarray(tau, float); out = np.zeros_like(tau)
    a0, a1 = CW - CE, CW + CE; ax = np.abs(tau)
    out[ax <= a0] = 1.0
    e = (ax > a0) & (ax < a1); u = (ax[e] - a0) / (a1 - a0)
    f0 = np.exp(-1.0 / np.clip(1 - u, 1e-300, None))
    f1 = np.exp(-1.0 / np.clip(u, 1e-300, None))
    out[e] = f0 / (f0 + f1)
    return out


def profile(w, w0=W0, sig=S_SIGMA):
    return np.exp(-0.5 * ((w - w0) / sig) ** 2)


def eikonal_S(w, kind, w0=W0):
    return (w - w0) if kind == "linear" else 0.5 * (w - w0) ** 2


def doppler(v):
    """Right-mover frequency factor on a boosted inertial worldline: D=gamma(1-v)."""
    g = 1.0 / np.sqrt(1.0 - v * v)
    return g * (1.0 - v)


# ---- inertial-boost amplitude A(Omega): eta = -w D tau  (right-mover, origin)
def amplitude_boost(Omega, eps, S_kind, v, w0=W0, n_w=None, n_tau=None):
    D = doppler(v)
    base = int(np.ceil(1.0 / np.sqrt(eps)))
    n_w = n_w or min(2400, max(300, 40 * base))
    n_tau = n_tau or min(4000, max(400, 60 * base))
    wlo, whi = w0 - 6 * S_SIGMA, w0 + 6 * S_SIGMA
    w, ww = gl(n_w, wlo, whi)
    tau, wt = gl(n_tau, -(CW + CE), CW + CE)
    a = profile(w); S = eikonal_S(w, S_kind)
    amp = a * np.exp(1j * S / eps) * ww
    c = chi(tau)
    # inner psi_B(tau) = sum_w amp e^{-i w D tau}
    psi = np.exp(-1j * D * np.outer(tau, w)) @ amp
    A = np.sum(wt * c * np.exp(1j * Omega * tau) * psi)
    return -1j * A


def amplitude_mpmath(Omega, eps, S_kind, v, w0=W0):
    import mpmath as mp
    mp.mp.dps = 18
    D = doppler(v); sig = mp.mpf(S_SIGMA)
    wlo, whi = w0 - 6 * S_SIGMA, w0 + 6 * S_SIGMA
    a0, a1 = CW - CE, CW + CE

    def chi_mp(t):
        ax = abs(t)
        if ax <= a0: return mp.mpf(1)
        if ax >= a1: return mp.mpf(0)
        u = (ax - a0) / (a1 - a0)
        f0 = mp.e ** (-1 / (1 - u)); f1 = mp.e ** (-1 / u)
        return f0 / (f0 + f1)

    def S_mp(w):
        return (w - w0) if S_kind == "linear" else mp.mpf("0.5") * (w - w0) ** 2

    def inner(t):
        def f(w):
            a = mp.e ** (mp.mpf("-0.5") * ((w - w0) / sig) ** 2)
            return a * mp.e ** (1j * (S_mp(w) / eps - w * D * t))
        return mp.quad(f, [wlo, whi])

    def g(t):
        return chi_mp(t) * mp.e ** (1j * Omega * t) * inner(t)

    return complex(-1j * mp.quad(g, [-(CW + CE), CW + CE]))


# ------------------------------------------------------------------ control
def run_control(rows):
    print("[control] numeric quadrature vs mpmath ...")
    grid = [(e, s) for e in [0.1, 0.05, 0.02] for s in ["linear", "quadratic"]][:N_CONTROL]
    max_rel = 0.0; npass = 0
    for eps, S_kind in grid:
        Om = doppler(0.0) * W0
        ref = abs(amplitude_mpmath(Om, eps, S_kind, 0.0))
        got = abs(amplitude_boost(Om, eps, S_kind, 0.0))
        rel = abs(got - ref) / max(ref, 1e-300); absd = abs(got - ref)
        if ref < NEARZERO:        # non-stationary near-zero tail: judge absolutely
            ok = absd <= MPMATH_RTOL; crit = "abs"
        else:
            ok = rel <= MPMATH_RTOL; crit = "rel"; max_rel = max(max_rel, rel)
        npass += int(ok)
        rows.append({"block": "control", "eps": eps, "S": S_kind, "ref_abs": ref,
                     "num_abs": got, "rel_err": rel, "abs_err": absd,
                     "criterion": crit, "pass": int(ok)})
        print(f"   eps={eps} {S_kind:9s} ref={ref:.3e} num={got:.3e} "
              f"rel={rel:.2e} [{crit}] {'PASS' if ok else 'CHECK'}")
    print(f"[control] max rel (stationary) = {max_rel:.2e}; "
          f"points pass = {npass}/{len(grid)}")
    return max_rel, npass, len(grid)


# ---------------------------------------------------------------------- V1
def run_V1(rows):
    print("[V1] off-incidence localization ...")
    # detune Omega far above the populated band: no frequency-matching tau in supp chi
    Om = 2.6 * doppler(0.0) * W0
    vals = []
    for eps in EPS_MESH:
        A = amplitude_boost(Om, eps, "quadratic", 0.0)
        vals.append((eps, abs(A) ** 2))
        rows.append({"block": "V1", "eps": eps, "absA2": abs(A) ** 2, "Omega": Om})
    exps = []
    for i in range(1, len(vals)):
        (e0, a0), (e1, a1) = vals[i - 1], vals[i]
        if a0 > 0 and a1 > 0:
            p = (np.log(a1) - np.log(a0)) / (np.log(e1) - np.log(e0))
            exps.append((e1, p))
            rows.append({"block": "V1_exp", "eps": e1, "p_eff": p})
    ps = [p for _, p in exps]
    grows = bool(len(ps) >= 2 and (ps[-1] - ps[0]) >= 1.0 and
                 sum(ps[i] > ps[i-1]-1e-9 for i in range(1, len(ps))) >= len(ps)-2)
    print("   p_eff:", ", ".join(f"{e:.3g}:{p:.2f}" for e, p in exps),
          "-> grows:", grows)
    return exps, grows


# ---------------------------------------------------------------------- V2
def eta_accel(w, tau):
    # accelerated worldline: eta = (w/a)(e^{-a tau}-1); d2_tau eta|_* = a*Omega != 0
    return (w[None, :] / ACCEL) * (np.exp(-ACCEL * np.asarray(tau)[:, None]) - 1.0)


def run_V2(rows):
    print("[V2] sqrt(eps) coherent temporal width (accelerated worldline) ...")
    Om = 0.8
    widths = []
    for eps in EPS_MESH:
        n_w = min(3000, max(600, int(60 / np.sqrt(eps))))
        w, ww = gl(n_w, W0 - 6 * 0.30, W0 + 6 * 0.30)   # broad profile for V2
        a = np.exp(-0.5 * ((w - W0) / 0.30) ** 2)
        S = 0.5 * (w - W0) ** 2
        taus = np.linspace(-6, 6, 20001)
        et = eta_accel(w, taus)
        g = np.exp(1j * ((Om * taus[:, None] + et + S[None, :]) / eps)) @ (a * ww)
        ph = np.unwrap(np.angle(g)); f = np.gradient(ph, taus)
        i0 = int(np.argmin(np.abs(f)))
        dev = ph - ph[i0] - f[i0] * (taus - taus[i0])
        m = np.abs(dev) < np.pi
        lo = i0
        while lo > 0 and m[lo - 1]: lo -= 1
        hi = i0
        while hi < len(taus) - 1 and m[hi + 1]: hi += 1
        wdt = taus[hi] - taus[lo]
        widths.append((eps, wdt))
        rows.append({"block": "V2", "eps": eps, "width": wdt, "tau_star": taus[i0]})
    e = np.array([x[0] for x in widths]); wd = np.array([x[1] for x in widths])
    q, b = np.polyfit(np.log(e), np.log(wd), 1)
    resid = np.log(wd) - (q * np.log(e) + b); qs = []
    for _ in range(4000):
        rs = resid[rng.integers(0, len(resid), len(resid))]
        qq, _ = np.polyfit(np.log(e), (q * np.log(e) + b) + rs, 1); qs.append(qq)
    qs = np.sort(qs); lo, hi = qs[100], qs[3899]
    for r in rows:
        if r["block"] == "V2":
            r["q"] = q; r["q_lo"] = lo; r["q_hi"] = hi
    contains = bool(lo <= Q_TARGET <= hi)
    print(f"   q={q:.3f}  95% CI [{lo:.3f},{hi:.3f}]  contains 0.5: {contains}")
    return q, (lo, hi), contains


# ---------------------------------------------------------------------- V3
def run_V3(rows, eps=0.01):
    print("[V3] net Jacobian power beta (THE GATE) -- exact broadband Parseval ...")
    # Per-event broadband weight = flat-Omega sum of |A|^2 = 2pi eps int dtau |psi_B|^2
    # (Parseval). For the isolated event we take the full packet (window wide enough),
    # which is the per-event weight prop:stationary-weight is about, free of the
    # window-truncation artifact. P([Z]) ~ |phi|^2 J^beta ; J = w0 D.
    n_w = 1400
    w, ww = gl(n_w, W0 - 6 * S_SIGMA, W0 + 6 * S_SIGMA)
    a = profile(w); S = eikonal_S(w, "quadratic")
    amp = a * np.exp(1j * S / eps) * ww
    norm = np.sum(a ** 2 * ww)
    Js, Ws = [], []
    for v in V_LIST:
        D = doppler(v)
        half = max(CW + CE, 8 * (6.0 / S_SIGMA) / D)   # contain full packet
        taus = np.linspace(-half, half, 16000)
        psi = np.exp(-1j * D * np.outer(taus, w)) @ amp
        W = np.trapezoid(np.abs(psi) ** 2, taus)
        Js.append(W0 * D); Ws.append(W)
        rows.append({"block": "V3", "v": v, "J_AB": W0 * D, "weight": W,
                     "a2": 1.0, "eps": eps})
    Js = np.array(Js); Ws = np.array(Ws)
    y = np.log(Ws / norm); x = np.log(Js)
    beta, c = np.polyfit(x, y, 1)
    resid = y - (beta * x + c); bs = []
    for _ in range(4000):
        rs = resid[rng.integers(0, len(resid), len(resid))]
        bb, _ = np.polyfit(x, (beta * x + c) + rs, 1); bs.append(bb)
    bs = np.sort(bs); lo, hi = bs[100], bs[3899]
    for r in rows:
        if r["block"] == "V3":
            r["beta"] = beta; r["beta_lo"] = lo; r["beta_hi"] = hi
    gate = bool(abs(beta - 1.0) <= BETA_GATE)
    # invariance check: W*J/|phi|^2 should be constant (== 2 pi) iff beta == -1
    inv = Ws / norm * Js / W0
    print(f"   beta={beta:.4f}  95% CI [{lo:.3f},{hi:.3f}]  |beta-1|<= {BETA_GATE}: {gate}")
    print(f"   W*J/|phi|^2 (const 2pi <=> beta=-1): {np.array2string(inv, precision=3)}")
    return beta, (lo, hi), gate


# ---------------------------------------------------------------------- V3X
def run_V3X(rows, eps=0.01):
    # EXPLORATORY, NOT preregistered: derivative (field-strength) coupling weights
    # the integrand by w (supplies the sqrt(omega)-per-photon the scalar proxy of
    # rem:scalar-proxy drops). Tests whether the J^2 deficit lives in the proxy.
    print("[V3X] exploratory: derivative-coupling beta (outside preregistered battery) ...")
    n_w = 1400
    w, ww = gl(n_w, W0 - 6 * S_SIGMA, W0 + 6 * S_SIGMA)
    a = profile(w); S = eikonal_S(w, "quadratic")
    norm = np.sum(a ** 2 * ww)                          # normalize by |phi|^2
    Js, Ws = [], []
    for v in V_LIST:
        D = doppler(v)
        # proper-time derivative of the pullback: d/dtau e^{-i w D tau} = -i (w D),
        # i.e. the coupling sees the Doppler-shifted (detector-frame) frequency w D,
        # supplying the sqrt(omega)-per-photon the scalar proxy of rem:scalar-proxy
        # drops. amp carries the extra (w D) factor.
        amp = (w * D * a) * np.exp(1j * S / eps) * ww
        half = max(CW + CE, 8 * (6.0 / S_SIGMA) / D)
        taus = np.linspace(-half, half, 16000)
        psi = np.exp(-1j * D * np.outer(taus, w)) @ amp
        Ws.append(np.trapezoid(np.abs(psi) ** 2, taus)); Js.append(W0 * D)
        rows.append({"block": "V3X", "v": v, "J_AB": W0 * D,
                     "weight": Ws[-1], "eps": eps})
    Js = np.array(Js); Ws = np.array(Ws)
    betaX, _ = np.polyfit(np.log(Js), np.log(Ws / norm), 1)
    for r in rows:
        if r["block"] == "V3X":
            r["betaX"] = betaX
    print(f"   beta_X (derivative coupling) = {betaX:.4f}")
    return betaX


# ---------------------------------------------------------------------- V4
def run_V4(rows):
    print("[V4] non-interference between two separated events ...")
    Om = 0.8
    ratios = []
    for eps in EPS_MESH:
        n_w = min(3000, max(600, int(60 / np.sqrt(eps))))
        # two spectral lobes -> two arrival events on accelerated worldline
        w, ww = gl(n_w, 0.5, 1.6)
        a1 = np.exp(-0.5 * ((w - 0.8) / 0.07) ** 2)
        a2 = np.exp(-0.5 * ((w - 1.25) / 0.07) ** 2)
        S = 8.0 * 0.5 * (w - 1.0) ** 2          # strong chirp -> separated arrivals
        taus = np.linspace(-6, 6, 8000)
        et = eta_accel(w, taus)
        base = np.exp(1j * ((Om * taus[:, None] + et + S[None, :]) / eps))
        c2 = chi(taus) ** 2
        psi1 = base @ (a1 * ww); psi2 = base @ (a2 * ww)
        W1 = np.trapezoid(c2 * np.abs(psi1) ** 2, taus)
        W2 = np.trapezoid(c2 * np.abs(psi2) ** 2, taus)
        Wc = np.trapezoid(c2 * np.abs(psi1 + psi2) ** 2, taus)
        cross = Wc - (W1 + W2); ratio = abs(cross) / max(W1 + W2, 1e-300)
        ratios.append((eps, ratio))
        rows.append({"block": "V4", "eps": eps, "diag": W1 + W2,
                     "cross": cross, "ratio": ratio})
    decays = all(ratios[i][1] <= ratios[i-1][1] + 1e-9 for i in range(1, len(ratios)))
    print("   cross/diag:", ", ".join(f"{e:.3g}:{r:.1e}" for e, r in ratios),
          "-> decays:", decays)
    return ratios, decays


# --------------------------------------------------------------------- figs
def make_figures(v1_exps, v2_rows, v3_rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLU, ORG = "#2b6cb0", "#dd6b20"
    # V1
    fig, ax = plt.subplots(figsize=(5, 3.6))
    if v1_exps:
        e = [x[0] for x in v1_exps]; p = [x[1] for x in v1_exps]
        ax.plot(e, p, "o-", color=BLU); ax.set_xscale("log")
        ax.set_xlabel(r"$\epsilon$"); ax.set_ylabel(r"$p_{\rm eff}$")
        ax.set_title("V1: off-incidence super-algebraic decay"); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/b3b_fig_localizacion.png", dpi=140); plt.close(fig)
    # V2
    fig, ax = plt.subplots(figsize=(5, 3.6))
    r = [x for x in v2_rows if x["block"] == "V2"]
    e = np.array([x["eps"] for x in r]); wd = np.array([x["width"] for x in r])
    q = r[0]["q"]; b = np.log(wd).mean() - q * np.log(e).mean()
    ax.loglog(e, wd, "o", color=BLU, label="coherent width")
    xx = np.array([e.min(), e.max()])
    ax.loglog(xx, np.exp(q * np.log(xx) + b), "-", color=ORG, label=fr"fit $q={q:.3f}$")
    ax.loglog(xx, wd[0] * (xx / e[0]) ** 0.5, ":", color="gray", label=r"$q=1/2$")
    ax.set_xlabel(r"$\epsilon$"); ax.set_ylabel(r"$w(\epsilon)$")
    ax.set_title(r"V2: $O(\sqrt{\epsilon})$ temporal width"); ax.legend(fontsize=8)
    ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/b3b_fig_anchura.png", dpi=140); plt.close(fig)
    # V3
    fig, ax = plt.subplots(figsize=(5, 3.6))
    r = [x for x in v3_rows if x["block"] == "V3"]
    J = np.array([x["J_AB"] for x in r]); W = np.array([x["weight"] for x in r])
    beta = r[0]["beta"]; c = np.log(W).mean() - beta * np.log(J).mean()
    ax.loglog(J, W, "o", color=BLU, label="broadband weight")
    xx = np.array([J.min(), J.max()])
    ax.loglog(xx, np.exp(beta * np.log(xx) + c), "-", color=ORG, label=fr"fit $\beta={beta:.3f}$")
    ax.loglog(xx, W[0] * (xx / J[0]) ** 1.0, ":", color="green", label=r"asserted $\beta=+1$")
    ax.set_xlabel(r"$J_{AB}$"); ax.set_ylabel(r"$P/|a(k_Z)|^2$")
    ax.set_title(r"V3: net Jacobian power $\beta$"); ax.legend(fontsize=8)
    ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/b3b_fig_jacobiano.png", dpi=140); plt.close(fig)
    print("[fig] wrote b3b_fig_localizacion.png, b3b_fig_anchura.png, b3b_fig_jacobiano.png")


# --------------------------------------------------------------------- main
def main():
    t0 = time.time()
    print("=" * 70)
    print(f"SIM-V3-1  stationary-phase weight validation (seed {SEED})")
    print("=" * 70)
    rows = []
    max_rel, cpass, ctot = run_control(rows)
    v1_exps, v1_grows = run_V1(rows)
    q, (qlo, qhi), v2_pass = run_V2(rows)
    beta, (blo, bhi), v3_pass = run_V3(rows, eps=0.01)
    betaX = run_V3X(rows, eps=0.01)
    v4_ratios, v4_decays = run_V4(rows)

    # CSV
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(f"{OUTDIR}/b3b_eikonal_results.csv", "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=keys); wri.writeheader()
        for r in rows: wri.writerow(r)
    make_figures(v1_exps, rows, rows)

    gates = [
        ("V1_localization", "p_eff grows w/o saturation",
         "grows" if v1_grows else "saturates", int(v1_grows)),
        ("V2_width", "q=0.5 in 95% CI",
         f"q={q:.3f} CI[{qlo:.3f},{qhi:.3f}]", int(v2_pass)),
        ("V3_jacobian", "|beta-1|<=0.05",
         f"beta={beta:.4f} CI[{blo:.3f},{bhi:.3f}]", int(v3_pass)),
        ("V4_noninterf", "cross/diag decays",
         "decays" if v4_decays else "persists", int(v4_decays)),
        ("control_mpmath", "control points pass",
         f"{cpass}/{ctot}; max_rel={max_rel:.1e}", int(cpass == ctot)),
    ]
    with open(f"{OUTDIR}/b3b_verdict.csv", "w", newline="") as f:
        wri = csv.writer(f); wri.writerow(["gate", "criterion", "measured", "pass"])
        for g in gates: wri.writerow(g)
        wri.writerow(["V3X_exploratory", "derivative-coupling beta (not preregistered)",
                      f"betaX={betaX:.4f}", ""])

    print("=" * 70)
    overall = all(g[3] for g in gates)
    print("VERDICT:", "PASS" if overall else "DOES NOT PASS AS A WHOLE (see table)")
    for g in gates:
        print(f"  {g[0]:18s} {'PASS' if g[3] else 'FAIL':4s}  {g[2]}")
    print(f"  V3X_exploratory      ----  betaX={betaX:.4f} (derivative coupling)")
    print("exclusions:", EXCLUSIONS if EXCLUSIONS else "none")
    print(f"elapsed: {time.time()-t0:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
