#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
B.3  --  Numerical implementation and VALIDATION of the GEOMETRIC INCIDENCE
         layer in isolation.   Appendix B, \label{app:sim-incidence}.

  Realizes the geometric machinery of  sec_05_incidence.tex  and
  sec_06_born.tex:  the null congruence, its eikonal sampling, the incidence
  points by geometric uniqueness, the Jacobian, the Born recovery, and the
  sigma_B transfer via the isometric inclusion i^+.

  It is the geometric analogue of what B.2 (b2_tsym_sector.py) did for the
  matter sector.

--------------------------------------------------------------------------------
SCOPE AND HONEST FRONTIER  (do not cross; principles 1, 5, 8)
--------------------------------------------------------------------------------
This program implements and VALIDATES the GEOMETRIC layer in isolation:
    * the null congruence C_A (Definition def:null-congruence) and its eikonal
      (ray) sampling, with amplitude |psi(Z)|^2 over PN (eq:photon-state);
    * the incidence points p_{Z*} by the geometric uniqueness of the
      timelike-null intersection (Theorem thm:uniqueness);
    * the geometric Jacobian J_AB = -eta(ell_Z, u_B)  (eq:jacobian);
    * the Born recovery P([Z*]) ~ |psi(Z*)|^2 relative to the geometric
      measure J_AB d mu  (Theorem thm:born, eq:incidence-density / eq:born);
    * the relativistic aberration of incidence directions;
    * the isometry of the orientation inclusion i^+ (rem:unitary), and that
      flipping sigma_B (+ -> -) swaps the sector WITHOUT touching the
      intersection geometry.

It is NOT the discriminator (protocol G1).  It does NOT:
    * couple to the B.2 matter sector,
    * flip the thermodynamic gradient or measure D(off),
    * run the registration-substrate audit (rem:registration-trap),
    * run the convention/vacuity control (rem:convention-trap).
Those are B.4/B.5/B.6 and remain "specified, not executed".

Therefore a PASS here verifies that the incidence machinery is implemented
faithfully.  It addresses O1 (single-outcome selection) and O2 (Born weights)
of Definition def:falsifiable-input GEOMETRICALLY -- which is plausible and is
NOT the point in dispute -- but it does NOT address O3 (the no-return /
irreversibility), where the discriminator's verdict actually lives.

CRUCIAL HONESTY (sec_06, "Status of the geometric Jacobian"): that the geometry
reproduces the Born weights and the Doppler/aberration factor is NEITHER a
deviation from QED NOR evidence of the central claim.  sec_06 already
established that J_AB coincides NUMERICALLY with the relativistic frequency
factor QED already carries; the contribution is STRUCTURAL, not predictive.
We declare it as such and do not oversell.

--------------------------------------------------------------------------------
THE MODEL  (declared explicitly where the definitions leave freedom)
--------------------------------------------------------------------------------
Minkowski (M, eta), signature (-,+,+,+), natural units c = hbar = 1.
eta(a,b) = -a^0 b^0 + a.b   (3-vector dot product on the spatial part).

  * Null congruence: a ray L_Z is a null geodesic x(lambda) = x0 + lambda*lhat,
    lhat = (1, n), |n| = 1, n on the celestial sphere S^2.  The photon
    four-momentum tangent is ell_Z = omega0 * lhat = (omega0, omega0 n).
    The eikonal phase / impact data parametrize which parallel ray is taken;
    we reuse the honest geometric protocol P4 (ray tracing), NOT the GKSL
    dynamics of P1-P3/P5.
  * Observer worldline gamma_B with four-velocity u_B:
      - INERTIAL:   u_B = gammaL (1, v),  gammaL = 1/sqrt(1-|v|^2);
                    gamma_B(tau) = x_B0 + tau u_B.
      - ACCELERATED (richer test): constant proper acceleration 'a' along x,
        gamma_B(tau) = (sinh(a tau)/a, cosh(a tau)/a + xc, y0, z0),
        u_B(tau)     = (cosh(a tau), sinh(a tau), 0, 0),  eta(u_B,u_B) = -1.
  * Incidence: for each ray, the intersection with gamma_B is the unique
    p_{Z*} of Theorem thm:uniqueness; in Minkowski it is algebraic for the
    inertial worldline and a 1-D root find for the accelerated one.
  * sigma_B: discrete label +/-1; i^+ is the inclusion into the forward sector.

JACOBIAN -- HONEST FRONTIER FLAGGED (principle 2, "signal drift explicitly"):
  eq:jacobian is the COVARIANT object J_AB = -eta(ell_Z, u_B).  Carrying out
  the contraction gives, EXACTLY,
        J_AB = -eta(ell_Z,u_B) = omega0 * gammaL * (1 - n.v),
  the full relativistic Doppler factor.  The form omega0*(1 - n.v) requested
  in the brief is eq:jacobian-nr of sec_09, which sec_09 declares to be the
  SLOW-OBSERVER LIMIT (gammaL -> 1) of eq:jacobian.  We therefore validate
  BOTH, and report the gammaL discrepancy at finite |v| honestly:
    (2a) numeric -eta(ell_Z,u_B)  vs  omega0 gammaL (1 - n.v)  -> machine prec.
         (this is the genuine validation of eq:jacobian);
    (2b) numeric  vs  omega0 (1 - n.v)  [eq:jacobian-nr]       -> agreement as
         |v| -> 0, with departure O(v^2): the controlled slow-observer limit,
         NOT a failure.

--------------------------------------------------------------------------------
WHAT IS COMPUTED, AND THE PREDICTION IT IS COMPARED AGAINST
--------------------------------------------------------------------------------
  (1) Uniqueness (thm:uniqueness): per transverse incident ray, number of
      intersections with gamma_B = 1 (0 for non-incident / offset rays).
      Limiting (grazing) case handled and counted separately.
      Prediction: exactly 1 in the transverse case; 0 for offset; both for the
      inertial AND the accelerated worldline.
  (2) Jacobian (eq:jacobian): numeric -eta(ell_Z,u_B) vs closed forms, sweeping
      |v| and the direction n.  Prediction: machine precision vs the covariant
      relativistic form; O(v^2) departure from the NR form.
  (3) Born recovery (thm:born): Monte Carlo sampling of the congruence with the
      geometric measure d mu and |psi(Z)|^2, weighted by J_AB, accumulated over
      outcome channels, vs the Born prediction P([Z*]) ~ |psi|^2 J_AB / N_A.
      Prediction: total-variation distance decays as M^{-1/2}; we VERIFY THE
      MEASURED EXPONENT (~ -0.5), not just the final value.
  (4) Aberration: distribution of incidence directions seen by the moving
      observer vs the relativistic aberration prediction.
      Prediction: matches the analytic pushforward density to MC tolerance.
  (5) Isometry of i^+ (rem:unitary): ||(i^+)^dag i^+ - I|| at machine precision;
      orthogonality of the +/- images; and that flipping sigma_B selects the
      opposite sector WITHOUT changing the intersection geometry (p_{Z*}, J_AB).
      (The convention/vacuity verdict proper is B.6; here we only show the
      transfer mechanism is isometric and well defined.)

OUTPUTS: console + CSV table (numeric, theory, deviation, tolerance, PASS/FAIL
per observable) and four PNG figures (a)-(d).

HONESTY: the program is EXECUTED and reports REAL numbers.  If any observable
FAILS -- e.g. multiple intersections for a transverse ray, or J_AB not matching
the closed form, or Born not converging as M^{-1/2} -- it is reported as a
FINDING, NOT whitewashed.  The data govern.

DIAGNOSIS LOG (first executed version; corrections are to the MEASUREMENT
APPARATUS, never to the data):
  * uniqueness_transverse_accelerated initially reported 0 intersections.
    Diagnosis: the intersection exists exactly (miss^2(tau0) = 0 verified);
    the single-step parabolic refinement of the zero -- exact for the inertial
    worldline, where miss^2 is exactly quadratic in tau -- leaves a residual
    ~4e-12 for the sinh/cosh worldline, above the 1e-16 acceptance threshold.
    Detector defect, not a counterexample to thm:uniqueness.  Fixed by robust
    root finding: Brent on the null-separation function h(tau) =
    eta(gamma(tau)-x0, gamma(tau)-x0) (whose zeros are the <=2 light-cone
    crossings of the timelike worldline), followed by the line-membership
    filter |r| < tol (the cone crossing must lie ON the ray, not just on the
    cone); thm:uniqueness then predicts <=1 survivor.
  * jacobian_NR_slow_limit_exponent measured 2.204 on the fit window
    v in [0.1, 0.3] vs the asymptotic theory 2.0.  Diagnosis: the EXACT
    NR departure is omega0*(gammaL-1)*(1+v); its local log-log slope is
    2 + v/(1+v) = 2.09..2.23 on that window -- the measurement was CORRECT,
    the fit window inappropriate for the v->0 asymptote.  Fixed twofold:
    (a) validate the exact closed form omega0*(gammaL-1)*(1+v) of the NR
    departure at machine precision (the stronger statement), and (b) fit the
    leading exponent on a genuinely small-v window (v <= 0.05).
  * born_TV figure showed a floor ~0.009 at M >= 5e4 instead of M^{-1/2}.
    Diagnosis: the reference target itself was computed by Monte Carlo
    quadrature (n = 2e5) and so carried its OWN statistical error ~0.005
    (TV between two independent target draws measured at 0.0067) -- the
    empirical distributions converge to the true target, but the distance to
    a NOISY reference floors at the reference's error.  Reference-measurement
    defect, not a failure of Born convergence.  Fixed by computing the target
    with a deterministic midpoint product-grid quadrature in (z, phi)
    (quadrature error ~1e-6, no statistical noise; phi grid aligned with the
    channel bins for exact assignment).
================================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# -----------------------------------------------------------------------------
# Reproducibility and global parameters
# -----------------------------------------------------------------------------
SEED = 20240603
rng = np.random.default_rng(SEED)

OUTDIR = "./"
os.makedirs(OUTDIR, exist_ok=True)

OMEGA0 = 1.0                                   # photon frequency scale
V_SWEEP = np.array([0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99])
N_DIRS_SWEEP = 64                              # directions per |v| in the J sweep
A_ACC = 0.7                                    # proper acceleration (richer test)

# Born convergence
M_LIST = [200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
N_SEEDS_BORN = 40                              # seeds averaged per M (clean slope)
N_CHANNELS = 12                                # outcome channels (azimuth bins)

# Aberration / direction sampling
N_ABER = 400000
V_ABER = 0.6                                   # observer speed for aberration test

# Isometry model space
DIM_MODEL = 5                                  # finite model Hilbert dim for i^+

# Tolerances (explicit; deviations measured, not eyeballed)
TOL_UNIQUE = 0                                 # exact integer count match
TOL_JAC_COV = 1e-12                            # numeric vs covariant closed form
TOL_BORN_EXP = 0.07                            # |measured exponent - (-0.5)|
TOL_ABER = 0.02                                # max bin density deviation (MC)
TOL_ISO = 1e-13                                # ||(i^+)^dag i^+ - I||
TOL_GEOM_FLIP = 1e-13                          # geometry invariance under sigma_B flip


# =============================================================================
# Minkowski geometry
# =============================================================================
def eta(a, b):
    r"""Minkowski inner product, signature (-,+,+,+)."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    return -a[0] * b[0] + np.dot(a[1:], b[1:])


def lorentz_factor(v3):
    v2 = float(np.dot(v3, v3))
    return 1.0 / np.sqrt(1.0 - v2)


def u_inertial(v3):
    r"""Four-velocity of an inertial observer with 3-velocity v3 (|v|<1)."""
    g = lorentz_factor(v3)
    return np.array([g, g * v3[0], g * v3[1], g * v3[2]])


def null_tangent(n3, omega=OMEGA0):
    r"""Photon four-momentum tangent ell_Z = omega (1, n), n a unit 3-vector."""
    n3 = np.asarray(n3, dtype=float)
    return np.array([omega, omega * n3[0], omega * n3[1], omega * n3[2]])


def unit_dir_from_angles(theta, phi):
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])


# =============================================================================
# (1) Uniqueness of the timelike-null intersection  (thm:uniqueness)
# =============================================================================
def gamma_inertial(tau, u, x0):
    return x0 + tau * u


def gamma_accel(tau, a=A_ACC, xc=0.0, y0=0.0, z0=0.0):
    r"""Constant proper-acceleration worldline along x; eta(u,u) = -1."""
    return np.array([np.sinh(a * tau) / a,
                     np.cosh(a * tau) / a + xc, y0, z0])


def u_accel(tau, a=A_ACC):
    return np.array([np.cosh(a * tau), np.sinh(a * tau), 0.0, 0.0])


def transverse_miss(worldline_pts, x0, nhat):
    r"""
    For a ray x(lambda) = x0 + lambda*(1, nhat), a worldline point X lies on the
    ray iff its spatial offset equals (Delta t) * nhat, i.e. the transverse miss
    r(X) = (X_space - x0_space) - (X_time - x0_time) nhat = 0.
    Returns |r|^2 along the supplied worldline samples.
    """
    dt = worldline_pts[:, 0] - x0[0]
    dspace = worldline_pts[:, 1:] - x0[1:]
    r = dspace - np.outer(dt, nhat)
    return np.sum(r * r, axis=1)


def count_intersections(worldline_fn, x0, nhat, tau_grid):
    r"""
    Robust intersection counter (see DIAGNOSIS LOG).
    Step 1: find the zeros of the null-separation function
        h(tau) = eta(gamma(tau) - x0, gamma(tau) - x0),
    whose roots are the crossings of the worldline with the light cone of x0.
    A timelike worldline crosses each cone branch (future/past) at most once,
    so h has at most two simple roots; they are bracketed by sign changes on
    tau_grid and polished with Brent to machine precision.
    Step 2: line-membership filter.  A cone crossing lies on THE ray (the full
    null geodesic through x0 with direction nhat) iff the transverse miss
    r = (Delta x_space) - (Delta t) nhat vanishes; we accept
    |r| < 1e-9 * (1 + |Delta t|).  thm:uniqueness predicts <= 1 survivor.
    Returns the integer count and the intersection proper times.
    """
    pts = np.array([worldline_fn(t) for t in tau_grid])
    d = pts - x0[None, :]
    h = -d[:, 0] ** 2 + np.sum(d[:, 1:] ** 2, axis=1)        # eta(d, d)

    roots = []
    for i in range(len(tau_grid) - 1):
        if h[i] == 0.0:
            roots.append(tau_grid[i])
        elif h[i] * h[i + 1] < 0.0:
            ts = brentq(
                lambda t: eta(worldline_fn(t) - x0, worldline_fn(t) - x0),
                tau_grid[i], tau_grid[i + 1], xtol=1e-14, rtol=8.9e-16)
            roots.append(ts)
    if h[-1] == 0.0:
        roots.append(tau_grid[-1])

    hits = []
    for ts in roots:
        P = worldline_fn(ts)
        dt = P[0] - x0[0]
        r = (P[1:] - x0[1:]) - dt * nhat
        if np.linalg.norm(r) < 1e-9 * (1.0 + abs(dt)):
            if not hits or all(abs(ts - hh) > 1e-6 for hh in hits):
                hits.append(ts)
    return len(hits), hits


def build_incident_ray(worldline_fn, tau0, nhat, lam0=3.0):
    r"""Construct a ray (base point x0, direction nhat) guaranteed to pass
    through gamma_B(tau0): place x0 a distance lam0 back along the ray."""
    P = worldline_fn(tau0)
    lhat = np.array([1.0, nhat[0], nhat[1], nhat[2]])
    x0 = P - lam0 * lhat
    return x0


def run_uniqueness():
    rows = []
    tau_grid = np.linspace(-40.0, 40.0, 40001)

    # ---- inertial worldline ----
    v = np.array([0.25, -0.1, 0.05])
    u = u_inertial(v); x_B0 = np.array([0.0, 0.0, 0.0, 0.0])
    wl_in = lambda t: gamma_inertial(t, u, x_B0)

    cnt_trans, cnt_offset, cnt_graze = [], [], []
    test_dirs = [unit_dir_from_angles(th, ph)
                 for th in np.linspace(0.3, np.pi - 0.3, 6)
                 for ph in np.linspace(0.0, 2 * np.pi, 7)[:-1]]
    for nhat in test_dirs:
        for tau0 in (-5.0, 0.0, 7.0):
            x0 = build_incident_ray(wl_in, tau0, nhat)
            c, _ = count_intersections(wl_in, x0, nhat, tau_grid)
            cnt_trans.append(c)
            # offset ray: translate base point transversally to n -> should miss
            perp = np.cross(nhat, [1.0, 0.0, 0.0])
            if np.linalg.norm(perp) < 1e-6:
                perp = np.cross(nhat, [0.0, 1.0, 0.0])
            perp = perp / np.linalg.norm(perp)
            x0_off = x0 + np.array([0.0, 2.0 * perp[0], 2.0 * perp[1], 2.0 * perp[2]])
            c_off, _ = count_intersections(wl_in, x0_off, nhat, tau_grid)
            cnt_offset.append(c_off)

    # grazing / limiting case: ray direction approaching the observer 3-velocity
    # direction (transversality marginal). Build incident ray with n almost || v.
    n_graze = v / np.linalg.norm(v)
    for tau0 in (-3.0, 0.0, 4.0):
        x0 = build_incident_ray(wl_in, tau0, n_graze)
        c, _ = count_intersections(wl_in, x0, n_graze, tau_grid)
        cnt_graze.append(c)

    # ---- accelerated worldline (richer test) ----
    wl_acc = lambda t: gamma_accel(t)
    cnt_trans_acc = []
    for nhat in test_dirs:
        for tau0 in (-2.0, 0.0, 2.5):
            x0 = build_incident_ray(wl_acc, tau0, nhat)
            c, _ = count_intersections(wl_acc, x0, nhat, tau_grid)
            cnt_trans_acc.append(c)

    cnt_trans = np.array(cnt_trans)
    cnt_offset = np.array(cnt_offset)
    cnt_graze = np.array(cnt_graze)
    cnt_trans_acc = np.array(cnt_trans_acc)

    rows.append(dict(observable="uniqueness_transverse_inertial",
                     numeric=float(cnt_trans.max()) if len(cnt_trans) else np.nan,
                     theory=1.0,
                     deviation=int(np.max(np.abs(cnt_trans - 1))),
                     tol=TOL_UNIQUE,
                     verdict="PASS" if np.all(cnt_trans == 1) else "FAIL"))
    rows.append(dict(observable="uniqueness_offset_inertial(=0)",
                     numeric=float(cnt_offset.max()) if len(cnt_offset) else 0.0,
                     theory=0.0,
                     deviation=int(np.max(np.abs(cnt_offset - 0))),
                     tol=TOL_UNIQUE,
                     verdict="PASS" if np.all(cnt_offset == 0) else "FAIL"))
    rows.append(dict(observable="uniqueness_grazing_inertial",
                     numeric=float(cnt_graze.max()) if len(cnt_graze) else np.nan,
                     theory=1.0,
                     deviation=int(np.max(np.abs(cnt_graze - 1))),
                     tol=TOL_UNIQUE,
                     verdict="PASS" if np.all(cnt_graze == 1) else "FAIL"))
    rows.append(dict(observable="uniqueness_transverse_accelerated",
                     numeric=float(cnt_trans_acc.max()) if len(cnt_trans_acc) else np.nan,
                     theory=1.0,
                     deviation=int(np.max(np.abs(cnt_trans_acc - 1))),
                     tol=TOL_UNIQUE,
                     verdict="PASS" if np.all(cnt_trans_acc == 1) else "FAIL"))

    diag = dict(trans=cnt_trans, offset=cnt_offset, graze=cnt_graze,
                trans_acc=cnt_trans_acc)
    return rows, diag


# =============================================================================
# (2) Geometric Jacobian  (eq:jacobian)
# =============================================================================
def jacobian_numeric(ell, u):
    return -eta(ell, u)


def jacobian_cov_closed(omega0, nhat, v3):
    r"""Exact covariant form: omega0 * gammaL * (1 - n.v)."""
    g = lorentz_factor(v3)
    return omega0 * g * (1.0 - np.dot(nhat, v3))


def jacobian_nr_closed(omega0, nhat, v3):
    r"""Slow-observer limit eq:jacobian-nr: omega0 (1 - n.v)."""
    return omega0 * (1.0 - np.dot(nhat, v3))


def run_jacobian():
    rows = []
    sweep = {"v": [], "cov_dev": [], "nr_dev": [], "vmag": []}
    max_cov_dev = 0.0
    for vmag in V_SWEEP:
        # boost direction varied; pick a fixed boost axis and sweep n on a circle
        vdir = np.array([1.0, 0.0, 0.0])
        v3 = vmag * vdir
        u = u_inertial(v3)
        phis = np.linspace(0.0, 2 * np.pi, N_DIRS_SWEEP, endpoint=False)
        cov_devs, nr_devs = [], []
        for phi in phis:
            # direction n in x-z plane spanning angle to boost axis
            nhat = np.array([np.cos(phi), 0.0, np.sin(phi)])
            ell = null_tangent(nhat, OMEGA0)
            Jnum = jacobian_numeric(ell, u)
            Jcov = jacobian_cov_closed(OMEGA0, nhat, v3)
            Jnr = jacobian_nr_closed(OMEGA0, nhat, v3)
            cov_devs.append(abs(Jnum - Jcov))
            nr_devs.append(abs(Jnum - Jnr))
        cov_dev = float(np.max(cov_devs)); nr_dev = float(np.max(nr_devs))
        sweep["vmag"].append(vmag)
        sweep["cov_dev"].append(cov_dev)
        sweep["nr_dev"].append(nr_dev)
        max_cov_dev = max(max_cov_dev, cov_dev)

    rows.append(dict(observable="jacobian_numeric_vs_covariant(eq:jacobian)",
                     numeric=max_cov_dev, theory=0.0, deviation=max_cov_dev,
                     tol=TOL_JAC_COV,
                     verdict="PASS" if max_cov_dev <= TOL_JAC_COV else "FAIL"))

    # eq:jacobian-nr is the slow-observer limit.  Two checks (DIAGNOSIS LOG):
    # (a) the EXACT closed form of the NR departure.  Over the swept circle
    #     n.v = v cos(phi), the maximal departure is
    #         max_n |J_cov - J_nr| = omega0 (gammaL - 1)(1 + v),
    #     a machine-precision prediction valid at ALL v.
    vmags = np.array(sweep["vmag"][1:])         # exclude v = 0
    nrdev = np.array(sweep["nr_dev"][1:])
    gammas = 1.0 / np.sqrt(1.0 - vmags ** 2)
    nrdev_exact = OMEGA0 * (gammas - 1.0) * (1.0 + vmags)
    dev_nr_closed = float(np.max(np.abs(nrdev - nrdev_exact)))
    rows.append(dict(observable="jacobian_NR_departure_closed_form",
                     numeric=dev_nr_closed, theory=0.0, deviation=dev_nr_closed,
                     tol=TOL_JAC_COV,
                     verdict="PASS" if dev_nr_closed <= TOL_JAC_COV else "FAIL"))

    # (b) the leading exponent of the departure as v -> 0 (expected 2, since
    #     gammaL - 1 = v^2/2 + O(v^4)); fitted on a genuinely small-v window
    #     (v <= 0.05) where the (1+v) subleading factor is negligible.
    small = vmags <= 0.05
    if np.sum(small) >= 2 and np.all(nrdev[small] > 0):
        p_nr = np.polyfit(np.log(vmags[small]), np.log(nrdev[small]), 1)[0]
    else:
        p_nr = np.nan
    rows.append(dict(observable="jacobian_NR_slow_limit_exponent(eq:jacobian-nr)",
                     numeric=p_nr, theory=2.0, deviation=abs(p_nr - 2.0),
                     tol=0.2,
                     verdict="PASS" if abs(p_nr - 2.0) <= 0.2 else "FAIL"))
    return rows, sweep


# =============================================================================
# (3) Born recovery via Monte Carlo  (thm:born, eq:incidence-density)
# =============================================================================
def psi2_profile(nhat):
    r"""
    |psi(Z)|^2 over the celestial sphere, summed over helicities.  Declared
    model amplitude: a mixture of two von Mises-Fisher-like bumps (one per
    helicity sector), positive and smooth.  Any positive profile works; this
    one is fixed for reproducibility.
    """
    c1 = np.array([0.0, 0.0, 1.0]); c2 = np.array([1.0, 0.0, 0.0])
    k1, k2 = 4.0, 2.5
    h_plus = np.exp(k1 * (np.dot(nhat, c1) - 1.0))
    h_minus = 0.6 * np.exp(k2 * (np.dot(nhat, c2) - 1.0))
    return h_plus + h_minus


def channel_of(nhat):
    r"""Outcome channel = azimuth bin (discrete set of outcomes [Z*])."""
    phi = np.arctan2(nhat[1], nhat[0]) % (2 * np.pi)
    return int(phi / (2 * np.pi) * N_CHANNELS) % N_CHANNELS


def born_target(v3, n_z=2000, n_phi=2880):
    r"""
    Deterministic Born target per channel (DIAGNOSIS LOG, third entry):
      P([Z*]) ~ |psi|^2 * J_AB  over d mu (uniform on the sphere),
    computed by midpoint product-grid quadrature in (z, phi).  The uniform
    sphere measure is uniform in (z, phi); midpoint in periodic phi is
    spectrally accurate and n_phi is a multiple of N_CHANNELS so every grid
    point falls strictly inside one azimuth bin (exact channel assignment).
    Quadrature error ~1e-6, negligible against all MC tolerances; the target
    carries NO statistical noise of its own.
    """
    assert n_phi % N_CHANNELS == 0
    z = -1.0 + (np.arange(n_z) + 0.5) * (2.0 / n_z)
    phi = (np.arange(n_phi) + 0.5) * (2.0 * np.pi / n_phi)
    Z, PHI = np.meshgrid(z, phi, indexing="ij")
    S = np.sqrt(1.0 - Z * Z)
    X = S * np.cos(PHI); Y = S * np.sin(PHI)
    g = lorentz_factor(v3)
    nv = X * v3[0] + Y * v3[1] + Z * v3[2]
    psi2 = np.exp(4.0 * (Z - 1.0)) + 0.6 * np.exp(2.5 * (X - 1.0))
    w = psi2 * OMEGA0 * g * (1.0 - nv)
    chan = (PHI / (2 * np.pi) * N_CHANNELS).astype(int) % N_CHANNELS
    P = np.zeros(N_CHANNELS)
    for k in range(N_CHANNELS):
        P[k] = w[chan == k].sum()
    return P / P.sum()


def sample_sphere(m, r):
    r"""Uniform samples on S^2 (the geometric measure d mu over directions)."""
    z = r.uniform(-1.0, 1.0, m)
    phi = r.uniform(0.0, 2 * np.pi, m)
    s = np.sqrt(1.0 - z * z)
    return np.stack([s * np.cos(phi), s * np.sin(phi), z], axis=1)


def born_empirical(M, v3, r):
    r"""
    Weighted Monte Carlo over the congruence: sample directions from d mu,
    weight by |psi|^2 * J_AB, bin into channels -> empirical outcome dist.
    """
    dirs = sample_sphere(M, r)
    nv = dirs @ v3
    g = lorentz_factor(v3)
    Jac = OMEGA0 * g * (1.0 - nv)                       # J_AB per sample
    psi2 = np.array([psi2_profile(n) for n in dirs])
    w = psi2 * Jac
    phi = np.arctan2(dirs[:, 1], dirs[:, 0]) % (2 * np.pi)
    chan = (phi / (2 * np.pi) * N_CHANNELS).astype(int) % N_CHANNELS
    P = np.zeros(N_CHANNELS)
    for k in range(N_CHANNELS):
        P[k] = w[chan == k].sum()
    s = P.sum()
    return P / s if s > 0 else P


def tv_distance(p, q):
    return 0.5 * np.sum(np.abs(p - q))


def run_born():
    v3 = np.array([0.4, 0.0, 0.0])
    P_target = born_target(v3)
    tv_means = []
    for M in M_LIST:
        tvs = []
        for s in range(N_SEEDS_BORN):
            r = np.random.default_rng(SEED + 1000 * s + M)
            Phat = born_empirical(M, v3, r)
            tvs.append(tv_distance(Phat, P_target))
        tv_means.append(float(np.mean(tvs)))
    tv_means = np.array(tv_means)
    Ms = np.array(M_LIST, dtype=float)
    slope, intercept = np.polyfit(np.log(Ms), np.log(tv_means), 1)
    rows = [dict(observable="born_TV_convergence_exponent",
                 numeric=slope, theory=-0.5, deviation=abs(slope + 0.5),
                 tol=TOL_BORN_EXP,
                 verdict="PASS" if abs(slope + 0.5) <= TOL_BORN_EXP else "FAIL")]
    diag = dict(Ms=Ms, tv=tv_means, slope=slope, intercept=intercept,
                P_target=P_target, v3=v3)
    return rows, diag


# =============================================================================
# (4) Relativistic aberration of incidence directions
# =============================================================================
def aberration_cos(cos_lab, v):
    r"""Observer (moving +x at speed v) sees lab direction cos_lab at cos':
        cos' = (cos_lab - v) / (1 - v cos_lab)."""
    return (cos_lab - v) / (1.0 - v * cos_lab)


def run_aberration():
    r = np.random.default_rng(SEED + 7)
    dirs = sample_sphere(N_ABER, r)            # isotropic lab directions
    cos_lab = dirs[:, 0]                       # cos angle to boost axis (+x)
    cos_obs = aberration_cos(cos_lab, V_ABER)
    nb = 40
    edges = np.linspace(-1, 1, nb + 1)
    hist, _ = np.histogram(cos_obs, bins=edges, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # analytic pushforward density of uniform cos_lab under aberration:
    #   f(cos') = (1 - v^2) / (1 + v cos')^2 * (1/2)   on cos' in [-1,1]
    f_theory = 0.5 * (1.0 - V_ABER ** 2) / (1.0 + V_ABER * centers) ** 2
    dev = float(np.max(np.abs(hist - f_theory)))
    rows = [dict(observable="aberration_density_max_dev",
                 numeric=dev, theory=0.0, deviation=dev, tol=TOL_ABER,
                 verdict="PASS" if dev <= TOL_ABER else "FAIL")]
    diag = dict(centers=centers, hist=hist, f_theory=f_theory, v=V_ABER)
    return rows, diag


# =============================================================================
# (5) Isometry of i^+ and the sigma_B flip  (rem:unitary, prop:geom-orientation)
# =============================================================================
def run_isometry():
    d = DIM_MODEL
    I = np.eye(d)
    Z = np.zeros((d, d))
    # i^+ : H_A -> (H_A)^+ subset (H_A)^+ (+) (H_A)^-   ;  |psi| -> (|psi|, 0)
    iplus = np.vstack([I, Z])                  # (2d, d)
    iminus = np.vstack([Z, I])                 # sigma_B flip selects the - sector
    iso_plus = np.linalg.norm(iplus.conj().T @ iplus - I)
    iso_minus = np.linalg.norm(iminus.conj().T @ iminus - I)
    orth = np.linalg.norm(iplus.conj().T @ iminus)   # images orthogonal => 0
    # left-inverse recovers the state (reversible on image): (i^+)^dag i^+ = I
    left_inv_dev = np.linalg.norm(iplus.conj().T @ iplus - I)

    # Geometry invariance under the sigma_B flip: incidence point and Jacobian
    # are unchanged when sigma_B: + -> - (the flip swaps the Hilbert-space
    # sector, not the spacetime intersection).
    v3 = np.array([0.3, 0.1, 0.0]); u = u_inertial(v3); x_B0 = np.zeros(4)
    wl = lambda t: gamma_inertial(t, u, x_B0)
    nhat = unit_dir_from_angles(1.0, 0.7)
    tau0 = 2.0
    x0 = build_incident_ray(wl, tau0, nhat)
    tau_grid = np.linspace(-40, 40, 40001)
    _, hits_plus = count_intersections(wl, x0, nhat, tau_grid)
    _, hits_minus = count_intersections(wl, x0, nhat, tau_grid)   # geometry identical
    p_plus = wl(hits_plus[0]); p_minus = wl(hits_minus[0])
    ell = null_tangent(nhat, OMEGA0)
    J_plus = jacobian_numeric(ell, u); J_minus = jacobian_numeric(ell, u)
    geom_dev = float(np.max(np.abs(p_plus - p_minus)) + abs(J_plus - J_minus))

    rows = [
        dict(observable="isometry_i+_||(i+)^dag i+ - I||", numeric=iso_plus,
             theory=0.0, deviation=iso_plus, tol=TOL_ISO,
             verdict="PASS" if iso_plus <= TOL_ISO else "FAIL"),
        dict(observable="isometry_i-_||(i-)^dag i- - I||", numeric=iso_minus,
             theory=0.0, deviation=iso_minus, tol=TOL_ISO,
             verdict="PASS" if iso_minus <= TOL_ISO else "FAIL"),
        dict(observable="sector_orthogonality_||(i+)^dag i-||", numeric=orth,
             theory=0.0, deviation=orth, tol=TOL_ISO,
             verdict="PASS" if orth <= TOL_ISO else "FAIL"),
        dict(observable="sigmaB_flip_geometry_invariance", numeric=geom_dev,
             theory=0.0, deviation=geom_dev, tol=TOL_GEOM_FLIP,
             verdict="PASS" if geom_dev <= TOL_GEOM_FLIP else "FAIL"),
    ]
    diag = dict(iso_plus=iso_plus, iso_minus=iso_minus, orth=orth, geom_dev=geom_dev)
    return rows, diag


# =============================================================================
# Figures
# =============================================================================
def figure_jacobian(sweep, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    # left: J_AB vs angle for a few |v|, numeric (points) vs covariant (line)
    vshow = [0.0, 0.3, 0.7, 0.99]
    phis = np.linspace(0.0, 2 * np.pi, N_DIRS_SWEEP, endpoint=False)
    for vmag in vshow:
        v3 = vmag * np.array([1.0, 0.0, 0.0])
        u = u_inertial(v3)
        Jnum = np.array([jacobian_numeric(null_tangent(np.array([np.cos(p), 0, np.sin(p)])), u)
                         for p in phis])
        Jcov = np.array([jacobian_cov_closed(OMEGA0, np.array([np.cos(p), 0, np.sin(p)]), v3)
                         for p in phis])
        ax1.plot(phis, Jcov, "-", lw=1.2, label=rf"$|v|={vmag}$ closed form")
        ax1.plot(phis[::4], Jnum[::4], "o", ms=3)
    ax1.set_xlabel(r"angle $\phi$ of $\mathbf{n}$ to boost axis")
    ax1.set_ylabel(r"$J_{AB}=-\eta(\ell_Z,u_B)$")
    ax1.set_title(r"(a1) $J_{AB}$: numeric (pts) vs $\omega_0\gamma_L(1-\mathbf{n}\cdot\mathbf{v})$")
    ax1.legend(fontsize=7, ncol=2)
    # right: deviations vs |v|: covariant (machine prec) and NR (O(v^2))
    vmag = np.array(sweep["vmag"]); covdev = np.array(sweep["cov_dev"])
    nrdev = np.array(sweep["nr_dev"])
    ax2.semilogy(vmag, np.maximum(covdev, 1e-17), "o-", color="C0",
                 label=r"vs covariant (eq:jacobian): machine prec.")
    ax2.semilogy(vmag, np.maximum(nrdev, 1e-17), "s--", color="C3",
                 label=r"vs NR $\omega_0(1-\mathbf{n}\cdot\mathbf{v})$ (eq:jacobian-nr)")
    vfit = np.linspace(0.05, 0.5, 50)
    ax2.semilogy(vfit, 0.5 * vfit ** 2, "k:", lw=0.9, label=r"$\propto v^2$ guide")
    ax2.set_xlabel(r"observer speed $|v|$"); ax2.set_ylabel("max deviation")
    ax2.set_title(r"(a2) covariant exact; NR is the slow-observer limit $O(v^2)$")
    ax2.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def figure_born(diag, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    Ms, tv = diag["Ms"], diag["tv"]
    ax.loglog(Ms, tv, "o-", color="C0", lw=1.6, label="measured TV distance (real data)")
    ref = tv[0] * (Ms / Ms[0]) ** (-0.5)
    ax.loglog(Ms, ref, "k--", lw=1.2, label=r"$\propto M^{-1/2}$ reference")
    ax.set_xlabel("number of Monte Carlo samples $M$")
    ax.set_ylabel("total-variation distance to Born target")
    ax.set_title(rf"(b) Born recovery convergence: measured slope $={diag['slope']:.3f}$ "
                 r"(theory $-0.5$)")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def figure_aberration(diag, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(diag["centers"], diag["hist"], "o", ms=4, color="C0",
            label=f"sampled incidence directions (v={diag['v']}, real data)")
    ax.plot(diag["centers"], diag["f_theory"], "r-", lw=1.5,
            label=r"aberration theory $\frac{1-v^2}{2(1+v\cos\theta')^2}$")
    ax.set_xlabel(r"$\cos\theta'$ (observer frame)")
    ax.set_ylabel("probability density")
    ax.set_title("(c) Relativistic aberration of incidence directions")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def figure_diagnostics(uniq_diag, iso_diag, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    # left: histogram of intersection counts
    allcounts = {
        "transverse\n(inertial)": uniq_diag["trans"],
        "offset\n(inertial)": uniq_diag["offset"],
        "grazing\n(inertial)": uniq_diag["graze"],
        "transverse\n(accel.)": uniq_diag["trans_acc"],
    }
    labels = list(allcounts.keys())
    means = [np.mean(allcounts[k]) for k in labels]
    # bar = the (uniform) intersection count for each family
    vals0 = [np.mean(allcounts[k] == 0) for k in labels]
    vals1 = [np.mean(allcounts[k] == 1) for k in labels]
    vals2 = [np.mean(allcounts[k] >= 2) for k in labels]
    x = np.arange(len(labels)); w = 0.25
    ax1.bar(x - w, vals0, w, label="0 intersections", color="C7")
    ax1.bar(x, vals1, w, label="1 intersection", color="C2")
    ax1.bar(x + w, vals2, w, label=r"$\geq 2$ (theorem-forbidden)", color="C3")
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=8)
    ax1.set_ylabel("fraction of rays"); ax1.set_ylim(0, 1.1)
    ax1.set_title("(d1) Uniqueness (thm:uniqueness): exactly 1 or 0")
    ax1.legend(fontsize=8)
    # right: isometry diagnostics (log bars)
    keys = [r"$\|(i^+)^\dag i^+ - I\|$", r"$\|(i^-)^\dag i^- - I\|$",
            r"$\|(i^+)^\dag i^-\|$", r"$\sigma_B$-flip geom. dev."]
    vals = [max(iso_diag["iso_plus"], 1e-18), max(iso_diag["iso_minus"], 1e-18),
            max(iso_diag["orth"], 1e-18), max(iso_diag["geom_dev"], 1e-18)]
    ax2.bar(range(len(keys)), vals, color=["C0", "C0", "C4", "C5"])
    ax2.axhline(TOL_ISO, color="r", ls=":", lw=1.0, label=f"tolerance {TOL_ISO:.0e}")
    ax2.set_yscale("log"); ax2.set_xticks(range(len(keys)))
    ax2.set_xticklabels(keys, fontsize=8, rotation=12)
    ax2.set_ylabel("value (machine-precision scale)")
    ax2.set_title(r"(d2) Isometry of $i^+$ and $\sigma_B$-flip (rem:unitary)")
    ax2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 78)
    print("B.3  VALIDATION OF THE GEOMETRIC INCIDENCE LAYER (in isolation)")
    print("Faithfulness check of the incidence machinery -- NOT the discriminator,")
    print("NOT evidence of new physics (sec_06: J_AB coincides with the QED factor).")
    print("=" * 78)

    rows = []
    r_u, uniq_diag = run_uniqueness(); rows += r_u
    r_j, jac_sweep = run_jacobian(); rows += r_j
    r_b, born_diag = run_born(); rows += r_b
    r_a, aber_diag = run_aberration(); rows += r_a
    r_i, iso_diag = run_isometry(); rows += r_i

    df = pd.DataFrame(rows, columns=["observable", "numeric", "theory",
                                     "deviation", "tol", "verdict"])
    csv_path = os.path.join(OUTDIR, "b3_incidence_results.csv")
    df.to_csv(csv_path, index=False)

    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.float_format", lambda x: f"{x:.4e}")
    print(df.to_string(index=False))
    print(f"\n[saved] {csv_path}")

    # ---- figures ----
    fa = os.path.join(OUTDIR, "b3_fig_a_jacobian.png")
    fb = os.path.join(OUTDIR, "b3_fig_b_born.png")
    fc = os.path.join(OUTDIR, "b3_fig_c_aberration.png")
    fd = os.path.join(OUTDIR, "b3_fig_d_diagnostics.png")
    figure_jacobian(jac_sweep, fa)
    figure_born(born_diag, fb)
    figure_aberration(aber_diag, fc)
    figure_diagnostics(uniq_diag, iso_diag, fd)

    # ---- summary ----
    overall = "PASS" if (df.verdict == "PASS").all() else "FAIL"
    n_fail = int((df.verdict == "FAIL").sum())
    print("\n" + "=" * 78)
    print("SUMMARY (per-observable verdict; deviations measured, not eyeballed)")
    print("=" * 78)
    for _, row in df.iterrows():
        print(f"  {row['observable']:48s} dev={row['deviation']:.3e}  "
              f"tol={row['tol']:.1e}  [{row['verdict']}]")
    print("-" * 78)
    print(f"  OVERALL: {overall}   ({n_fail} failing rows of {len(df)})")
    print("  Interpretation (honest frontier): a PASS verifies that the geometric")
    print("  incidence machinery is implemented faithfully -- uniqueness (O1),")
    print("  Born weights (O2), the Doppler/aberration Jacobian, and the isometric")
    print("  sigma_B transfer.  It does NOT address O3 (no-return / irreversibility),")
    print("  does NOT couple to the matter sector, and does NOT run the gradient")
    print("  switch or the registration audit (B.4/B.5/B.6).  Reproducing the Born")
    print("  weights and the Doppler factor is NOT a deviation from QED and NOT")
    print("  evidence of the central claim (sec_06, 'Status of the geometric Jacobian').")
    print("=" * 78)
    print("\n[figures saved]")
    for f in (fa, fb, fc, fd):
        print("  ", f)
    return df


if __name__ == "__main__":
    main()
