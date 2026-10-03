"""Independent analytic model: Gaussian beams, mode overlap and 1-dB tolerances.

This module is the project's independent check on OpticStudio. It must never
import from the Zemax path (coupling_analysis, ZOS-API, pythonnet): only the
standard library, numpy and scipy. Every POP result is compared against it.

Contents
  Gaussian-beam ABCD propagation (complex q), refraction at spherical surfaces
  numerical power overlap of a curved Gaussian with a flat Gaussian receiver
  closed-form A0 overlap (Stage 5) and the B-family prediction (Stage 9/10)
  tolerance extraction from sampled curves (spline + brentq) and by root-finding
"""
import math

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

ONE_DB = math.log(10) / 10          # ln(10^0.1): eta falls 1 dB when ln(eta) drops by this
LAM = 1.31                          # um
N_SI = 3.5039127043661025           # OpticStudio SILICON_1310 at 1.31 um (Stage 7)
W0 = 4.6                            # um, grating-coupler output waist (Mangal 2021)
ZR_SI = math.pi * W0 ** 2 * N_SI / LAM
IDEAL_PRODUCT = ONE_DB * LAM / math.pi  # 0.096015 um rad: exact floor of lateral x angular (matched flat Gaussians; "0.0733 lambda" rounded)


def loss_db(eta):
    return -10 * np.log10(eta)


# ===================================================================== Gaussian beams
def q_param(w, R, n, lam):
    """Complex beam parameter from 1/e^2 radius w, wavefront radius R, index n."""
    return 1 / ((0 if math.isinf(R) else 1 / R) - 1j * lam / (math.pi * n * w * w))


def w_R(q, n, lam):
    iq = 1 / q
    return math.sqrt(-lam / (math.pi * n * iq.imag)), (math.inf if abs(iq.real) < 1e-18 else 1 / iq.real)


def refract(q, n1, n2, radius):
    """Refraction at a spherical surface (OpticStudio sign convention for radius)."""
    c = 0 if math.isinf(radius) else (n1 - n2) / (radius * n2)
    return q / (c * q + n1 / n2)


def beam_at(surfaces, gap_um, n_after, w0, lam, gap_surface):
    """Pilot-beam (w, R) in um at the last surface for a sequential surface list.

    surfaces: run_config "surfaces" list. The source waist sits on surface 1 and
    propagates through each thickness; refraction at each later surface uses
    n_after[i], the index of the medium after surface i.
    """
    q = q_param(w0, math.inf, n_after[1], lam)
    for i in range(1, len(surfaces) - 1):
        if i > 1:
            r = surfaces[i].get("radius_mm", "inf")
            q = refract(q, n_after[i - 1], n_after[i], math.inf if r == "inf" else float(r) * 1000)
        t = gap_um if i == gap_surface else float(surfaces[i]["thickness_mm"]) * 1000
        q = q + t
    return w_R(q, n_after[len(surfaces) - 2], lam)


def overlap(w1, R1, w2, lam, d=0.0, theta_deg=0.0, n=1.0, half=None, npts=40001):
    """Numerical power overlap of a curved Gaussian (w1, R1) with a flat Gaussian
    receiver (w2) decentred by d and tilted by theta. Separable in x and y; the
    y factor has neither offset nor tilt."""
    half = half or 8 * max(w1, w2) + abs(d)
    x = np.linspace(-half, half, npts)
    k = 2 * math.pi * n / lam
    curv = 0 if math.isinf(R1) else k / (2 * R1)

    def one(dd, th):
        e1 = np.exp(-x ** 2 / w1 ** 2 - 1j * curv * x ** 2)
        e2 = np.exp(-(x - dd) ** 2 / w2 ** 2 + 1j * k * math.sin(math.radians(th)) * x)
        num = abs(np.trapezoid(e1 * np.conj(e2), x)) ** 2
        return num / (np.trapezoid(abs(e1) ** 2, x) * np.trapezoid(abs(e2) ** 2, x))
    return one(d, theta_deg) * one(0.0, 0.0)


# ====================================================================== A0 closed form
def a0_closed_form(gap_um, w0=W0, w2=W0, lam=LAM):
    """A0: waist w0 after gap_um of air onto a flat receiver w2 (Stage 5).

    eta(d) = eta0 exp(-c d^2). Returns beam radius w1 and wavefront R at the
    fibre, eta0, c, the exact 1-dB lateral tolerance (incl. curvature), its
    absolute-convention value and the simple 0.339 sqrt(w1^2 + w2^2) estimate.
    """
    zr = math.pi * w0 ** 2 / lam
    w1 = w0 * math.sqrt(1 + (gap_um / zr) ** 2)
    R1 = math.inf if gap_um == 0 else gap_um * (1 + (zr / gap_um) ** 2)
    a1 = 1 / w1 ** 2 + (0 if math.isinf(R1) else 1j * (2 * math.pi / lam) / (2 * R1))
    a2 = 1 / w2 ** 2
    s = a1 + a2
    eta0 = 4 / (w1 ** 2 * w2 ** 2 * abs(s) ** 2)
    c = 2 * (a1 * a2 / s).real
    return dict(w1=w1, R1=R1, eta0=eta0, c=c, d_exact=math.sqrt(ONE_DB / c),
                d_abs=math.sqrt(math.log(eta0 / 10 ** -0.1) / c) if eta0 > 10 ** -0.1 else float("nan"),
                d_simple=0.339 * math.sqrt(w1 ** 2 + w2 ** 2))


# ============================================================ backside (B family) model
def r_match_um(t_um, n=N_SI):
    """Lens radius magnitude that collimates the Gaussian after t_um of silicon."""
    return (n - 1) / n * t_um * (1 + (ZR_SI / t_um) ** 2)


def w_exit_um(t_um):
    return W0 * math.sqrt(1 + (t_um / ZR_SI) ** 2)


def gauss_prediction(t_um=630.0, radius_um=-480.0, gap_um=20.0, mfd_um=34.0, n_si=N_SI):
    """ABCD + overlap prediction of loss and all tolerances (aberration-free Gaussian physics)."""
    surf = [{"thickness_mm": "inf"}, {"thickness_mm": t_um / 1000}, {"radius_mm": radius_um / 1000}, {}]
    n_after = {0: 1.0, 1: n_si, 2: 1.0}
    wf = mfd_um / 2

    def beam(g):
        return beam_at(surf, g, n_after, W0, LAM, 2)

    w, R = beam(gap_um)
    lat = lambda d: float(loss_db(overlap(w, R, wf, LAM, d=d)))
    ang = lambda a: float(loss_db(overlap(w, R, wf, LAM, theta_deg=a)))
    lon = lambda dz: float(loss_db(overlap(*beam(gap_um + dz), wf, LAM)))
    l0 = lat(0.0)
    out = dict(loss_dB=l0, w_um=w, R_um=R)
    for conv, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        for name, f, s0, hmax in (("lateral", lat, 5.0, 200.0), ("angular", ang, 0.4, 20.0),
                                  ("longitudinal", lon, 200.0, 20000.0)):
            out["%s_%s" % (name, conv)] = level_crossing(f, level, s0, hmax, l0) if level > l0 else float("nan")
        out["product_" + conv] = out["lateral_" + conv] * math.radians(out["angular_" + conv])
    return out


def aperture_prediction(diam_um, w_um=16.935, wf=17.0):
    """Hard circular aperture on a Gaussian of radius w, receiver wf, diffraction over the gap ignored."""
    a = diam_um / 2
    c = 1 / w_um ** 2 + 1 / wf ** 2
    eta = 4 * (1 - math.exp(-c * a * a)) ** 2 / (c * c * w_um ** 2 * wf ** 2)
    return dict(loss_dB=float(loss_db(eta)), S=1 - math.exp(-2 * a * a / w_um ** 2))


# ============================================================ tolerance extraction
def level_crossing(f, level, step0, hi_max, f0, xtol=1e-4):
    """Smallest d > 0 with f(d) = level, for f rising (eventually) from f(0) = f0 < level.
    Used for root-finding directly on a model or on POPD."""
    lo, hi = 0.0, step0
    fhi = f(hi)
    while fhi < level:
        lo, hi = hi, hi * 1.6
        if hi > hi_max:
            return float("nan")
        fhi = f(hi)
    return float(brentq(lambda d: f(d) - level, lo, hi, xtol=xtol))


def tolerance(x, loss, x0):
    """1-dB tolerance either side of x0 from a sampled curve, both conventions.

    rel: loss rises 1 dB above its value at x0.   abs: total loss reaches 1.000 dB.
    Cubic spline + brentq. Returns distances from x0 (positive), nan where the level
    is not reached inside the sweep or, for abs, where L(x0) is already >= 1 dB.
    """
    x, loss = np.asarray(x, float), np.asarray(loss, float)
    order = np.argsort(x)
    x, loss = x[order], loss[order]
    spl = CubicSpline(x, loss)
    l0 = float(spl(x0))
    out = dict(x0=float(x0), L0=l0)
    for name, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        for side in ("plus", "minus"):
            out[name + "_" + side] = float("nan")
            if level <= l0:
                continue
            xs = x[x > x0 + 1e-12] if side == "plus" else x[x < x0 - 1e-12][::-1]
            hit = np.nonzero(spl(xs) >= level)[0]
            if len(hit):
                j = hit[0]
                prev = x0 if j == 0 else xs[j - 1]
                root = brentq(lambda v: float(spl(v)) - level, min(prev, xs[j]), max(prev, xs[j]), xtol=1e-10)
                out[name + "_" + side] = abs(root - x0)
        vals = [out[name + "_plus"], out[name + "_minus"]]
        out[name] = float(np.nanmean(vals)) if not all(math.isnan(v) for v in vals) else float("nan")
    return out


def peak(x, loss):
    """Interpolated position and value of minimum loss."""
    x, loss = np.asarray(x, float), np.asarray(loss, float)
    order = np.argsort(x)
    x, loss = x[order], loss[order]
    i = int(np.argmin(loss))
    spl = CubicSpline(x, loss)
    lo, hi = x[max(i - 1, 0)], x[min(i + 1, len(x) - 1)]
    if lo == hi:
        return float(x[i]), float(loss[i])
    r = minimize_scalar(lambda v: float(spl(v)), bounds=(lo, hi), method="bounded", options=dict(xatol=1e-9))
    return float(r.x), float(r.fun)


def selftest():
    """The extractor must recover the analytic A0 tolerance from analytic samples (gap 20 um)."""
    a = a0_closed_form(20.0)
    d = np.round(np.arange(-6, 6.0001, 0.1), 6)
    got = tolerance(d, loss_db(a["eta0"] * np.exp(-a["c"] * d ** 2)), 0.0)["rel"]
    if abs(got - a["d_exact"]) > 1e-4:
        raise RuntimeError("tolerance extractor self-test failed: %.6f vs %.6f" % (got, a["d_exact"]))
    return got
