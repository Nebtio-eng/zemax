"""Stage 9: one-at-a-time parameter study of configuration B (backside micro-lens).

For every parameter value the script records nominal POPD 0/1/2 and the 1-dB
lateral, angular and longitudinal tolerances in both conventions (rel: 1 dB above
the on-axis value; abs: total loss 1.000 dB), plus the lateral x angular product.
Tolerances are found by root-finding directly on POPD (brentq, 1e-4 um / 1e-4 deg),
not by sampling. Each point is also predicted independently (ABCD + numerical
Gaussian overlap) where the physics is Gaussian.

Sweeps: thickness with the radius re-matched (receiver fixed, and receiver matched),
receiver MFD, lens radius alone (plus a radius x gap loss surface), gap, aperture,
conic constant, incidence angle. Convergence is spot-checked at both ends of
every sweep.

Usage:
    python parameter_study.py --mode standalone [--only receiver_mfd,conic]
"""
import argparse, csv, json, math, time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize

import coupling_analysis as C
from coupling_analysis import ROOT, PopSession, SettingsMismatch, loss_db

B_CONFIG = ROOT / "zemax" / "microlens" / "B" / "run_config.json"
OUT = ROOT / "results" / "parameter_studies"
LAM = 1.31
N_SI = 3.5039127043661025           # OpticStudio SILICON_1310 at 1.31 um (Stage 7)
ZR_SI = math.pi * 4.6 ** 2 * N_SI / LAM
IDEAL_PRODUCT = 0.0733 * LAM        # um rad, minimum for matched flat Gaussian modes
INK, INK2, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834"


def log(m):
    print(m, flush=True)
    with open(OUT / "progress.log", "a") as f:
        f.write(time.strftime("%H:%M:%S ") + m + "\n")


def r_match_um(t_um, n=N_SI):
    """Lens radius magnitude that collimates the Gaussian after t_um of silicon."""
    return (n - 1) / n * t_um * (1 + (ZR_SI / t_um) ** 2)


def w_exit_um(t_um):
    return 4.6 * math.sqrt(1 + (t_um / ZR_SI) ** 2)


# ============================================================================ analytics
def gauss_prediction(t_um=630.0, radius_um=-480.0, gap_um=20.0, mfd_um=34.0, n_si=N_SI):
    """ABCD + overlap prediction of loss and all tolerances (aberration-free Gaussian physics)."""
    surf = [{"thickness_mm": "inf"}, {"thickness_mm": t_um / 1000}, {"radius_mm": radius_um / 1000}, {}]
    n_after = {0: 1.0, 1: n_si, 2: 1.0}
    wf = mfd_um / 2

    def beam(g):
        return C.beam_at(surf, g, n_after, 4.6, LAM, 2)

    w, R = beam(gap_um)
    lat = lambda d: float(loss_db(C.overlap(w, R, wf, LAM, d=d)))
    ang = lambda a: float(loss_db(C.overlap(w, R, wf, LAM, theta_deg=a)))
    lon = lambda dz: float(loss_db(C.overlap(*beam(gap_um + dz), wf, LAM)))
    l0 = lat(0.0)
    out = dict(loss_dB=l0, w_um=w, R_um=R)
    for conv, level in (("rel", l0 + 1.0), ("abs", 1.0)):
        for name, f, s0, hmax in (("lateral", lat, 5.0, 200.0), ("angular", ang, 0.4, 20.0), ("longitudinal", lon, 200.0, 20000.0)):
            out["%s_%s" % (name, conv)] = level_crossing(f, level, s0, hmax, l0) if level > l0 else float("nan")
        out["product_" + conv] = out["lateral_" + conv] * math.radians(out["angular_" + conv])
    return out


def aperture_prediction(diam_um, w_um=16.935, wf=17.0):
    """Hard circular aperture on a Gaussian of radius w, matched receiver, diffraction over the gap ignored."""
    a = diam_um / 2
    c = 1 / w_um ** 2 + 1 / wf ** 2
    eta = 4 * (1 - math.exp(-c * a * a)) ** 2 / (c * c * w_um ** 2 * wf ** 2)
    return dict(loss_dB=float(loss_db(eta)), S=1 - math.exp(-2 * a * a / w_um ** 2))


# =========================================================================== root finding
def level_crossing(f, level, step0, hi_max, f0, xtol=1e-4):
    """Smallest d > 0 with f(d) = level, for f rising (eventually) from f(0) = f0 < level."""
    lo, hi = 0.0, step0
    fhi = f(hi)
    while fhi < level:
        lo, hi = hi, hi * 1.6
        if hi > hi_max:
            return float("nan")
        fhi = f(hi)
    return float(brentq(lambda d: f(d) - level, lo, hi, xtol=xtol))


# ================================================================================ bench
class Bench:
    """Configuration B in OpticStudio with setters for every studied parameter."""

    def __init__(self, app, zos, cfgd=None, incidence=False):
        self.app, self.zos = app, zos
        self.cfgd = cfgd or json.loads(B_CONFIG.read_text())
        self.cfgd = dict(self.cfgd)
        if incidence:
            self._add_incidence_cb()
        self.ps = PopSession(app, self.cfgd, zos)
        self.lde = self.ps.sys.LDE
        self.lens = 3 if incidence else 2
        self.gap0 = self.cfgd["parameters"]["gap_um"]
        self.x0 = self.th0 = 0.0
        self.lat_token = "POP_FPARAM4" if incidence else "POP_FPARAM3"

    def _add_incidence_cb(self):
        """Variant with a coordinate break (tilt about X) between the silicon and the lens."""
        c = json.loads(json.dumps(self.cfgd))
        c["surfaces"] = c["surfaces"][:2] + [{"comment": "CB: beam incidence on lens", "type": "CoordinateBreak",
                                              "thickness_mm": 0.0}] + c["surfaces"][2:]
        c["gap_surface"] = 3
        c["pop_tokens"]["POP_END"] = 4
        c["surface_settings"] = {"_comment": c["surface_settings"].get("_comment", ""),
                                 "1": c["surface_settings"]["1"], "3": c["surface_settings"]["2"]}
        c["id"], c["model"] = "B_incidence", "zemax/microlens/B_incidence/B_incidence.zmx"
        C.build_model(self.app, c, self.zos)
        self.cfgd = c

    # -- setters -------------------------------------------------------------------
    def set_radius_um(self, r):
        s = self.lde.GetSurfaceAt(self.lens)
        s.Radius = float("inf") if math.isinf(r) else r / 1000
        got = s.Radius
        if (math.isinf(r) and abs(got) < 1e9) or (not math.isinf(r) and abs(got - r / 1000) > 1e-12):
            raise SettingsMismatch("lens radius did not read back: %r" % got)

    def set_thickness_um(self, t):
        self.ps.set_thickness_um(1, t)

    def set_gap_um(self, g):
        self.ps.set_gap_um(g)

    def set_receiver_mfd(self, mfd):
        for k in ("POP_FPARAM1", "POP_FPARAM2"):
            self.ps.set(k, round(mfd / 2, 6) / 1000)

    def set_conic(self, k):
        s = self.lde.GetSurfaceAt(self.lens)
        s.Conic = k
        if s.Conic != k:
            raise SettingsMismatch("conic did not read back")

    def set_aperture(self, diam_um):
        ad = self.lde.GetSurfaceAt(self.lens).ApertureData
        T = self.zos.Editors.LDE.SurfaceApertureTypes
        if diam_um is None:
            ad.ChangeApertureTypeSettings(ad.CreateApertureTypeSettings(getattr(T, "None")))
            return
        st = ad.CreateApertureTypeSettings(T.CircularAperture)
        st._S_CircularAperture.MinimumRadius = 0.0
        st._S_CircularAperture.MaximumRadius = diam_um / 2000
        ad.ChangeApertureTypeSettings(st)
        got = ad.CurrentTypeSettings._S_CircularAperture.MaximumRadius
        if abs(got - diam_um / 2000) > 1e-12:
            raise SettingsMismatch("aperture radius did not read back: %r" % got)

    def set_incidence_deg(self, a):
        cb = self.lde.GetSurfaceAt(2)
        cell = cb.GetSurfaceCell(self.zos.Editors.LDE.SurfaceColumn.Par3)      # tilt about X
        cell.DoubleValue = a
        if abs(cb.GetSurfaceCell(self.zos.Editors.LDE.SurfaceColumn.Par3).DoubleValue - a) > 1e-12:
            raise SettingsMismatch("incidence tilt did not read back")

    # -- evaluation ----------------------------------------------------------------
    def popd(self, x=None, th=None, gap=None):
        self.ps.set(self.lat_token, (self.x0 if x is None else x) / 1000)
        self.ps.set("POP_TILTX", self.th0 if th is None else th)
        self.set_gap_um(self.gap0 if gap is None else gap)
        return self.ps.popd()

    def loss(self, **kw):
        return float(loss_db(self.popd(**kw)[0]))

    def metrics(self, guess=None, both_sides=False, longitudinal=True):
        """Nominal POPD and 1-dB tolerances (rel and abs) by root-finding on POPD."""
        guess = guess or {}
        e = self.popd()
        l0 = float(loss_db(e[0]))
        out = dict(eta=e[0], S=e[1], T=e[2], loss_dB=l0, fiber_x_um=self.x0, fiber_tilt_deg=self.th0)
        axes = [("lateral", lambda d, s: self.loss(x=self.x0 + s * d), guess.get("lateral", 8.0), 80.0),
                ("angular", lambda a, s: self.loss(th=self.th0 + s * a), guess.get("angular", 0.7), 10.0)]
        if longitudinal:
            axes.append(("longitudinal", lambda dz, s: self.loss(gap=self.gap0 + dz), guess.get("longitudinal", 700.0), 8000.0))
        for conv, level in (("rel", l0 + 1.0), ("abs", 1.0)):
            for name, f, g0, hmax in axes:
                sides = (1, -1) if (both_sides and name != "longitudinal") else (1,)
                vals = [level_crossing(lambda v: f(v, s), level, 0.8 * g0, hmax, l0) if level > l0 else float("nan")
                        for s in sides]
                out["%s_%s" % (name, conv)] = float(np.mean(vals))
                if len(vals) > 1:
                    out["%s_%s_sides" % (name, conv)] = vals
            out["product_" + conv] = out["lateral_" + conv] * math.radians(out["angular_" + conv])
        self.popd()
        return out

    def convergence(self, label, extra=None):
        """Grid x2 and grid x2 + window x2 at nominal and at one off-axis point."""
        s0, w0 = self.ps.tokens["POP_SAMPX"], self.ps.tokens["POP_WIDEX"]
        res = {}
        for pname, kw in (("nominal", {}), ("off-axis", extra or dict(x=self.x0 + 8.0))):
            etas = {}
            for name, s_, w_ in (("base", s0, w0), ("grid x2", s0 + 1, w0), ("grid x2, window x2", s0 + 1, 2 * w0)):
                for k in ("POP_SAMPX", "POP_SAMPY"):
                    self.ps.set(k, s_)
                for k in ("POP_WIDEX", "POP_WIDEY"):
                    self.ps.set(k, w_)
                etas[name] = self.popd(**kw)[0]
            res[pname] = {k: 100 * (v - etas["base"]) / etas["base"] for k, v in etas.items() if k != "base"}
        for k in ("POP_SAMPX", "POP_SAMPY"):
            self.ps.set(k, s0)
        for k in ("POP_WIDEX", "POP_WIDEY"):
            self.ps.set(k, w0)
        self.popd()
        worst = max(abs(v) for d in res.values() for v in d.values())
        log("  convergence @ %s: %s (worst %.4f%%)" % (label, res, worst))
        return dict(at=label, change_pct=res, worst_pct=worst, passed=worst < 0.5)

    def restore(self):
        p = self.cfgd
        self.x0 = self.th0 = 0.0
        self.set_radius_um(float(p["surfaces"][self.lens]["radius_mm"]) * 1000)
        self.set_thickness_um(float(p["surfaces"][1]["thickness_mm"]) * 1000)
        self.set_receiver_mfd(p["parameters"]["receiver_mfd_um"])
        self.set_conic(0.0)
        self.set_aperture(None)
        self.popd()

    def close(self):
        self.ps.close()


# ================================================================================ sweeps
FIELDS = ["value", "loss_dB", "eta", "S", "T", "lateral_rel", "lateral_abs", "angular_rel", "angular_abs",
          "longitudinal_rel", "longitudinal_abs", "product_rel", "product_abs"]


def run_sweep(b, name, values, apply, predict=None, guess_from_pred=True, unit="", note="", both_sides=False,
              longitudinal=True, extreme_extra=None):
    """Apply each value, measure metrics, compare with prediction, spot-check convergence at both ends."""
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    rows, preds, conv = [], [], []
    t0 = time.time()
    for i, v in enumerate(values):
        b.restore()
        apply(v)
        p = predict(v) if predict else {}
        guess = {k: p.get(k + "_rel") for k in ("lateral", "angular", "longitudinal")
                 if guess_from_pred and p.get(k + "_rel") and not math.isnan(p.get(k + "_rel"))}
        m = b.metrics(guess=guess, both_sides=both_sides, longitudinal=longitudinal)
        m["value"] = v
        rows.append(m)
        preds.append(p)
        log("%s = %s%s: loss %.4f dB (pred %s) | lat %.3f/%.3f um | ang %.4f/%.4f deg | long %s/%s um | prod %.5f  [%.0fs]" % (
            name, v, unit, m["loss_dB"], ("%.4f" % p["loss_dB"]) if "loss_dB" in p else "-",
            m["lateral_rel"], m["lateral_abs"], m["angular_rel"], m["angular_abs"],
            "%.1f" % m.get("longitudinal_rel", float("nan")), "%.1f" % m.get("longitudinal_abs", float("nan")),
            m["product_rel"], time.time() - t0))
        if i in (0, len(values) - 1):
            conv.append(b.convergence("%s = %s%s" % (name, v, unit),
                                      extra=extreme_extra(m) if extreme_extra else dict(x=b.x0 + m["lateral_rel"])))
    b.restore()
    with open(d / (name + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS + ["pred_" + k for k in FIELDS[1:] if k not in ("eta", "S", "T")] + ["fiber_x_um", "fiber_tilt_deg"])
        for m, p in zip(rows, preds):
            w.writerow([m.get(k, "") for k in FIELDS] +
                       [p.get(k, "") for k in FIELDS[1:] if k not in ("eta", "S", "T")] + [m["fiber_x_um"], m["fiber_tilt_deg"]])
    cfg = dict(stage="9", sweep=name, unit=unit, note=note, base_config=str(B_CONFIG.relative_to(ROOT)),
               values=values, opticstudio=dict(build=str(b.app.OpticStudioVersion), license=str(b.app.LicenseStatus)),
               pop_tokens=b.ps.tokens, surface_settings=b.cfgd.get("surface_settings"),
               method="tolerances by brentq on POPD (xtol 1e-4); rel = +1 dB from on-axis, abs = 1.000 dB total",
               settings_write_retries=getattr(b.ps, "retries", 0),
               convergence=conv, settings_readbacks_passed=b.ps.nchecks, results=rows, predictions=preds,
               wall_time_s=time.time() - t0)
    (d / "run_config.json").write_text(json.dumps(cfg, indent=2, default=str))
    plot_sweep(name, values, rows, preds, unit, d / (name + ".png"))
    if not all(c["passed"] for c in conv):
        raise RuntimeError("%s: convergence failed at an extreme: %s" % (name, conv))
    return rows, preds, conv


def run_radius_gap_surface(b, radii, gaps):
    """Loss over (lens radius, gap) and the radius error that costs 1 dB at each gap."""
    d = OUT / "radius_alone"
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    b.restore()
    grid = np.zeros((len(gaps), len(radii)))
    for j, r in enumerate(radii):
        b.set_radius_um(-r)
        for i, g in enumerate(gaps):
            grid[i, j] = b.loss(gap=g)
    tol = []
    for g in gaps:
        b.set_radius_um(-480.0)
        l0 = b.loss(gap=g)
        row = dict(gap_um=g, loss_at_480_dB=l0)
        for conv, level in (("rel", l0 + 1.0), ("abs", 1.0)):
            for side, s in (("smaller", -1), ("larger", 1)):
                def f(dr, s=s):
                    b.set_radius_um(-(480.0 + s * dr))
                    return b.loss(gap=g)
                row["%s_%s_um" % (conv, side)] = level_crossing(f, level, 20.0, 470.0 if s < 0 else 4000.0, l0) \
                    if level > l0 else float("nan")
        tol.append(row)
        log("radius tolerance @ gap %g um: L(480) %.3f dB; rel -%.1f/+%.1f um; abs -%s/+%s um  [%.0fs]" % (
            g, l0, row["rel_smaller_um"], row["rel_larger_um"], round(row["abs_smaller_um"], 1),
            round(row["abs_larger_um"], 1), time.time() - t0))
    b.restore()
    with open(d / "radius_gap_loss_surface.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gap_um \\ radius_um"] + radii)
        for g, row in zip(gaps, grid):
            w.writerow([g] + ["%.6f" % v for v in row])
    with open(d / "radius_tolerance_vs_gap.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(tol[0].keys()))
        w.writeheader()
        w.writerows(tol)
    plot_surface(radii, gaps, grid, tol, d / "radius_gap_loss_surface.png")
    return grid, tol


def optimise_fiber(b, guess_x, guess_th):
    """Re-point the receiver (X decenter, tilt) for best coupling: Nelder-Mead on POPD."""
    f = lambda p: -b.popd(x=p[0], th=p[1])[0]
    r = minimize(f, [guess_x, guess_th], method="Nelder-Mead",
                 options=dict(xatol=1e-3, fatol=1e-9, initial_simplex=[[guess_x, guess_th], [guess_x + 2, guess_th],
                                                                       [guess_x, guess_th + 0.3]]))
    return float(r.x[0]), float(r.x[1]), -float(r.fun), int(r.nfev)


# ================================================================================ plots
def _style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)


def plot_sweep(name, values, rows, preds, unit, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array(values, float)
    panels = [("loss_dB", "Nominal loss (dB)"), ("lateral", "Lateral 1-dB (um)"), ("angular", "Angular 1-dB (deg)"),
              ("longitudinal", "Longitudinal 1-dB (um)"), ("product", "Lateral x angular (um rad)")]
    fig, axes = plt.subplots(1, 5, figsize=(16, 3.4), dpi=150)
    for ax, (key, title) in zip(axes, panels):
        if key == "loss_dB":
            series = [("loss_dB", "Zemax", BLUE, "-")]
        else:
            series = [(key + "_rel", "Zemax rel", BLUE, "-"), (key + "_abs", "Zemax abs", ORANGE, "-")]
        for k, lab, col, ls in series:
            y = np.array([r.get(k, np.nan) for r in rows], float)
            ax.plot(x, y, color=col, lw=1.8, marker="o", ms=4, label=lab)
            yp = np.array([p.get(k, np.nan) for p in preds], float)
            if np.isfinite(yp).any():
                ax.plot(x, yp, color=INK2, lw=1, ls="--", label="predicted" if k.endswith(("rel", "dB")) else None)
        if key == "product":
            ax.axhline(IDEAL_PRODUCT, color=INK2, lw=0.8, ls=":")
            ax.annotate("minimum 0.0960", (x.min(), IDEAL_PRODUCT), fontsize=7, color=INK2, va="bottom")
        ax.set_title(title, loc="left", fontsize=9, color=INK)
        ax.set_xlabel("%s %s" % (name.replace("_", " "), "(%s)" % unit.strip() if unit.strip() else ""), fontsize=8, color=INK)
        _style(ax)
    axes[1].legend(frameon=False, fontsize=7)
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Stage 9, configuration B: %s" % name.replace("_", " "), x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_surface(radii, gaps, grid, tol, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), dpi=150)
    ax = axes[0]
    im = ax.pcolormesh(radii, gaps, grid, cmap="Blues", vmin=0, vmax=6, shading="nearest")
    cs = ax.contour(radii, gaps, grid, levels=[1.0, 3.0], colors=[ORANGE, INK], linewidths=[2, 1])
    ax.clabel(cs, fmt="%.0f dB", fontsize=8)
    ax.axvline(480, color=INK2, ls="--", lw=1)
    ax.set_xlabel("Lens radius |R| (um), substrate fixed at 630 um", color=INK)
    ax.set_ylabel("Air gap (um)", color=INK)
    ax.set_title("Coupling loss (dB)", loc="left", fontsize=10, color=INK)
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax = axes[1]
    g = [t["gap_um"] for t in tol]
    ax.plot(g, [t["rel_smaller_um"] for t in tol], color=BLUE, marker="o", ms=4, label="rel, radius smaller")
    ax.plot(g, [t["rel_larger_um"] for t in tol], color=BLUE, ls="--", marker="o", ms=4, label="rel, radius larger")
    ax.plot(g, [t["abs_smaller_um"] for t in tol], color=ORANGE, marker="s", ms=4, label="abs, radius smaller")
    ax.plot(g, [t["abs_larger_um"] for t in tol], color=ORANGE, ls="--", marker="s", ms=4, label="abs, radius larger")
    ax.set_yscale("log")
    ax.set_xlabel("Air gap (um)", color=INK)
    ax.set_ylabel("Radius error for 1 dB (um)", color=INK)
    ax.set_title("Lens-radius manufacturing tolerance vs gap", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ================================================================================ driver
def run_all(app, zos, only=None):
    OUT.mkdir(parents=True, exist_ok=True)
    want = lambda n: only is None or n in only
    summary = {}
    b = Bench(app, zos)
    try:
        # validation of the root-finder against Stage 7 sweep values
        m = b.metrics()
        log("root-finder check vs Stage 7: lat %.4f (8.1399) ang %.4f (0.6758) long %.1f (713.3) loss %.4f (0.0004)" % (
            m["lateral_rel"], m["angular_rel"], m["longitudinal_rel"], m["loss_dB"]))
        for k, ref, tol in (("lateral_rel", 8.1399, 1e-3), ("angular_rel", 0.6758, 1e-3), ("longitudinal_rel", 713.3, 0.2)):
            if abs(m[k] - ref) > tol:
                raise RuntimeError("root-finder disagrees with Stage 7 on %s: %.5f vs %.5f" % (k, m[k], ref))
        summary["baseline"] = m
        if want("thickness_rematched"):
            ts = [400.0, 450.0, 500.0, 550.0, 592.0, 630.0, 700.0, 800.0, 900.0, 1000.0]

            def ap(t):
                b.set_thickness_um(t)
                b.set_radius_um(-r_match_um(t))
            summary["thickness_rematched"] = run_sweep(
                b, "thickness_rematched", ts, ap, lambda t: gauss_prediction(t, -r_match_um(t)), unit=" um",
                note="lens radius re-matched for collimation at each thickness; receiver fixed at 34 um MFD")
        if want("thickness_rematched_receiver_matched"):
            ts = [400.0, 450.0, 500.0, 550.0, 592.0, 630.0, 700.0, 800.0, 900.0, 1000.0]

            def ap(t):
                b.set_thickness_um(t)
                b.set_radius_um(-r_match_um(t))

            def ap2(t):
                ap(t)
                b.set_receiver_mfd(2 * w_exit_um(t))
            summary["thickness_rematched_receiver_matched"] = run_sweep(
                b, "thickness_rematched_receiver_matched", ts, ap2,
                lambda t: gauss_prediction(t, -r_match_um(t), mfd_um=2 * w_exit_um(t)), unit=" um",
                note="lens radius re-matched AND receiver MFD = 2 x beam radius at each thickness (scaling law)")
        if want("receiver_mfd"):
            mfds = [20.0, 24.0, 28.0, 30.0, 32.0, 34.0, 38.0, 42.0, 48.0]
            summary["receiver_mfd"] = run_sweep(b, "receiver_mfd", mfds, b.set_receiver_mfd,
                                                lambda v: gauss_prediction(mfd_um=v), unit=" um")
        if want("radius_alone"):
            radii = [300.0, 340.0, 380.0, 420.0, 450.0, 480.0, 510.0, 550.0, 600.0, 650.0, 700.0]
            summary["radius_alone"] = run_sweep(b, "radius_alone", radii, lambda r: b.set_radius_um(-r),
                                                lambda r: gauss_prediction(radius_um=-r), unit=" um",
                                                note="substrate fixed at 630 um: de-collimation, the manufacturing tolerance")
            summary["radius_gap_surface"] = run_radius_gap_surface(
                b, [float(r) for r in range(300, 801, 25)], [20.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 850.0, 1000.0])
        if want("gap"):
            gaps = [0.0, 20.0, 50.0, 100.0, 200.0, 400.0, 700.0]

            def apg(g):
                b.gap0 = g
            try:
                summary["gap"] = run_sweep(b, "gap", gaps, apg, lambda g: gauss_prediction(gap_um=g), unit=" um",
                                           note="longitudinal tolerance measured from each gap, gap increasing")
            finally:
                b.gap0 = b.cfgd["parameters"]["gap_um"]
        if want("conic"):
            ks = [-10.0, -3.0, -1.0, 0.0, 1.0, 3.0, 10.0]
            summary["conic"] = run_sweep(b, "conic", ks, b.set_conic, lambda k: gauss_prediction(), unit="",
                                         note="prediction is the k = 0 Gaussian value: conic changes sag by ~k*0.0015 um at 2w")
        if want("aperture"):
            ds = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 102.0, 150.0, 200.0]
            summary["aperture"] = run_sweep(b, "aperture", ds, b.set_aperture, aperture_prediction, unit=" um",
                                            guess_from_pred=False, note="hard circular aperture on the lens surface (diameter)")
        b.restore()
    finally:
        b.close()
    if want("incidence"):
        summary["incidence"] = run_incidence(app, zos)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def run_incidence(app, zos, angles=(0.0, 1.0, 2.0, 3.0, 4.0, 6.0)):
    """Beam incidence on the lens (angle in silicon). Si path set to 630/cos(a); receiver re-pointed."""
    b = Bench(app, zos, incidence=True)
    rows = []
    try:
        def ap(a):
            b.set_incidence_deg(a)
            b.set_thickness_um(630.0 / math.cos(math.radians(a)))
            b.x0, b.th0 = 0.0, 0.0
            unrepointed = b.popd()[0]
            exit_deg = math.degrees(math.asin(min(0.99, N_SI * math.sin(math.radians(a)))))
            # POP defines the receiver relative to the beam's chief ray (verified: at 4 deg the optimum fiber tilt is
            # -0.003 deg although the beam leaves at 14.2 deg), so the optimiser starts on-axis and only removes the
            # small aberration-induced shift. The physical fibre tilt relative to the chip normal is the Snell angle.
            x, th, eta, nfev = optimise_fiber(b, 0.0, 0.0)
            b.x0, b.th0 = x, th
            rows.append(dict(angle_deg=a, exit_angle_deg_snell=exit_deg, eta_fiber_on_chief_ray=unrepointed,
                             eta_repointed=eta, fiber_y_um=x, fiber_tilt_deg=th, nfev=nfev))
            log("incidence %.1f deg (Snell exit %.2f deg): eta on chief ray %.6f; optimised eta %.6f at y %.3f um, "
                "tilt %.4f deg (%d evals)" % (a, exit_deg, unrepointed, eta, x, th, nfev))

        res = run_sweep(b, "incidence", list(angles), ap, None, unit=" deg", both_sides=True, longitudinal=False,
                        note="coordinate break tilts the lens relative to the beam; Si path 630/cos(a); "
                             "receiver re-pointed (X decenter + tilt, Nelder-Mead) before measuring tolerances")
        (OUT / "incidence" / "repointing.json").write_text(json.dumps(rows, indent=2))
        return res
    finally:
        b.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    zos = C.load_zosapi()
    conn, app = C.connect(args.mode, zos)
    try:
        run_all(app, zos, only=set(args.only.split(",")) if args.only else None)
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
