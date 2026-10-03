"""1-D alignment sweeps with 1-dB tolerance extraction (Stages 5, 7 and 8).

Stage 5 (A0): receiver X decenter at several air gaps.            run(app)
Stage 7 (A1, B): lateral, angular and longitudinal sweeps, lens-sign test,
                 propagator comparison and convergence.              run_stage7(app)
Stage 8: A0 angular and longitudinal; B rebuilt on the 2021 pair.    run_stage8(app)

POPD 0, 1 and 2 are logged at every point; tolerances by spline interpolation
(analytic.tolerance), never the nearest sample, in both conventions.

Everything about the POP setup comes from zemax/baseline/run_config.json
("pop_tokens"). The POP .CFG is generated from that file in code, and every
setting is read back out of the CFG and checked; any mismatch raises.

Usage:
    python alignment_sweep.py --mode standalone      # long sweeps (headless)
    python alignment_sweep.py --mode extension       # needs Interactive Extension armed
or import and call run(app) with an already-connected application.
"""
import argparse, csv, json, math, time
from pathlib import Path

import numpy as np

import analytic as AN
import coupling_analysis as C
from analytic import loss_db, a0_closed_form, tolerance, selftest
from coupling_analysis import ROOT, SettingsMismatch, PopSession, load_zosapi, write_csv, frange

CONFIG = ROOT / "zemax" / "baseline" / "run_config.json"
MODEL = ROOT / "zemax" / "baseline" / "A0_conventional.zmx"
OUT = ROOT / "results" / "coupling_curves"


def log(m):
    C.progress("sweep", m)


# ------------------------------------------------------------------ OpticStudio side
def a0_report_checks(ps, gap_um):
    """Generic report checks (PopSession.report) plus the A0-specific ones."""
    rep = ps.report()
    t = ps.tokens
    a = a0_closed_form(gap_um, w0=t["POP_PARAM1"] * 1000, w2=t["POP_FPARAM1"] * 1000, lam=ps.cfgd["wavelength_um"])
    for cond, msg in ((abs(rep["pilot_waist_um"] / 1000 - t["POP_PARAM1"]) / t["POP_PARAM1"] < 1e-4, "source waist"),
                      (abs(rep["pilot_pos_um"] - gap_um) < 1e-4, "pilot position != gap"),
                      (abs(rep["pilot_size_um"] - a["w1"]) / a["w1"] < 1e-4, "beam radius at fiber != analytic"),
                      (t["POP_WIDEX"] * 1000 >= 4 * a["w1"], "window narrower than 4x beam radius")):
        if not cond:
            raise SettingsMismatch("report: " + msg)
    return dict(grid=rep["grid"], pilot_size_um=rep["pilot_size_um"], pilot_waist_um=rep["pilot_waist_um"],
                pilot_pos_um=rep["pilot_pos_um"])


# ------------------------------------------------------------------------- sweeps
def sweep_gap(ps, gap_um, step_um, range_um, log):
    ps.set_gap_um(gap_um)
    rep = a0_report_checks(ps, gap_um)
    ps.set("POP_FPARAM3", 0.0)
    a = a0_closed_form(gap_um)
    eta0 = ps.popd()
    if abs(eta0[0] - a["eta0"]) > 1e-4:
        raise SettingsMismatch("gap %g: on-axis eta %.6f vs analytic %.6f (receiver waist?)" % (gap_um, eta0[0], a["eta0"]))
    pts = []
    for d in np.round(np.arange(-range_um, range_um + step_um / 2, step_um), 6):
        ps.set("POP_FPARAM3", float(d) / 1000.0)
        eta, s, t = ps.popd()
        pts.append((float(d), eta, s, t, -10 * math.log10(eta)))
    ps.set("POP_FPARAM3", 0.0)
    arr = np.array(pts)
    # symmetry: +d and -d must agree
    asym = float(np.max(np.abs(arr[:, 1] - arr[::-1, 1])))
    if asym > 1e-6:
        raise RuntimeError("gap %g: eta(+d) != eta(-d), max diff %.2e" % (gap_um, asym))
    tol = tolerance(arr[:, 0], arr[:, 4], 0.0)
    # convergence (shared project test) at d = 0 and at the predicted tolerance
    res = C.convergence(ps, [(0.0, lambda: ps.set("POP_FPARAM3", 0.0)),
                             (round(a["d_exact"], 3), lambda: ps.set("POP_FPARAM3", round(a["d_exact"], 3) / 1000))],
                        after=lambda: ps.set("POP_FPARAM3", 0.0))
    conv = {c["point"]: max(c["change_pct"].values(), key=abs) for c in res}
    if not all(c["passed"] for c in res):
        raise RuntimeError("gap %g: convergence failed: %s" % (gap_um, res))
    return arr, tol, conv, rep, asym, a


def run(app, out=OUT, zos=None):
    zos = zos or load_zosapi()
    cfgd = json.loads(CONFIG.read_text())
    sw = cfgd["sweep_lateral"]
    log("extractor self-test (analytic data, gap 20 um): %.6f um" % selftest())
    out.mkdir(parents=True, exist_ok=True)
    ps = PopSession(app, cfgd, zos)
    nominal_loss = None
    rows, t0 = [], time.time()
    try:
        # reproducibility: the generated CFG must reproduce the committed Stage 3 value
        ps.set_gap_um(sw["nominal_gap_um"])
        ref = float(list(csv.DictReader(open(ROOT / "zemax" / "baseline" / "a0_baseline_popd.csv")))[0]["popd0_eta_total"])
        got = ps.popd()[0]
        if abs(got - ref) > 1e-6:
            raise SettingsMismatch("generated CFG gives eta %.9f, Stage 3 committed %.9f" % (got, ref))
        log("generated CFG reproduces Stage 3 eta = %.9f" % got)
        for gap in sw["gaps_um"]:
            step = sw["nominal_step_um"] if gap == sw["nominal_gap_um"] else sw["step_um"]
            arr, tol, conv, rep, asym, a = sweep_gap(ps, gap, step, sw["range_um"], log)
            with open(out / ("a0_lateral_gap%03dum.csv" % gap), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["x_decenter_um", "popd0_eta_total", "popd1_S_system", "popd2_T_receiver", "loss_dB"])
                for r in arr:
                    w.writerow(["%.4f" % r[0]] + ["%.9f" % v for v in r[1:4]] + ["%.6f" % r[4]])
            rows.append(dict(gap_um=gap, step_um=step, eta0=arr[len(arr) // 2, 1], L0_dB=tol["L0"],
                             d1dB_rel_um=tol["rel"], d1dB_rel_minus_um=tol["rel_minus"], d1dB_rel_plus_um=tol["rel_plus"],
                             d1dB_abs_um=tol["abs"], pred_exact_um=a["d_exact"], pred_simple_um=a["d_simple"],
                             pred_abs_um=a["d_abs"], beam_radius_um=a["w1"],
                             diff_vs_exact_pct=100 * (tol["rel"] - a["d_exact"]) / a["d_exact"],
                             diff_vs_simple_pct=100 * (tol["rel"] - a["d_simple"]) / a["d_simple"],
                             conv_pct_d0=conv[0.0], conv_pct_dtol=[v for k, v in conv.items() if k != 0.0][0],
                             asym=asym))
            log("gap %3g um: d1dB = %.4f um (exact %.4f, simple %.4f)  S min %.6f  conv %.4f%%/%.4f%%  [%.0fs]" % (
                gap, tol["rel"], a["d_exact"], a["d_simple"], arr[:, 2].min(), rows[-1]["conv_pct_d0"],
                rows[-1]["conv_pct_dtol"], time.time() - t0))
        # interpolation resolution check at the nominal gap: thin the 0.05 um data to 0.2 um
        nom = np.loadtxt(out / ("a0_lateral_gap%03dum.csv" % sw["nominal_gap_um"]), delimiter=",", skiprows=1)
        thin = nom[::4] if abs(nom[::4][len(nom[::4]) // 2, 0]) < 1e-9 else nom[::2]
        res_check = tolerance(thin[:, 0], thin[:, 4], 0.0)["rel"] - tolerance(nom[:, 0], nom[:, 4], 0.0)["rel"]
        log("resolution check (0.05 um vs thinned step): change %.2e um" % res_check)
    finally:
        ps.set_gap_um(sw["nominal_gap_um"])
        ps.set("POP_FPARAM3", 0.0)
        ps.commit()
        ps.close()
    keys = list(rows[0].keys())
    with open(out / "a0_tolerance_vs_gap.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (("%.6f" % v) if isinstance(v, float) else v) for k, v in r.items()})
    plot_tolerance(rows, out / "a0_tolerance_vs_gap.png")
    meta = dict(stage="5", description="A0 lateral X decenter of the receiver", date=time.strftime("%Y-%m-%d"),
                opticstudio=dict(build=str(app.OpticStudioVersion), license=str(app.LicenseStatus), mode=str(app.Mode)),
                source_config=str(CONFIG.relative_to(ROOT)), model=str(MODEL.relative_to(ROOT)),
                pop_tokens=dict(ps.tokens),
                sweep=sw, convergence_rule="double grid, <0.5% change in eta, checked at d=0 and d=tolerance for every gap",
                settings_readbacks_passed=ps.nchecks, resample_after_refraction="False on all surfaces (no refracting surface in A0)", surface_settings=cfgd.get("surface_settings"),
                tolerance_method="cubic spline of loss(d), brentq root; both sides, mean reported",
                tolerance_conventions={"rel": "loss rises 1 dB above on-axis value (primary)", "abs": "total loss reaches 1.000 dB"},
                interpolation_resolution_check_um=res_check, wall_time_s=time.time() - t0)
    (out / "run_config_a0_lateral.json").write_text(json.dumps(meta, indent=2, default=str))
    return rows


# --------------------------------------------------------------------------- plot
def plot_tolerance(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3de"
    blue, orange = "#2a78d6", "#eb6834"
    g = np.linspace(0, 105, 300)
    ex = [a0_closed_form(x)["d_exact"] for x in g]
    si = [a0_closed_form(x)["d_simple"] for x in g]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    ax.plot(g, ex, color=orange, lw=1.8, label="Exact Gaussian overlap (with wavefront curvature)")
    ax.plot(g, si, color=ink2, lw=1.4, ls="--", label="0.339 sqrt(w1^2 + w2^2)")
    ax.plot([r["gap_um"] for r in rows], [r["d1dB_rel_um"] for r in rows], color=blue, lw=0, marker="o", ms=7,
            mec="white", mew=1.5, label="Zemax POP, interpolated", zorder=5)
    for r in rows:
        ax.annotate("%.2f" % r["d1dB_rel_um"], (r["gap_um"], r["d1dB_rel_um"]), textcoords="offset points",
                    xytext=(0, -15), ha="center", fontsize=8, color=ink2)
    ax.set_xlabel("Air gap, grating plane to fiber facet (um)", color=ink)
    ax.set_ylabel("Lateral 1-dB tolerance, one side (um)", color=ink)
    ax.set_title("A0 lateral 1-dB tolerance versus air gap", loc="left", fontsize=11, color=ink)
    ax.set_xlim(-3, 108)
    ax.set_ylim(2.0, None)
    ax.grid(color=grid, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(ink2)
    ax.tick_params(colors=ink2)
    ax.legend(frameon=False, loc="upper left", fontsize=8.5, labelcolor=ink)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================== Stage 7
STAGE7 = [ROOT / "zemax" / "microlens" / "A1" / "run_config.json",
          ROOT / "zemax" / "microlens" / "B" / "run_config.json"]
MODEL_KEYS = ("wavelength_um", "catalogs", "surfaces", "gap_surface", "parameters", "pop_tokens", "surface_settings")
INK, INK2, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834"


def one_cell_difference(a, b):
    """A1 and B must differ in exactly one model cell: the exit-surface radius."""
    diffs = []
    for k in MODEL_KEYS:
        if k == "surfaces":
            for i, (sa, sb) in enumerate(zip(a[k], b[k])):
                diffs += ["surfaces[%d].%s" % (i, f) for f in sorted(set(sa) | set(sb)) if sa.get(f) != sb.get(f)]
        elif a.get(k) != b.get(k):
            diffs.append(k)
    if diffs != ["surfaces[2].radius_mm"]:
        raise SettingsMismatch("A1 and B must differ only in surfaces[2].radius_mm; they differ in %s" % diffs)
    return diffs


def predicted(ps, gap_um, n_after, wf_um=None):
    """Independent ABCD + numerical-overlap prediction (w, R, eta) at this gap."""
    c = ps.cfgd
    w, R = AN.beam_at(c["surfaces"], gap_um, n_after, c["parameters"]["source_waist_um"], c["wavelength_um"], ps.gap_surface)
    return w, R, AN.overlap(w, R, wf_um or ps.tokens["POP_FPARAM1"] * 1000, c["wavelength_um"])


def sweep(ps, kind, values):
    """Sweep receiver X decenter (um), receiver tilt (deg) or gap (um). Rows: value, eta, S, T, loss."""
    rows = []
    for v in values:
        if kind == "lateral":
            ps.set("POP_FPARAM3", v / 1000)
        elif kind == "angular":
            ps.set("POP_TILTX", v)
        else:
            ps.set_gap_um(v)
        eta, s_, t_ = ps.popd()
        rows.append((v, eta, s_, t_, float(loss_db(eta))))
    ps.set("POP_FPARAM3", 0.0)
    ps.set("POP_TILTX", 0.0)
    return rows


def sign_test(ps, n_after):
    """Pilot beam radius at the fibre for both signs of the exit-surface radius."""
    s2 = ps.sys.LDE.GetSurfaceAt(2)
    nominal = s2.Radius
    out = {}
    for r in (-abs(nominal), abs(nominal)):
        s2.Radius = r
        rep = ps.report(strict=False)          # diagnostic: record what POP did, do not abort
        surf = [dict(x) for x in ps.cfgd["surfaces"]]
        surf[2]["radius_mm"] = r
        w, R = AN.beam_at(surf, ps.cfgd["parameters"]["gap_um"], n_after, ps.cfgd["parameters"]["source_waist_um"],
                         ps.cfgd["wavelength_um"], ps.gap_surface)
        out["R = %+.0f um" % (r * 1000)] = dict(pop_pilot_radius_um=rep["pilot_size_um"], abcd_radius_um=w,
                                               reported_window_mm=rep["width_mm"],
                                               abcd_wavefront_R_um=R, eta=ps.popd()[0])
    s2.Radius = nominal
    C.check_model(ps.sys, ps.cfgd)
    return out


def propagator_test(ps):
    """Nominal POPD and the window POP actually used, with each propagator."""
    out = {}
    keys = [k for k in ps.cfgd["surface_settings"] if not k.startswith("_")]
    for asp in (True, False):
        for k in keys:
            ps.sys.LDE.GetSurfaceAt(int(k)).PhysicalOpticsData.UseAngularSpectrumPropagator = asp
        rep = ps.report(strict=False)
        out["angular_spectrum" if asp else "default"] = dict(popd=ps.popd(), reported_window_mm=rep["width_mm"],
                                                              point_spacing_um=rep["width_mm"] * 1000 / rep["grid"])
    ps.apply_surface_settings(ps.cfgd["surface_settings"])
    return out


def convergence(ps, points):
    """Shared project convergence test; restores the nominal state after each point."""
    def after():
        ps.set("POP_FPARAM3", 0.0)
        ps.set("POP_TILTX", 0.0)
        ps.set_gap_um(ps.cfgd["parameters"]["gap_um"])
    return C.convergence(ps, points, after=after)


def run_stage7(app, zos, build=True):
    """Build A1 and B from their run_config.json files and measure them."""
    cfgs = [json.loads(p.read_text()) for p in STAGE7]
    one_cell_difference(*cfgs)
    OUT.mkdir(parents=True, exist_ok=True)
    meta, curves, t0 = {}, {}, time.time()
    for cfgd in cfgs:
        cid, gap0, lam = cfgd["id"], cfgd["parameters"]["gap_um"], cfgd["wavelength_um"]
        if build:
            C.build_model(app, cfgd, zos)
        ps = PopSession(app, cfgd, zos)
        try:
            n_after = {i: ps.index_after(i) for i in range(len(cfgd["surfaces"]) - 1)}
            rep = ps.report()
            w_p, R_p, eta_p = predicted(ps, gap0, n_after)
            if abs(rep["pilot_size_um"] - w_p) / w_p > 1e-3:
                raise SettingsMismatch("%s: POP pilot radius %.4f um != ABCD %.4f um" % (cid, rep["pilot_size_um"], w_p))
            nominal = ps.popd()
            m = dict(n_after=n_after, nominal_popd=nominal, nominal_loss_dB=float(loss_db(nominal[0])),
                     predicted=dict(w_um=w_p, R_um=R_p, eta=eta_p, loss_dB=float(loss_db(eta_p))),
                     pilot={k: rep[k] for k in ("pilot_size_um", "pilot_waist_um", "pilot_pos_um", "width_mm", "grid")},
                     pop_messages=rep["messages"])
            log("%s nominal: eta %.6f (S %.6f, T %.6f) = %.3f dB; predicted %.3f dB | w %.3f um (ABCD %.3f), R %.1f um" % (
                cid, nominal[0], nominal[1], nominal[2], loss_db(nominal[0]), loss_db(eta_p), rep["pilot_size_um"], w_p, R_p))
            if cid == "B":
                m["sign_test"] = sign_test(ps, n_after)
                log("sign test: " + json.dumps(m["sign_test"]))
            m["propagators"] = propagator_test(ps)
            log("propagators: " + json.dumps(m["propagators"]))
            gmax = cfgd["sweeps"]["gap_um"][1]
            m["convergence"] = convergence(ps, [
                ("nominal", lambda: None),
                ("lateral 8 um", lambda: ps.set("POP_FPARAM3", 0.008)),
                ("tilt 0.6 deg", lambda: ps.set("POP_TILTX", 0.6)),
                ("gap %g um (largest beam)" % gmax, lambda: ps.set_gap_um(gmax))])
            for c in m["convergence"]:
                log("convergence %-24s %s" % (c["point"], {k: "%+.4f%%" % v for k, v in c["change_pct"].items()}))
                if any(abs(v) > 0.5 for v in c["change_pct"].values()):
                    raise RuntimeError("%s: convergence failed at %s" % (cid, c["point"]))
            # receiver variants at nominal: SMF control, and the best fibre for this beam
            wbest = float(max(np.arange(2.0, 40.0, 0.1), key=lambda x: AN.overlap(w_p, R_p, x, lam)))
            m["receiver_variants"] = {}
            for label, wf in (("SMF, 9.2 um MFD", 4.6), ("best fibre, %.1f um radius" % wbest, wbest)):
                for k in ("POP_FPARAM1", "POP_FPARAM2"):
                    ps.set(k, round(wf, 4) / 1000)
                e = ps.popd()
                v = dict(waist_um=wf, popd=e, loss_dB=float(loss_db(e[0])), predicted_eta=AN.overlap(w_p, R_p, wf, lam))
                if cid == "B" and wf == 4.6:
                    rows = sweep(ps, "lateral", frange(-24, 24, 0.25))
                    write_csv(OUT / "B_SMF_lateral.csv", ["offset_um", "popd0_eta_total", "popd1_S_system",
                                                          "popd2_T_receiver", "loss_dB"], rows)
                    arr = np.array(rows)
                    v["lateral"] = AN.tolerance(arr[:, 0], arr[:, 4], 0.0)
                m["receiver_variants"][label] = v
                for k in ("POP_FPARAM1", "POP_FPARAM2"):
                    ps.set(k, cfgd["pop_tokens"][k])
            log("receiver variants: " + json.dumps({k: round(v["loss_dB"], 3) for k, v in m["receiver_variants"].items()}))
            # longitudinal first: it fixes the re-optimised gap
            sw = cfgd["sweeps"]
            gaps = frange(*sw["gap_um"])
            rows = sweep(ps, "longitudinal", gaps)
            ps.set_gap_um(gap0)
            arr = np.array(rows)
            pred = np.array([predicted(ps, g, n_after)[2] for g in gaps])
            zopt, lmin = AN.peak(arr[:, 0], arr[:, 4])
            res = dict(z_opt_um=zopt, L_min_dB=lmin,
                       longitudinal_nominal=AN.tolerance(arr[:, 0], arr[:, 4], gap0),
                       longitudinal_zopt=AN.tolerance(arr[:, 0], arr[:, 4], zopt),
                       longitudinal_nominal_predicted=AN.tolerance(gaps, loss_db(pred), gap0),
                       longitudinal_zopt_predicted=AN.tolerance(gaps, loss_db(pred), AN.peak(gaps, loss_db(pred))[0]))
            write_csv(OUT / ("%s_longitudinal.csv" % cid), ["gap_um", "popd0_eta_total", "popd1_S_system",
                      "popd2_T_receiver", "loss_dB", "predicted_loss_dB"],
                      [r + (float(loss_db(p)),) for r, p in zip(rows, pred)])
            curves[(cid, "longitudinal")] = (arr[:, 0], arr[:, 4], loss_db(pred))
            log("%s longitudinal: z_opt %.2f um (L %.3f dB) | from nominal rel %.1f / abs %.1f um | from z_opt rel %.1f um" % (
                cid, zopt, lmin, res["longitudinal_nominal"]["rel_plus"], res["longitudinal_nominal"]["abs_plus"],
                res["longitudinal_zopt"]["rel_plus"]))
            # lateral and angular, at the nominal gap and at the re-optimised gap
            for zlabel, z in (("nominal", gap0), ("zopt", max(zopt, 0.0))):
                ps.set_gap_um(z)
                w_z, R_z, _ = predicted(ps, z, n_after)
                wf = ps.tokens["POP_FPARAM1"] * 1000
                for kind, rng in (("lateral", sw["lateral_um"]), ("angular", sw["angular_deg"])):
                    vals = frange(*rng)
                    rows = sweep(ps, kind, vals)
                    arr = np.array(rows)
                    if np.max(np.abs(arr[:, 1] - arr[::-1, 1])) > 1e-6:
                        raise RuntimeError("%s %s: sweep not symmetric" % (cid, kind))
                    kw = (lambda v: dict(d=v)) if kind == "lateral" else (lambda v: dict(theta_deg=v))
                    pred = np.array([AN.overlap(w_z, R_z, wf, lam, **kw(v)) for v in vals])
                    key = "%s_%s" % (kind, zlabel)
                    res[key] = AN.tolerance(arr[:, 0], arr[:, 4], 0.0)
                    res[key + "_predicted"] = AN.tolerance(vals, loss_db(pred), 0.0)
                    write_csv(OUT / ("%s_%s%s.csv" % (cid, kind, "" if zlabel == "nominal" else "_zopt")),
                              ["offset_um" if kind == "lateral" else "tilt_deg", "popd0_eta_total", "popd1_S_system",
                               "popd2_T_receiver", "loss_dB", "predicted_loss_dB"],
                              [r + (float(loss_db(p)),) for r, p in zip(rows, pred)])
                    if zlabel == "nominal":
                        curves[(cid, kind)] = (arr[:, 0], arr[:, 4], loss_db(pred))
                    log("%s %-8s @%-7s (gap %.2f um): rel %.4f  abs %.4f | predicted rel %.4f  [%.0fs]" % (
                        cid, kind, zlabel, z, res[key]["rel"], res[key]["abs"], res[key + "_predicted"]["rel"], time.time() - t0))
            ps.set_gap_um(gap0)
            m["tolerances"] = res
            m["settings_readbacks_passed"] = ps.nchecks
            meta[cid] = m
            cfg_path = ROOT / "zemax" / "microlens" / cid / "run_config.json"
            c2 = json.loads(cfg_path.read_text())
            c2["opticstudio"] = dict(version="2026 R1.00", build=str(app.OpticStudioVersion),
                                     license=str(app.LicenseStatus), mode=str(app.Mode))
            c2["results_nominal"] = dict(popd0_eta_total=nominal[0], popd1_S_system=nominal[1],
                                         popd2_T_receiver=nominal[2], loss_dB=float(loss_db(nominal[0])),
                                         pilot=m["pilot"], index_after_surface=n_after)
            cfg_path.write_text(json.dumps(c2, indent=2))
        finally:
            ps.close()
    (OUT / "run_config_stage7.json").write_text(json.dumps(dict(
        stage="7", configs=[str(p.relative_to(ROOT)) for p in STAGE7], date=time.strftime("%Y-%m-%d"),
        opticstudio=dict(build=str(app.OpticStudioVersion), license=str(app.LicenseStatus), mode=str(app.Mode)),
        tolerance_method="cubic spline of loss, brentq root. rel: 1 dB above the reference point; abs: total loss 1.000 dB",
        results=meta, wall_time_s=time.time() - t0), indent=2, default=str))
    plot_stage7(curves, OUT / "stage7_A1_vs_B.png")
    return meta


def plot_stage7(curves, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0), dpi=150, sharey=True)
    labels = dict(lateral="Lateral offset (um)", angular="Fibre tilt (deg)", longitudinal="Air gap (um)")
    for ax, kind in zip(axes, ("lateral", "angular", "longitudinal")):
        for cid, col in (("A1", ORANGE), ("B", BLUE)):
            x, y, yp = curves[(cid, kind)]
            ax.plot(x, yp, color=INK2, lw=1, ls="--", label="analytic (ABCD + overlap)" if cid == "A1" else None)
            ax.plot(x, y, color=col, lw=2, label="%s, Zemax POP" % cid)
        ax.set_xlabel(labels[kind], color=INK)
        ax.set_title(kind.capitalize(), loc="left", fontsize=10, color=INK)
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2)
    axes[0].set_ylabel("Coupling loss (dB)", color=INK)
    axes[0].set_ylim(0, 12)
    axes[0].legend(frameon=False, fontsize=8, loc="upper center")
    fig.suptitle("A1 (flat exit) vs B (micro-lens), 34 um MFD TEC receiver, nominal 20 um gap",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================== Stage 8
B592_CONFIG = ROOT / "zemax" / "microlens" / "B592" / "run_config.json"


def run_stage8(app, zos):
    """Fill the data gaps: A0 angular and longitudinal; B rebuilt with the 2021-paper pair (B592)."""
    OUT.mkdir(parents=True, exist_ok=True)
    out, t0 = {}, time.time()
    # ---- A0 angular and longitudinal ---------------------------------------------------
    cfgd = json.loads(CONFIG.read_text())
    ps = PopSession(app, cfgd, zos)
    try:
        gap0, lam = cfgd["sweep_lateral"]["nominal_gap_um"], cfgd["wavelength_um"]
        ps.set_gap_um(gap0)
        if abs(ps.popd()[0] - 0.962617754) > 1e-6:
            raise SettingsMismatch("A0 does not reproduce the Stage 3 nominal eta")
        _a = a0_closed_form(gap0)
        w1, R1 = _a["w1"], _a["R1"]
        th = frange(-5.0, 5.0, 0.02)
        rows = sweep(ps, "angular", th)
        arr = np.array(rows)
        if np.max(np.abs(arr[:, 1] - arr[::-1, 1])) > 1e-6:
            raise RuntimeError("A0 angular sweep not symmetric")
        pred = np.array([AN.overlap(w1, R1, 4.6, lam, theta_deg=v) for v in th])
        write_csv(OUT / "A0_angular.csv", ["tilt_deg", "popd0_eta_total", "popd1_S_system", "popd2_T_receiver",
                                           "loss_dB", "predicted_loss_dB"],
                  [r + (float(loss_db(p)),) for r, p in zip(rows, pred)])
        ang = AN.tolerance(arr[:, 0], arr[:, 4], 0.0)
        ang_p = AN.tolerance(th, loss_db(pred), 0.0)
        gaps = frange(0.0, 150.0, 0.5)
        rows = sweep(ps, "longitudinal", gaps)
        ps.set_gap_um(gap0)
        arr = np.array(rows)
        pred = np.array([1 / (1 + (g / (2 * math.pi * 4.6 ** 2 / lam)) ** 2) for g in gaps])
        write_csv(OUT / "A0_longitudinal.csv", ["gap_um", "popd0_eta_total", "popd1_S_system", "popd2_T_receiver",
                                                "loss_dB", "predicted_loss_dB"],
                  [r + (float(loss_db(p)),) for r, p in zip(rows, pred)])
        zopt, lmin = AN.peak(arr[:, 0], arr[:, 4])
        conv = convergence_a0(ps, gap0)
        out["A0"] = dict(angular_nominal=ang, angular_nominal_predicted=ang_p,
                         longitudinal_nominal=AN.tolerance(arr[:, 0], arr[:, 4], gap0),
                         longitudinal_zopt=AN.tolerance(arr[:, 0], arr[:, 4], zopt), z_opt_um=zopt, L_min_dB=lmin,
                         longitudinal_nominal_predicted=AN.tolerance(gaps, loss_db(pred), gap0),
                         longitudinal_zopt_predicted=AN.tolerance(gaps, loss_db(pred), 0.0),
                         convergence=conv, settings_readbacks_passed=ps.nchecks)
        log("A0 angular: rel %.4f abs %.4f deg (predicted %.4f / %.4f) [%.0fs]" % (
            ang["rel"], ang["abs"], ang_p["rel"], ang_p["abs"], time.time() - t0))
        log("A0 longitudinal: z_opt %.2f um; from 20 um rel +%.2f abs +%.2f; from z_opt rel +%.2f um" % (
            zopt, out["A0"]["longitudinal_nominal"]["rel_plus"], out["A0"]["longitudinal_nominal"]["abs_plus"],
            out["A0"]["longitudinal_zopt"]["rel_plus"]))
        log("A0 convergence: %s" % conv)
    finally:
        ps.set_gap_um(cfgd["sweep_lateral"]["nominal_gap_um"])
        ps.set("POP_TILTX", 0.0)
        ps.commit()
        ps.close()
    # ---- B592: the 2021-paper literature pair -------------------------------------------
    cfgd = json.loads(B592_CONFIG.read_text())
    C.build_model(app, cfgd, zos)
    ps = PopSession(app, cfgd, zos)
    try:
        lam, gap0 = cfgd["wavelength_um"], cfgd["parameters"]["gap_um"]
        n_after = {i: ps.index_after(i) for i in range(len(cfgd["surfaces"]) - 1)}
        rep = ps.report()
        w_p, R_p, eta_p = predicted(ps, gap0, n_after)
        if abs(rep["pilot_size_um"] - w_p) / w_p > 1e-3:
            raise SettingsMismatch("B592: POP pilot radius %.4f um != ABCD %.4f um" % (rep["pilot_size_um"], w_p))
        nominal = ps.popd()
        conv = convergence(ps, [("nominal", lambda: None), ("lateral 8 um", lambda: ps.set("POP_FPARAM3", 0.008)),
                                ("tilt 0.6 deg", lambda: ps.set("POP_TILTX", 0.6))])
        if any(abs(v) > 0.5 for c in conv for v in c["change_pct"].values()):
            raise RuntimeError("B592 convergence failed: %s" % conv)
        res = dict(nominal_popd=nominal, nominal_loss_dB=float(loss_db(nominal[0])),
                   pilot_radius_um=rep["pilot_size_um"], abcd=dict(w_um=w_p, R_um=R_p, eta=eta_p),
                   convergence=conv, n_after=n_after)
        wf = ps.tokens["POP_FPARAM1"] * 1000
        for kind, rng in (("lateral", cfgd["sweeps"]["lateral_um"]), ("angular", cfgd["sweeps"]["angular_deg"])):
            vals = frange(*rng)
            rows = sweep(ps, kind, vals)
            arr = np.array(rows)
            kw = (lambda v: dict(d=v)) if kind == "lateral" else (lambda v: dict(theta_deg=v))
            pred = np.array([AN.overlap(w_p, R_p, wf, lam, **kw(v)) for v in vals])
            res[kind] = AN.tolerance(arr[:, 0], arr[:, 4], 0.0)
            res[kind + "_predicted"] = AN.tolerance(vals, loss_db(pred), 0.0)
            write_csv(OUT / ("B592_%s.csv" % kind), ["offset_um" if kind == "lateral" else "tilt_deg", "popd0_eta_total",
                      "popd1_S_system", "popd2_T_receiver", "loss_dB", "predicted_loss_dB"],
                      [r + (float(loss_db(p)),) for r, p in zip(rows, pred)])
        res["settings_readbacks_passed"] = ps.nchecks
        out["B592"] = res
        log("B592 nominal %.4f dB, w %.3f um (ABCD %.3f); lateral rel %.4f abs %.4f um; angular rel %.4f deg [%.0fs]" % (
            res["nominal_loss_dB"], rep["pilot_size_um"], w_p, res["lateral"]["rel"], res["lateral"]["abs"],
            res["angular"]["rel"], time.time() - t0))
        c2 = dict(cfgd)
        c2["opticstudio"] = dict(version="2026 R1.00", build=str(app.OpticStudioVersion),
                                 license=str(app.LicenseStatus), mode=str(app.Mode))
        c2["results_nominal"] = dict(popd0_eta_total=nominal[0], popd1_S_system=nominal[1], popd2_T_receiver=nominal[2],
                                     loss_dB=res["nominal_loss_dB"], pilot_radius_um=rep["pilot_size_um"],
                                     index_after_surface=n_after)
        B592_CONFIG.write_text(json.dumps(c2, indent=2))
    finally:
        ps.close()
    (OUT / "run_config_stage8.json").write_text(json.dumps(dict(
        stage="8", date=time.strftime("%Y-%m-%d"),
        opticstudio=dict(build=str(app.OpticStudioVersion), license=str(app.LicenseStatus), mode=str(app.Mode)),
        sweeps=dict(A0_angular_deg=[-5, 5, 0.02], A0_gap_um=[0, 150, 0.5]), results=out,
        wall_time_s=time.time() - t0), indent=2, default=str))
    return out


def convergence_a0(ps, gap0):
    """Shared convergence test at the A0 angular tolerance point and at the largest gap."""
    def after():
        ps.set("POP_TILTX", 0.0)
        ps.set_gap_um(gap0)
    res = C.convergence(ps, [("tilt 2.4 deg", lambda: ps.set("POP_TILTX", 2.4)),
                             ("gap 150 um", lambda: ps.set_gap_um(150.0))], after=after)
    if not all(c["passed"] for c in res):
        raise RuntimeError("A0 convergence failed: %s" % res)
    return {c["point"]: max(c["change_pct"].values(), key=abs) for c in res}


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("5", "7", "8"), default="8")
    ap.add_argument("--mode", choices=("standalone", "extension"), default="standalone")
    args = ap.parse_args()
    zos = load_zosapi()
    conn, app = C.connect(args.mode, zos)
    try:
        {"5": lambda: run(app, zos=zos), "7": lambda: run_stage7(app, zos), "8": lambda: run_stage8(app, zos)}[args.stage]()
    finally:
        if args.mode == "standalone":
            app.CloseApplication()


if __name__ == "__main__":
    main()
